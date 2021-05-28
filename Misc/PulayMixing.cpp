#include <iostream>     // std::cout
#include <algorithm>    // std::rotate
#include <vector>       // std::vector
#include <math.h>       // std::vector
#include <mpi.h>       // std::vector
#include "PulayMixing.h"
#include "blas.h"
#include "GlobalSums.h"
#include "transition.h"
#include "RmgParallelFft.h"
#include "RmgException.h"

PulayMixing::PulayMixing(size_t Nsize, int pulay_order, int refresh_steps, double mix_first, 
        double beta, MPI_Comm comm)
{
    this->Nsize = Nsize;
    this->pulay_order = pulay_order;

    if(pulay_order > this->max_order)
    {
        std::cout << " pulay_order is too large " << pulay_order << std::endl;
        exit(0);
    }

    this->refresh_steps = refresh_steps;
    this->mix_first = mix_first;
    this->beta = beta;
    this->comm = comm;
    this->hist = new double[Nsize * (size_t)(pulay_order) +1024];
    this->res_hist = new double[Nsize * (size_t)(pulay_order) + 1024];
    this->A_mat = new double[(this->max_order+1)*(this->max_order+1)];

    for(int i = 0; i < this->pulay_order;i++)
    {
        this->hist_ptr.push_back(&this->hist[Nsize * i]);
        this->res_hist_ptr.push_back(&this->res_hist[Nsize * (size_t)i]);
    }

    this->step = 0;

    this->need_precond = 0;
}

PulayMixing::~PulayMixing(void)
{
    delete [] this->hist;
    if(this->Gspace)
    {
        delete [] this->res_histG;
    }
    else
    {
        delete [] this->res_hist;
    }
    if(c_fm != NULL) delete [] c_fm;
}

void PulayMixing::SetPrecond(std::function<void(double*, int)> precond)
{ 
    this->need_precond = 1;
    this->Precond = precond;
}

void PulayMixing::Refresh(){ this->step = 0;}

void PulayMixing::Mixing(double *xm, double *fm)
{
    double A[(this->max_order+1) * (this->max_order+1)];
    double b[this->max_order+1];
    int ipvt[this->max_order+1];
    int ione = 1;
    int info;

    this->step = this->step % this->refresh_steps;
    int N = int(this->Nsize);

    int lda = this->max_order +1;
    double *A_mat = this->A_mat;

    if(this->need_precond) this->Precond(fm, this->nstates);
    if(this->pulay_order <=1)
    {

        if(this->drho_pre)
        {
            Precond_drho(fm);
        }
        daxpy(&N, &this->mix_first, fm, &ione, xm, &ione);

        return;
    }

    // copy the xm and fm to the last history pointer.
    int current_pos = std::min(this->step, this->pulay_order-1);
    dcopy(&N, xm, &ione, this->hist_ptr[current_pos], &ione);
    dcopy(&N, fm, &ione, this->res_hist_ptr[current_pos], &ione);
    if (this->step == 0)
    {
        A_mat[this->step * lda + this->step] = ddot(&N, fm, &ione, fm, &ione);
        //       if(this->need_precond) this->Precond(fm, this->nstates);

        if(this->drho_pre)
        {
            Precond_drho(fm);
        }
        daxpy(&N, &this->mix_first, fm, &ione, xm, &ione);

        this->step++;
        return;
    }

    //  remove the first row and column of A matrix and other matrix elements  will be used for next step
    //  A_mat = <fi|fj> for i, j being residul from previous steps
    if(this->step >= this->pulay_order)
    {
        for(int i = 0; i < this->pulay_order-1; i++)
            for(int j = 0; j < this->pulay_order-1; j++)
            {
                A_mat[ i * lda +j] = A_mat[(i+1)*lda +j+1];
            }
    }

    //  calculate the <fi|fm> 
    int num_prev_steps = std::min(this->step, this->pulay_order-1);
    for(int i = 0; i < num_prev_steps; i++)
    {
        double *fi = this->res_hist_ptr[i];
        A_mat[i * lda + num_prev_steps] = ddot(&N, fi, &ione, fm, &ione);

        A_mat[num_prev_steps * lda + i] = 
            A_mat[i * lda + num_prev_steps] ;
    }

    //  calculate <fm|fm>
    A_mat[num_prev_steps * lda + num_prev_steps] = ddot(&N, fm, &ione, fm, &ione);

    int s2 = (this->max_order+1)*(this->max_order+1);

    dcopy(&s2, A_mat, &ione, A, &ione);
    GlobalSums(A, s2, comm);
    GlobalSums(A, s2, pct.spin_comm);

    int size = num_prev_steps + 1; 
    int A_size = size +1;

    for (int i = 0; i < size; i++)
    {
        A[i * lda + size] = 1.0;
        A[size * lda + i] = 1.0;
        b[i] = 0.0;
    }
    b[size] = 1.0;
    A[size*lda + size] = 0.0;

    /*   b = A^(-1) * b     */
    dgesv(&A_size, &ione, A, &lda, ipvt, b, &A_size, &info);

    if(pct.gridpe == 0 && 0)
    {
        printf("\n");
        for (int i = 0; i < size; i++)
            std::cout << "Pulay_b:" << i <<"  "<< b[i]<<std::endl;
        printf("\n");
    }

    dscal(&N, &b[size-1], xm, &ione);
    for (int i = 0; i < size - 1; i++)
    {
        daxpy(&N, &b[i], this->hist_ptr[i], &ione, xm, &ione);
    }

    dscal(&N, &b[size-1], fm, &ione);
    for (int i = 0; i < size - 1; i++)
    {
        daxpy(&N, &b[i], this->res_hist_ptr[i], &ione, fm, &ione);
    }

    if(this->drho_pre)
    {
        Precond_drho(fm);
    }
    daxpy(&N, &this->beta, fm, &ione, xm, &ione);

    // if the this->step larger than pulay_order, rotate the hist_ptr so that 
    //the first pointer becomes the last pointer which will be used for next step xm and fm.
    // otherwise the history pointers don't need to rotate.
    if (this->step >= this->pulay_order -1) 
    {
        std::rotate(this->hist_ptr.begin(),this->hist_ptr.begin()+1,this->hist_ptr.end());
        std::rotate(this->res_hist_ptr.begin(),this->res_hist_ptr.begin()+1,this->res_hist_ptr.end());
    }

    this->step++;

}

void PulayMixing::SetGspace(bool drho_pre_in, bool Gspace_in, double q0_in)
{

    this->c_fm = new std::complex<double>[Nsize];
    this->q0 = q0_in;
    if(drho_pre_in) 
    {
        this->drho_pre = true;
    }
    if(Gspace_in) {
        this->drho_pre = true;
        this->Gspace = true;
        delete [] this->res_hist;
        this->res_histG = new std::complex<double>[Nsize * (size_t)(pulay_order) + 1024];

        for(int i = 0; i < this->pulay_order;i++)
        {
            this->res_histG_ptr.push_back(&this->res_histG[Nsize * (size_t)i]);
        }
        // seet q1 for scalar product
        double qmin = 2.0 * PI / Rmg_L.celldm[0];
        qmin = std::min(qmin, qmin/Rmg_L.celldm[1]);
        qmin = std::min(qmin, qmin/Rmg_L.celldm[2]);
        double hx = Rmg_G->get_hxgrid(ct.FG_RATIO) * Rmg_L.celldm[0];
        double qmax = 2.0 * PI/hx;

        double weight_scale = 20.0;
        //  f_q = (q^2 + q1^2)/q^2
        //<A|B> = sum_q f_q * A_q * B_q
        // f_q(qmin) = weight_scale * f_q(qmax) 

        this->q1 = (weight_scale - 1.0) * qmin * qmin * qmax * qmax /(qmax * qmax - weight_scale * qmin *qmin);

        //std::cout << "q1 = " << this->q1 << std::endl;
        if(this->q1 < 0.0) {
            throw RmgFatalException() << "Error! q1^2 negative." << " in " << __FILE__ << " at line " << __LINE__ << "\n";

        }   
        
        this->q1 = std::sqrt(this->q1);
    }
}


void PulayMixing::Mixing_rhoG(double *xm, double *fm)
{

    //store the resuduals in G space 
    // The scalar products are calculated in Gspace and scaled by (q^2+q1^2)/q^2 
    // q1 is chosen so that qmin and qmax's weight differs by 20 times.
    // see SetGspace


    double tpiba = 2.0 * PI / Rmg_L.celldm[0];
    double tpiba2 = tpiba * tpiba;

    double A[(this->max_order+1) * (this->max_order+1)];
    double b[this->max_order+1];
    int ipvt[this->max_order+1];
    int ione = 1;
    int info;
    size_t pbasis = fine_pwaves->pbasis;

    this->step = this->step % this->refresh_steps;
    int N = int(this->Nsize);

    int lda = this->max_order +1;
    double *A_mat = this->A_mat;

    if(this->pulay_order <=1)
    {

        Precond_drho(fm);
        daxpy(&N, &this->mix_first, fm, &ione, xm, &ione);

        return;
    }

    // copy the xm and fm to the last history pointer.
    int current_pos = std::min(this->step, this->pulay_order-1);
    dcopy(&N, xm, &ione, this->hist_ptr[current_pos], &ione);

    int nspin = N/pbasis;
    for(int ig=0; ig < N; ig++) c_fm[ig] = std::complex<double>(fm[ig], 0.0);
    for(int is = 0; is < nspin; is++)
    {
        fine_pwaves->FftForward(&c_fm[is*pbasis], &c_fm[is*pbasis] );
    }

    zcopy(&N, c_fm, &ione, this->res_histG_ptr[current_pos], &ione);

    int num_prev_steps = std::min(this->step, this->pulay_order-1);



    //  remove the first row and column of A matrix and other matrix elements  will be used for next step
    //  A_mat = <fi|fj> for i, j being residul from previous steps
    if(this->step >= this->pulay_order)
    {
        for(int i = 0; i < this->pulay_order-1; i++)
            for(int j = 0; j < this->pulay_order-1; j++)
            {
                A_mat[ i * lda +j] = A_mat[(i+1)*lda +j+1];
            }
    }

    //  calculate the <fi|fm> 
    for(int i = 0; i < num_prev_steps; i++)
    {
        std::complex<double> *c_fi = this->res_histG_ptr[i];
        A_mat[i * lda + num_prev_steps] = 0.0;
        for(int idx = 0; idx < N; idx++) {
            int ig = idx%pbasis;
            double g2 = fine_pwaves->gmags[ig] * tpiba2;
            if( g2 < 1.0e-5) g2 =  tpiba2;

            double f_q = (g2 + q1*q1)/g2;

            A_mat[i * lda + num_prev_steps] += f_q * std::real(std::conj(c_fi[idx]) * c_fm[idx]);

        }

        A_mat[num_prev_steps * lda + i] = 
            A_mat[i * lda + num_prev_steps] ;
    }

    //  calculate <fm|fm>

    A_mat[num_prev_steps * lda + num_prev_steps] = 0.0;

    for(int idx = 0; idx < N; idx++) {
        int ig = idx%pbasis;
        double g2 = fine_pwaves->gmags[ig] * tpiba2;
        if( g2 < 1.0e-5) g2 =  tpiba2;

        double f_q = (g2 + q1*q1)/g2;

        A_mat[num_prev_steps * lda + num_prev_steps] += f_q * std::real(std::conj(c_fm[idx]) * c_fm[idx]);
    }

    int s2 = (this->max_order+1)*(this->max_order+1);

    dcopy(&s2, A_mat, &ione, A, &ione);
    GlobalSums(A, s2, comm);
    GlobalSums(A, s2, pct.spin_comm);

    int size = num_prev_steps + 1; 
    int A_size = size +1;

    for (int i = 0; i < size; i++)
    {
        A[i * lda + size] = 1.0;
        A[size * lda + i] = 1.0;
        b[i] = 0.0;
    }
    b[size] = 1.0;
    A[size*lda + size] = 0.0;

    /*   b = A^(-1) * b     */
    dgesv(&A_size, &ione, A, &lda, ipvt, b, &A_size, &info);

    if(pct.gridpe == 0 && ct.verbose)
    {
        printf("\n");
        for (int i = 0; i < size; i++)
            std::cout << "Pulay_b:" << i <<"  "<< b[i]<<std::endl;
        printf("\n");
    }

    dscal(&N, &b[size-1], xm, &ione);
    for (int i = 0; i < size - 1; i++)
    {
        daxpy(&N, &b[i], this->hist_ptr[i], &ione, xm, &ione);
    }

    std::complex<double> b_c = b[size-1];
    zscal(&N, &b_c, c_fm, &ione);
    for (int i = 0; i < size - 1; i++)
    {
        b_c = b[i];
        zaxpy(&N, &b_c, this->res_histG_ptr[i], &ione, c_fm, &ione);
    }

    for(int idx = 0; idx < N; idx++) {
        int ig = idx%pbasis;
        double g2 = fine_pwaves->gmags[ig] * tpiba2;
        if( g2 < 1.0e-5) g2 =  tpiba2;
        double alpha = g2/(g2+ q0 * q0);
        c_fm[idx] = c_fm[idx] * alpha;
    }

    for(int is = 0; is < nspin; is++)
    {
        fine_pwaves->FftInverse(&c_fm[is*pbasis], &c_fm[is*pbasis] );
    }

    for(int i = 0; i < N; i++) fm[i] = std::real(c_fm[i])/(double)fine_pwaves->global_basis;

    daxpy(&N, &this->beta, fm, &ione, xm, &ione);

    // if the this->step larger than pulay_order, rotate the hist_ptr so that 
    //the first pointer becomes the last pointer which will be used for next step xm and fm.
    // otherwise the history pointers don't need to rotate.
    if (this->step >= this->pulay_order -1) 
    {
        std::rotate(this->hist_ptr.begin(),this->hist_ptr.begin()+1,this->hist_ptr.end());
        std::rotate(this->res_histG_ptr.begin(),this->res_histG_ptr.begin()+1,this->res_histG_ptr.end());
    }

    this->step++;

}

