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
    delete [] this->res_hist;
    if(this->Gspace) delete [] c_fm;
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
    size_t pbasis = fine_pwaves->pbasis;

    this->step = this->step % this->refresh_steps;
    int N = int(this->Nsize);
    
    int lda = this->max_order +1;
    double *A_mat = this->A_mat;

    if(this->need_precond) this->Precond(fm, this->nstates);
    if(this->Gspace)
    {
        double tpiba = 2.0 * PI / Rmg_L.celldm[0];
        double tpiba2 = tpiba * tpiba;
        for(int is = 0; is < ct.noncoll_factor * ct.noncoll_factor; is++)
        {
            for(size_t ig=0; ig < pbasis; ig++) c_fm[ig] = std::complex<double>(fm[is*pbasis + ig], 0.0);
            fine_pwaves->FftForward(c_fm, c_fm);

            for(size_t ig=0;ig < pbasis;ig++) {
                // fm is the charge density difference, sums of them always equal to zero
                // so there is no G= 0 term.
                double g2 = fine_pwaves->gmags[ig] * tpiba2;

                double alpha = std::max(0.0, g2/(g2+this->ktf * this->ktf));
                c_fm[ig] = c_fm[ig] * alpha;
                // if(!fine_pwaves->gmask[ig]) c_fm[ig]  = 0.0;
                //c_fm[ig] = c_fm[ig] * g2/(g2+this->ktf * this->ktf);
            }
            fine_pwaves->FftInverse(c_fm, c_fm);
            for(size_t i = 0;i < pbasis;i++) fm[is*pbasis + i] = std::real(c_fm[i])/(double)fine_pwaves->global_basis;
        }
    }
    if(this->pulay_order <=1)
    {

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

void PulayMixing::SetGspace()
{
    this->Gspace = true;
    this->c_fm = new std::complex<double>[Nsize];
}



