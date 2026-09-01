#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>

#include <complex>


#include "const.h"

#include "rmgtypedefs.h"
#include "typedefs.h"
#include "Kpoint.h"
#include "common_prototypes.h"
#include "common_prototypes1.h"
#include "transition.h"
#include "RmgException.h"
#include "AtomicInterpolate.h"


void InitWeight_xyz_One (SPECIES * sp, std::complex<double> *reptr, std::complex<double> *phaseptr, int ip, int l2mi, int l2mj);

/*This sets loop over species does forward fourier transform, finds and stores whatever is needed so that
 * only backwards Fourier transform is needed in the calculation*/
void SPECIES::InitLocalizedWeight (void)
{

    fftw_complex *betaptr;
    std::complex<double> *phaseptr;

    typedef struct {int species; int ip; int l; int m; int proj_index;} PROJ_INFO;
    PROJ_INFO proj;
    std::vector<PROJ_INFO> proj_iter;


    RmgTimer RT0("Weight");
    // get tot number of projectors and their information

    int nldim_max = ct.max_nldim;

    RmgTimer *RT1= new RmgTimer("Weight: phase and set");

    int size = this->nldim * this->nldim * this->nldim;
    this->phase = new fftw_complex[size * ct.num_kpts_pe];
    phaseptr = (std::complex<double> *)this->phase;
    // phase: multiply exp(-ik.r) in real space
    GetPhaseSpecies(this, phaseptr);
    /*Loop over all betas to calculate num of projectors for given species */
    int prjcount = 0;
    for (int ip = 0; ip < this->nbeta; ip++)
    {

        for(int m = 0; m < 2*this->llbeta[ip]+1; m++)
        {
            proj.ip = ip;
            proj.l = this->llbeta[ip];
            proj.m = m;
            proj.proj_index = prjcount;
            proj_iter.push_back(proj);
            prjcount++;
        }
    }

    size = this->nldim * this->nldim * this->nldim;
    this->num_projectors = prjcount;

    /*This array will store forward fourier transform on the coarse grid for all betas */
    if(this->forward_beta) fftw_free(this->forward_beta);
    this->forward_beta = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * this->num_projectors * size * ct.num_kpts_pe);

    if (this->forward_beta == NULL)
        throw RmgFatalException() << "cannot allocate mem "<< " at line " << __LINE__ << "\n";


    delete RT1;

    RmgTimer *RT3= new RmgTimer("Weight: proj cal");
    for(int iproj = pct.gridpe; iproj < this->num_projectors; iproj+=pct.grid_npes)
    {
        proj = proj_iter[iproj];

        // if this->nldim > get_NX_GRID, folding of neighbor cells are needed. 

        for(int kpt = 0; kpt <ct.num_kpts_pe; kpt++)
        {
            phaseptr = (std::complex<double> *) &this->phase[kpt * this->nldim * this->nldim * this->nldim];
            betaptr = &this->forward_beta[kpt *this->num_projectors *size + proj.proj_index * size];
            InitWeightOne(this, betaptr, phaseptr, proj.ip, proj.l, proj.m);
        }

    }                           /* end for */

    delete RT3;
    RmgTimer *RT4= new RmgTimer("Weight: bcast");

    int root;
    for(int iproj = 0; iproj < this->num_projectors; iproj++)
    {
        proj = proj_iter[iproj];
        root = iproj % pct.grid_npes;
        size = this->nldim * this->nldim * this->nldim;

        for(int kpt = 0; kpt <ct.num_kpts_pe; kpt++)
        {
            betaptr = &this->forward_beta[kpt *this->num_projectors *size + proj.proj_index * size];
            MPI_Bcast(betaptr, 2*size, MPI_DOUBLE, root, pct.grid_comm);
        }
    }

    delete [] this->phase;

    if(ct.stress || ct.LOPTICS || ct.forceflag == TDDFT)
    {
        this->InitLocalizedWeight_xyz();
    }

    delete RT4;

}
void SPECIES::InitLocalizedWeight_xyz (void)
{

    fftw_complex *betaptr;
    std::complex<double> *phaseptr;

    typedef struct {int species; int ip; int l; int m; int proj_index;} PROJ_INFO;
    PROJ_INFO proj;
    std::vector<PROJ_INFO> proj_iter;


    RmgTimer RT0("Weight_xyz");
    // get tot number of projectors and their information

    int nldim_max = ct.max_nldim;

    RmgTimer *RT1= new RmgTimer("Weight_xyz: phase and set");

    int size = this->nldim * this->nldim * this->nldim;
    this->phase = new fftw_complex[size * ct.num_kpts_pe];

    phaseptr = (std::complex<double> *)this->phase;
    // phase: multiply exp(-ik.r) in real space
    GetPhaseSpecies(this, phaseptr);
    /*Loop over all betas to calculate num of projectors for given species */
    int prjcount = 0;
    for (int ip = 0; ip < this->nbeta; ip++)
    {

        for(int m = 0; m < 2*this->llbeta[ip]+1; m++)
        {
            proj.ip = ip;
            proj.l = this->llbeta[ip];
            proj.m = m;
            proj.proj_index = prjcount;
            proj_iter.push_back(proj);
            prjcount++;
        }
    }

    size = this->nldim * this->nldim * this->nldim;
    this->num_projectors = prjcount;

    /*This array will store forward fourier transform on the coarse grid for all betas */
    if(this->forward_beta_r[0]) fftw_free(this->forward_beta_r[0]);
    if(this->forward_beta_r[1]) fftw_free(this->forward_beta_r[1]);
    if(this->forward_beta_r[2]) fftw_free(this->forward_beta_r[2]);
    this->forward_beta_r[0] = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * this->num_projectors * size * ct.num_kpts_pe);
    this->forward_beta_r[1] = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * this->num_projectors * size * ct.num_kpts_pe);
    this->forward_beta_r[2] = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * this->num_projectors * size * ct.num_kpts_pe);

    if (this->forward_beta_r[2] == NULL)
        throw RmgFatalException() << "cannot allocate mem "<< " at line " << __LINE__ << "\n";

    delete RT1;

    RmgTimer *RT3= new RmgTimer("Weight_xyz: proj cal");
    for(int iproj = pct.gridpe; iproj < this->num_projectors; iproj+=pct.grid_npes)
    {
        proj = proj_iter[iproj];

        // if this->nldim > get_NX_GRID, folding of neighbor cells are needed. 

        for(int kpt = 0; kpt <ct.num_kpts_pe; kpt++)
        {
            size_t index_ptr = kpt *this->num_projectors *size + proj.proj_index * size;
            std::complex<double> *betaptr_r[3];

            betaptr_r[0] = (std::complex<double> *)&this->forward_beta_r[0][index_ptr];
            betaptr_r[1] = (std::complex<double> *)&this->forward_beta_r[1][index_ptr];
            betaptr_r[2] = (std::complex<double> *)&this->forward_beta_r[2][index_ptr];
            for(int idx = 0;idx < size; idx++)
            {
                betaptr_r[0][idx] = 0.0;
                betaptr_r[1][idx] = 0.0;
                betaptr_r[2][idx] = 0.0;
            }
            phaseptr = (std::complex<double> *) &this->phase[kpt * this->nldim * this->nldim * this->nldim];
            int l2mi = proj.l * proj.l + proj.m;

            for(int l2mj = 1; l2mj < 4; l2mj++)     // index for cubic harmonics x, y, z
            {
                InitWeight_xyz_One(this, betaptr_r[l2mj-1], phaseptr, proj.ip, l2mi, l2mj);
            }
        }
    }


    delete RT3;
    RmgTimer *RT4= new RmgTimer("Weight: bcast");

    int root;
    for(int iproj = 0; iproj < this->num_projectors; iproj++)
    {
        proj = proj_iter[iproj];
        root = iproj % pct.grid_npes;
        size = this->nldim * this->nldim * this->nldim;

        for(int kpt = 0; kpt <ct.num_kpts_pe; kpt++)
        {
            betaptr = &this->forward_beta_r[0][kpt *this->num_projectors *size + proj.proj_index * size];
            MPI_Bcast(betaptr, 2*size, MPI_DOUBLE, root, pct.grid_comm);
            betaptr = &this->forward_beta_r[1][kpt *this->num_projectors *size + proj.proj_index * size];
            MPI_Bcast(betaptr, 2*size, MPI_DOUBLE, root, pct.grid_comm);
            betaptr = &this->forward_beta_r[2][kpt *this->num_projectors *size + proj.proj_index * size];
            MPI_Bcast(betaptr, 2*size, MPI_DOUBLE, root, pct.grid_comm);
        }
    }

    delete [] this->phase;


    delete RT4;

}


void InitWeight_xyz_One (SPECIES * sp, std::complex<double> *rtptr, std::complex<double> *phaseptr, int ip, int l2mi, int l2mj)
{

    double ax[3];
    double t1;
    double tpiba = 2.0 * PI / Rmg_L.celldm[0];
    double tpiba2 = tpiba * tpiba;
    double scale = 2.0 * PI / sp->prj_pwave->L->celldm[0];

    std::complex<double> I_t(0.0, 1.0);
    std::complex<double> *IL =  new std::complex<double>[ct.max_l+2];
    for(int L = 0; L <  ct.max_l+2; L++) IL[L] = std::pow(-I_t, L);
    int lmax = ct.max_l + 1;
    int num_lm = (lmax + 1) * (lmax + 1);
    int num_LM2 = (2*lmax + 1) * (2*lmax + 1);

    std::vector<int> lpx(num_lm * num_lm);
    std::vector<int> lpl(num_lm * num_lm  * num_LM2);
    std::vector<double> ap(num_LM2 * num_lm * num_lm);

    InitClebschGordan(lmax, ap.data(), lpx.data(), lpl.data());

    /* nl[xyz]fdim is the size of the non-local box in the high density grid */
    int size = sp->nldim * sp->nldim * sp->nldim;

    std::complex<double> *weptr = new std::complex<double>[size]();
    std::complex<double> *gwptr = new std::complex<double>[size];


    if(!weptr || !gwptr)
        throw RmgFatalException() << "cannot allocate mem "<< " at line " << __LINE__ << "\n";

    int ixx, iyy, izz;
    double gval, gcut;

    int xdim = sp->nldim;
    int ydim = sp->nldim;
    int zdim = sp->nldim;

    gcut = sqrt(coarse_pwaves->gcut*tpiba2); // pwave structures store the squared magnitude

    double vol = Rmg_L.get_omega() * (double)xdim * (double)ydim * (double)zdim /
        (double)Rmg_G->get_GLOBAL_BASIS(1);


    //  shift the atom to center with a phase 2pi *(xdim/2)/xdim

    std::complex<double> phase;

    phase = 2.0 * PI * ((xdim+1)/2)/xdim * I_t;
    phase = std::exp(phase);

    for (int ix = -xdim/2; ix < xdim/2+1; ix++)
    {
        ixx = (ix + xdim)%xdim;
        for (int iy = -ydim/2; iy < ydim/2+1; iy++)
        {
            iyy = (iy + ydim)%ydim;
            for (int iz = -zdim/2; iz < zdim/2+1; iz++)
            {
                izz = (iz + zdim)%zdim;
                int idx1 = ixx * ydim * zdim + iyy * zdim + izz;

                ax[0] = scale * sp->prj_pwave->g[idx1].a[0];
                ax[1] = scale * sp->prj_pwave->g[idx1].a[1];
                ax[2] = scale * sp->prj_pwave->g[idx1].a[2];
                gval = scale*sqrt(sp->prj_pwave->gmags[idx1]);

                if(gval >= gcut) continue;

                for (int LM = 0; LM < lpx[l2mi * num_lm + l2mj]; LM++)
                {
                    int L2M = lpl[(l2mi * num_lm + l2mj) * num_LM2 + LM];   // L*L + M for one LM harmonic function 

                    int L, M;
                    if(L2M == 0)
                        L = 0;
                    else if (L2M < 4)
                        L = 1;
                    else if (L2M < 9)
                        L = 2;
                    else
                        L = (int)sqrt(L2M + 0.1);

                    M = L2M - L * L;
                    double t2 = AtomicInterpolateInline_Ggrid(sp->rbeta_g[ip][L], gval);
                    weptr[idx1] += IL[L] * Ylm(L, M, ax) * t2 * ap[L2M * num_lm * num_lm + l2mi * num_lm + l2mj]/vol;
                }
            }
        }
    }

    // shift atom to the center instead of corner
    for(int ix = 0; ix < xdim; ix++)
        for(int iy = 0; iy < ydim; iy++)
            for(int iz = 0; iz < zdim; iz++)
            {
                int idx = ix * ydim * zdim + iy * zdim + iz;
                weptr[idx] *= std::pow(phase, ix+iy+iz);
                //        if( (ix + iy + iz) %2 ) weptr[idx] *=-1.0;
            }



    //    for(int ix = 0; ix < sp->nldim; ix++) printf("\n aaa %d %e %e", ix, weptr[ix]);

    RmgTimer *RT1 = new RmgTimer("weight fft_nldim");
    sp->prj_pwave->FftInverse(weptr, gwptr);
    delete RT1;

    RmgTimer *RT2 = new RmgTimer("weight fold");

    size = xdim * ydim * zdim;

    for(int idx = 0; idx < size; idx++) weptr[idx] = 0.0;
    for (int ix = 0; ix < sp->nldim; ix++)
    {
        int ixx = (ix-sp->nldim/2 + xdim/2 + 20 * xdim) % xdim;

        for (int iy = 0; iy < sp->nldim; iy++)
        {
            int iyy = (iy-sp->nldim/2 + ydim/2 + 20 * ydim) % ydim;
            for (int iz = 0; iz < sp->nldim; iz++)
            {
                int izz = (iz-sp->nldim/2 + zdim/2 + 20 * zdim) % zdim;
                int idx1 = ix * sp->nldim * sp->nldim + iy * sp->nldim + iz;
                int idx = ixx * ydim * zdim + iyy * zdim + izz;
                weptr[idx] += gwptr[idx1] * phaseptr[idx1]/double(size)/std::sqrt(3.0/fourPI);
            }
        }
    }


    delete RT2;
    RmgTimer *RT3 = new RmgTimer("weight fft_forward");
    sp->prj_pwave->FftForward(weptr, rtptr);
    delete RT3;

    // shift atom to the center instead of corner

    delete []gwptr;
    delete []weptr;
}

