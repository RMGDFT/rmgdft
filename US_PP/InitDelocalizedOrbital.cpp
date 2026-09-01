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
#include "RmgParallelFft.h"
#include "AtomicInterpolate.h"
#include "packfuncs.h"


/*This sets loop over species does forward fourier transofrm, finds and stores whatever is needed so that
 * only backwards Fourier transform is needed in the calculation*/
void SPECIES::InitDelocalizedOrbital (void)
{
    RmgTimer RT0("Orbital");
    double tpiba = 2.0 * PI / Rmg_L.celldm[0];
    double tpiba2 = tpiba * tpiba;


    PROJ_INFO proj;
    std::vector<PROJ_INFO> proj_iter;
    double ax[3];


    // get tot number of projectors and their information
    int tot_orbitals = 0;
    
    int dimx = Rmg_G->get_PX0_GRID(1);
    int dimy = Rmg_G->get_PY0_GRID(1);
    int dimz = Rmg_G->get_PZ0_GRID(1);
    int ixstart = Rmg_G->get_PX_OFFSET(1);
    int iystart = Rmg_G->get_PY_OFFSET(1);
    int izstart = Rmg_G->get_PZ_OFFSET(1);
    int pbasis = Rmg_G->get_P0_BASIS(1);
    double vol = Rmg_L.get_omega();

    std::complex<double> I_t(0.0, 1.0);
    std::complex<double> phase = PI * I_t;
    phase = std::exp(phase);


    std::complex<double> *weptr = new std::complex<double>[pbasis];

    RmgTimer *RT1= new RmgTimer("Orbital: phase and set");
    /*Loop over all radial atomic wavefunctions to calculate num of 3D orbitals for given species */
    int orbital_count = 0;
    for (int ip = 0; ip < this->num_atomic_waves; ip++)
    {

        // If the occupation condition changes here then it must change in
        // GetFdFactors as well.
        if(this->atomic_wave_oc[ip] > 0.0) {

            for(int m = 0; m < 2*this->atomic_wave_l[ip]+1; m++)
            {
                proj.ip = ip;
                proj.l = this->atomic_wave_l[ip];
                proj.m = m;
                proj.proj_index = orbital_count;
                proj_iter.push_back(proj);
                orbital_count++;
            }
        }
    }

    this->num_orbitals = orbital_count;
    tot_orbitals += orbital_count;

    /*This array will store forward fourier transform on the coarse grid for all atomic orbitals of this species */
    if(this->forward_orbital[0]) fftw_free(this->forward_orbital[0]);
    this->forward_orbital[0] = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * this->num_orbitals * pbasis * ct.num_kpts_pe);
    if(this->forward_orbital_gamma) fftw_free(this->forward_orbital_gamma);
    this->forward_orbital_gamma = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * this->num_orbitals * pbasis);

    if (this->forward_orbital[0] == NULL || this->forward_orbital_gamma == NULL)
        throw RmgFatalException() << "cannot allocate mem "<< " at line " << __LINE__ << "\n";



    double gcut = sqrt(coarse_pwaves->gcut*tpiba2);
    RmgTimer *RT3= new RmgTimer("Orbital: proj cal");
    for(int iproj = 0; iproj < tot_orbitals; iproj++)
    {
        proj = proj_iter[iproj];
        std::complex<double> IL = std::pow(-I_t, proj.l);

        for(int kpt = 0; kpt <ct.num_kpts_pe; kpt++)
        {
            int kpt1 = kpt + pct.kstart;
            std::complex<double> *betaptr = (std::complex<double> *)&this->forward_orbital[0][kpt *this->num_orbitals * pbasis + proj.proj_index * pbasis];
            for(int idx = 0;idx < pbasis;idx++)
            {
                if(!coarse_pwaves->gmask[idx]) continue;
                weptr[idx] = std::complex<double>(0.0,0.0);
                ax[0] = coarse_pwaves->g[idx].a[0] * tpiba;
                ax[1] = coarse_pwaves->g[idx].a[1] * tpiba;
                ax[2] = coarse_pwaves->g[idx].a[2] * tpiba;

                ax[0] += ct.kp[kpt1].kvec[0];
                ax[1] += ct.kp[kpt1].kvec[1];
                ax[2] += ct.kp[kpt1].kvec[2];

                double gval = sqrt(ax[0]*ax[0] + ax[1]*ax[1] + ax[2]*ax[2]);
                if(gval >= gcut) continue;
                double t1 = AtomicInterpolateInline_Ggrid(this->atomic_wave_g[proj.ip], gval);
                weptr[idx] = IL * Ylm(proj.l, proj.m, ax) * t1;

            }

            // Shift atom to the center instead of corner.
            for(int ix = 0; ix < dimx; ix++)
            {
                for(int iy = 0; iy < dimy; iy++)
                {
                    for(int iz = 0; iz < dimz; iz++)
                    {
                        int idx = ix * dimy * dimz + iy * dimz + iz;
                        weptr[idx] *= std::pow(phase, ix + ixstart + iy + iystart + iz + izstart);
                    }
                }
            }
            for(int idx=0;idx < pbasis;idx++) betaptr[idx] =  weptr[idx]/vol;

        }







        // initilize orbital for gamma point used in GetFdFactor to optimize the FD operator
        {
            std::complex<double> *betaptr = (std::complex<double> *)&this->forward_orbital_gamma[proj.proj_index * pbasis];
            for(int idx = 0;idx < pbasis;idx++)
            {
                if(!coarse_pwaves->gmask[idx]) continue;
                weptr[idx] = std::complex<double>(0.0,0.0);
                ax[0] = coarse_pwaves->g[idx].a[0] * tpiba;
                ax[1] = coarse_pwaves->g[idx].a[1] * tpiba;
                ax[2] = coarse_pwaves->g[idx].a[2] * tpiba;

                double gval = sqrt(ax[0]*ax[0] + ax[1]*ax[1] + ax[2]*ax[2]);
                if(gval >= gcut) continue;
                double t1 = AtomicInterpolateInline_Ggrid(this->atomic_wave_g[proj.ip], gval);
                weptr[idx] = IL * Ylm(proj.l, proj.m, ax) * t1;

            }

            // Shift atom to the center instead of corner.
            for(int ix = 0; ix < dimx; ix++)
            {
                for(int iy = 0; iy < dimy; iy++)
                {
                    for(int iz = 0; iz < dimz; iz++)
                    {
                        int idx = ix * dimy * dimz + iy * dimz + iz;
                        weptr[idx] *= std::pow(phase, ix + ixstart + iy + iystart + iz + izstart);
                    }
                }
            }
            for(int idx=0;idx < pbasis;idx++) betaptr[idx] =  weptr[idx]/vol;
        }

    }  // end for
    delete RT3;
    delete RT1;


    if(ct.stress || ct.LOPTICS || ct.forceflag == TDDFT)
    {

        RmgTimer *RT3= new RmgTimer("Orbital: r_orbital cal");
        /*This array will store forward fourier transform on the coarse grid for all atomic orbitals of this species */
        if(this->forward_orbital[1]) fftw_free(this->forward_orbital[1]);
        if(this->forward_orbital[2]) fftw_free(this->forward_orbital[2]);
        if(this->forward_orbital[3]) fftw_free(this->forward_orbital[3]);
        this->forward_orbital[1] = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * this->num_orbitals * pbasis * ct.num_kpts_pe);
        this->forward_orbital[2] = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * this->num_orbitals * pbasis * ct.num_kpts_pe);
        this->forward_orbital[3] = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * this->num_orbitals * pbasis * ct.num_kpts_pe);

        std::complex<double> *IL =  new std::complex<double>[ct.max_l+2];
        for(int L = 0; L <  ct.max_l+2; L++) IL[L] = std::pow(-I_t, L);
        std::complex<double> phase = PI * I_t;
        phase = std::exp(phase);
        int max_pbasis = 0;

        int lmax = ct.max_l + 1;
        int num_lm = (lmax + 1) * (lmax + 1);
        int num_LM2 = (2*lmax + 1) * (2*lmax + 1);

        std::vector<int> lpx(num_lm * num_lm);
        std::vector<int> lpl(num_lm * num_lm  * num_LM2);
        std::vector<double> ap(num_LM2 * num_lm * num_lm);

        InitClebschGordan(lmax, ap.data(), lpx.data(), lpl.data());

        double gcut = sqrt(coarse_pwaves->gcut*tpiba2);
        for(int iproj = 0; iproj < tot_orbitals; iproj++)
        {
            proj = proj_iter[iproj];

            for(int kpt = 0; kpt <ct.num_kpts_pe; kpt++)
            {
                int kpt1 = kpt + pct.kstart;
                std::complex<double> *betaptr[3];
                betaptr[0] = (std::complex<double> *)&this->forward_orbital[1][kpt *this->num_orbitals * pbasis + proj.proj_index * pbasis];
                betaptr[1] = (std::complex<double> *)&this->forward_orbital[2][kpt *this->num_orbitals * pbasis + proj.proj_index * pbasis];
                betaptr[2] = (std::complex<double> *)&this->forward_orbital[3][kpt *this->num_orbitals * pbasis + proj.proj_index * pbasis];
                std::fill(betaptr[0], betaptr[0] + pbasis, 0.0);
                std::fill(betaptr[1], betaptr[1] + pbasis, 0.0);
                std::fill(betaptr[2], betaptr[2] + pbasis, 0.0);

                for(int idx = 0;idx < pbasis;idx++)
                {
                    if(!coarse_pwaves->gmask[idx]) continue;
                    weptr[idx] = std::complex<double>(0.0,0.0);
                    ax[0] = coarse_pwaves->g[idx].a[0] * tpiba;
                    ax[1] = coarse_pwaves->g[idx].a[1] * tpiba;
                    ax[2] = coarse_pwaves->g[idx].a[2] * tpiba;

                    ax[0] += ct.kp[kpt1].kvec[0];
                    ax[1] += ct.kp[kpt1].kvec[1];
                    ax[2] += ct.kp[kpt1].kvec[2];

                    double gval = sqrt(ax[0]*ax[0] + ax[1]*ax[1] + ax[2]*ax[2]);
                    if(gval >= gcut) continue;

                    int l2mi = proj.l * proj.l + proj.m;
                    for(int l2mj = 1; l2mj < 4; l2mj++)     // index for cubic harmonics x, y, z
                    {
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
                            double t2 = AtomicInterpolateInline_Ggrid(this->r_atomic_wave_g[proj.ip][L].origin(), gval);
                            betaptr[l2mj-1][idx] += IL[L] * Ylm(L, M, ax) * t2 * ap[L2M * num_lm * num_lm + l2mi * num_lm + l2mj];
                        }

                    }
                }

                // Shift atom to the center instead of corner.
                for(int ix = 0; ix < dimx; ix++)
                {
                    for(int iy = 0; iy < dimy; iy++)
                    {
                        for(int iz = 0; iz < dimz; iz++)
                        {
                            int idx = ix * dimy * dimz + iy * dimz + iz;
                            std::complex<double> phaseshift =std::pow(phase, ix + ixstart + iy + iystart + iz + izstart);
                            betaptr[0][idx] *= phaseshift/vol / std::sqrt(3.0/fourPI);
                            betaptr[1][idx] *= phaseshift/vol / std::sqrt(3.0/fourPI);
                            betaptr[2][idx] *= phaseshift/vol / std::sqrt(3.0/fourPI);
                        }
                    }
                }

            }  // end for

        }

        delete RT3;
    }

    delete [] weptr;

} /* end InitDelocalizedOrbital */
