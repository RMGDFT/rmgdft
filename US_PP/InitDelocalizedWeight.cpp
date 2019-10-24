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
void InitDelocalizedWeight (void)
{
    RmgTimer RT0("Weight");


    PROJ_INFO proj;
    std::vector<PROJ_INFO> proj_iter;

    // get tot number of projectors and their information
    int tot_proj = 0;
    

    std::complex<double> I_t(0.0, 1.0);
    std::complex<double> *IL =  new std::complex<double>[ct.max_l+2];
    for(int L = 0; L <  ct.max_l+2; L++) IL[L] = std::pow(-I_t, L);
    std::complex<double> phase = PI * I_t;
    phase = std::exp(phase);
    int max_pbasis = 0;

    int num_lm = (ct.max_l + 1) * (ct.max_l+1);
    int num_LM2 = (2*ct.max_l + 1) * (2*ct.max_l+1);

    int *lpx = new int[num_lm * num_lm];
    int *lpl = new int[num_lm * num_lm  * num_LM2];
    double *ap = new double[num_LM2 * num_lm * num_lm];
    //ylm = new double[25];

    InitClebschGordan(ct.max_l, ap, lpx, lpl);

    for (int isp = 0; isp < ct.num_species; isp++)
    {
        /* Get species type */
        SPECIES *sp = &Species[isp];
        Pw *pwave = sp->prj_pwave;
        int pbasis = pwave->Grid->get_P0_BASIS(1);
        max_pbasis = std::max(max_pbasis, pbasis);

        /*Loop over all betas to calculate num of projectors for given species */
        int prjcount = 0;
        for (int ip = 0; ip < sp->nbeta; ip++)
        {
            for(int m = 0; m < 2*sp->llbeta[ip]+1; m++)
            {
                proj.species = isp;
                proj.ip = ip;
                proj.l = sp->llbeta[ip];
                proj.m = m;
                proj.proj_index = prjcount;
                proj_iter.push_back(proj);
                prjcount++;
            }
        }

        sp->num_projectors = prjcount;
        tot_proj += prjcount;


        /*This array will store forward fourier transform on the coarse grid for all betas of this species */
        sp->forward_beta = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * sp->num_projectors * pbasis * ct.num_kpts_pe);
        sp->forward_beta_r[0] = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * sp->num_projectors * pbasis * ct.num_kpts_pe);
        sp->forward_beta_r[1] = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * sp->num_projectors * pbasis * ct.num_kpts_pe);
        sp->forward_beta_r[2] = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * sp->num_projectors * pbasis * ct.num_kpts_pe);

        if (sp->forward_beta == NULL)
            throw RmgFatalException() << "cannot allocate mem "<< " at line " << __LINE__ << "\n";

    }


    for(int iproj = 0; iproj < tot_proj; iproj++)
    {
        proj = proj_iter[iproj];
        SPECIES *sp = &Species[proj.species];
        Pw *pwave = sp->prj_pwave;
        int dimx = pwave->Grid->get_PX0_GRID(1);
        int dimy = pwave->Grid->get_PY0_GRID(1);
        int dimz = pwave->Grid->get_PZ0_GRID(1);
        int ixstart = pwave->Grid->get_PX_OFFSET(1);
        int iystart = pwave->Grid->get_PY_OFFSET(1);
        int izstart = pwave->Grid->get_PZ_OFFSET(1);
        double vol = pwave->L->get_omega();
        int pbasis = pwave->Grid->get_P0_BASIS(1);
        double tpiba = 2.0 * PI / Rmg_L.celldm[0];
        double tpiba2 = tpiba * tpiba;
        double gcut = sqrt(ct.filter_factor*pwave->gcut*tpiba2);
        double ax[3];

        for(int kpt = 0; kpt <ct.num_kpts_pe; kpt++)
        {
            int kpt1 = kpt + pct.kstart;
            size_t index_ptr = kpt *sp->num_projectors * pbasis + proj.proj_index * pbasis;
            std::complex<double> *betaptr = (std::complex<double> *)&sp->forward_beta[index_ptr];
            std::complex<double> *betaptr_r[3];
        
            betaptr_r[0] = (std::complex<double> *)&sp->forward_beta_r[0][index_ptr];
            betaptr_r[1] = (std::complex<double> *)&sp->forward_beta_r[1][index_ptr];
            betaptr_r[2] = (std::complex<double> *)&sp->forward_beta_r[2][index_ptr];

            for(int idx = 0;idx < pbasis;idx++)
            {
                betaptr[idx] = 0.0;
                betaptr_r[0][idx] = 0.0;
                betaptr_r[1][idx] = 0.0;
                betaptr_r[2][idx] = 0.0;
            }

            for(int idx = 0;idx < pbasis;idx++)
            {
                if(!pwave->gmask[idx]) continue;
                ax[0] = pwave->g[idx].a[0] * tpiba;
                ax[1] = pwave->g[idx].a[1] * tpiba;
                ax[2] = pwave->g[idx].a[2] * tpiba;

                ax[0] += ct.kp[kpt1].kvec[0];
                ax[1] += ct.kp[kpt1].kvec[1];
                ax[2] += ct.kp[kpt1].kvec[2];

                double gval = sqrt(ax[0]*ax[0] + ax[1]*ax[1] + ax[2]*ax[2]);
                if(gval >= gcut) continue;
                double t1 = AtomicInterpolateInline_Ggrid(sp->beta_g[proj.ip], gval);
                betaptr[idx] = IL[proj.l] * Ylm(proj.l, proj.m, ax) * t1;

                // l2m_i: l*l + m for the first angular momentum
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

                        double t2 = AtomicInterpolateInline_Ggrid(sp->rbeta_g[proj.ip][L], gval);
                        betaptr_r[l2mj-1][idx] += IL[L] * Ylm(L, M, ax) * t2 * ap[L2M * num_lm * num_lm + l2mi * num_lm + l2mj];
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
                        betaptr[idx] *= phaseshift/vol;
                        betaptr_r[0][idx] *= phaseshift/vol / std::sqrt(3.0/fourPI);
                        betaptr_r[1][idx] *= phaseshift/vol / std::sqrt(3.0/fourPI);
                        betaptr_r[2][idx] *= phaseshift/vol / std::sqrt(3.0/fourPI);

                        // sqrt(3/4pi) is the normlaized constant for hamonic x, y, z and we only want non-normalized x, y, z
                    }
                }
            }

        }

    }  // end for

    delete [] IL;

} /* end InitDelocalizedWeight */
