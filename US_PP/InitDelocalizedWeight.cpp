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
    std::complex<double> phase = PI * I_t;
    phase = std::exp(phase);
    int max_pbasis = 0;


    RmgTimer *RT1= new RmgTimer("Weight: phase and set");
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

        if (sp->forward_beta == NULL)
            throw RmgFatalException() << "cannot allocate mem "<< " at line " << __LINE__ << "\n";

    }

    std::complex<double> *weptr = new std::complex<double>[max_pbasis];

    RmgTimer *RT3= new RmgTimer("Weight: proj cal");
    for(int iproj = 0; iproj < tot_proj; iproj++)
    {
        proj = proj_iter[iproj];
        std::complex<double> IL = std::pow(-I_t, proj.l);
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
            std::complex<double> *betaptr = (std::complex<double> *)&sp->forward_beta[kpt *sp->num_projectors * pbasis + proj.proj_index * pbasis];
            for(int idx = 0;idx < pbasis;idx++)
            {
                if(!pwave->gmask[idx]) continue;
                weptr[idx] = std::complex<double>(0.0,0.0);
                ax[0] = pwave->g[idx].a[0] * tpiba;
                ax[1] = pwave->g[idx].a[1] * tpiba;
                ax[2] = pwave->g[idx].a[2] * tpiba;

                ax[0] += ct.kp[kpt1].kvec[0];
                ax[1] += ct.kp[kpt1].kvec[1];
                ax[2] += ct.kp[kpt1].kvec[2];

                double gval = sqrt(ax[0]*ax[0] + ax[1]*ax[1] + ax[2]*ax[2]);
                if(gval >= gcut) continue;
                double t1 = AtomicInterpolateInline_Ggrid(sp->beta_g[proj.ip], gval);
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
    delete [] weptr;

} /* end InitDelocalizedWeight */


