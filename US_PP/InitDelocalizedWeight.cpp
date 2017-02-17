#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>

#include <complex>

#include "make_conf.h"
#include "const.h"
#include "grid.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "Kpoint.h"
#include "common_prototypes.h"
#include "common_prototypes1.h"
#include "transition.h"
#include "RmgException.h"
#include "RmgShm.h"
#include "RmgParallelFft.h"
#include "AtomicInterpolate.h"


/*This sets loop over species does forward fourier transofrm, finds and stores whatever is needed so that
 * only backwards Fourier transform is needed in the calculation*/
void InitWeightDelocalized (void)
{
    RmgTimer RT0("Weight");

    std::complex<double> *phaseptr;

    typedef struct {int species; int ip; int l; int m; int proj_index;} PROJ_INFO;
    PROJ_INFO proj;
    std::vector<PROJ_INFO> proj_iter;


    // get tot number of projectors and their information
    int tot_proj = 0;
    
    int dimx = get_PX0_GRID();
    int dimy = get_PY0_GRID();
    int dimz = get_PZ0_GRID();
    int ixstart = get_PX_OFFSET();
    int iystart = get_PY_OFFSET();
    int izstart = get_PZ_OFFSET();
    int NX_GRID = get_NX_GRID();
    int NY_GRID = get_NY_GRID();
    int NZ_GRID = get_NZ_GRID();
    int global_basis = NX_GRID * NY_GRID * NZ_GRID;
    int pbasis = get_P0_BASIS();
    double hxx = get_hxgrid() * get_xside();
    double hyy = get_hygrid() * get_yside();
    double hzz = get_hzgrid() * get_zside();
    double vol = hxx * hyy * hzz * global_basis;
    std::complex<double> I_t(0.0, 1.0);
    std::complex<double> phase = -PI * I_t;
    phase = std::exp(phase);


    std::complex<double> *weptr = new std::complex<double>[pbasis];

    RmgTimer *RT1= new RmgTimer("Weight: phase and set");
    for (int isp = 0; isp < ct.num_species; isp++)
    {
        /* Get species type */
        SPECIES *sp = &ct.sp[isp];

//        sp->phase = new fftw_complex[pbasis * ct.num_kpts_pe];
//        phaseptr = (std::complex<double> *)sp->phase;
//        GetPhaseSpecies(sp, phaseptr);
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

        RmgTimer *RT3= new RmgTimer("Weight: proj cal");
        for(int iproj = 0; iproj < sp->num_projectors; iproj++)
        {
            proj = proj_iter[iproj];
            std::complex<double> IL = std::pow(-I_t, proj.l);
            sp = &ct.sp[proj.species];

            for(int kpt = 0; kpt <ct.num_kpts_pe; kpt++)
            {
//              phaseptr = (std::complex<double> *) &sp->phase[kpt * pbasis];
                std::complex<double> *betaptr = (std::complex<double> *)&sp->forward_beta[kpt *sp->num_projectors * pbasis + proj.proj_index * pbasis];

                for(int idx = 0;idx < pbasis;idx++)
                {
                    double gval = coarse_pwaves->gmags[idx];
                    if(gval > 0.8*coarse_pwaves->gcut) continue;
                    double t1 = AtomicInterpolateInline_Ggrid(&sp->beta_g[proj.ip][0], sqrt(gval)/PI);
                    weptr[idx] = IL * Ylm(proj.l, proj.m, coarse_pwaves->g->a) * t1;
                }

                // shift atom to the center instead of corner
                for(int ix = 0; ix < dimx; ix++)
                {
                    for(int iy = 0; iy < dimy; iy++)
                    {
                        for(int iz = 0; iz < dimz; iz++)
                        {
                            int idx = ix * dimy * dimz + iy * dimz + iz;
                            //weptr[idx] *= std::pow(phase, ix + ixstart + iy + iystart + iz + izstart);
                            if( !(ix + iy + iz) %2 ) weptr[idx] *=-1.0;
                        }
                    }
                }

                for(int idx=0;idx < pbasis;idx++) betaptr[idx] =  weptr[idx]/vol;

            }

        }  // end for
        delete RT3;
    }

    delete RT1;
    delete [] weptr;

} /* end InitWeightDelocalized */
