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
#include "packfuncs.h"


/*This sets loop over species does forward fourier transofrm, finds and stores whatever is needed so that
 * only backwards Fourier transform is needed in the calculation*/
void InitDelocalizedWeight (void)
{
    RmgTimer RT0("Weight");
    Mgrid MG(&Rmg_L, Rmg_T);

    int density = 2;
    std::complex<double> *phaseptr;

    typedef struct {int species; int ip; int l; int m; int proj_index;} PROJ_INFO;
    PROJ_INFO proj;
    std::vector<PROJ_INFO> proj_iter;
    double ax[3];


    // get tot number of projectors and their information
    int tot_proj = 0;
    
    int dimx = Rmg_G->get_PX0_GRID(1);
    int dimy = Rmg_G->get_PY0_GRID(1);
    int dimz = Rmg_G->get_PZ0_GRID(1);
    int fdimx = Rmg_G->get_PX0_GRID(density);
    int fdimy = Rmg_G->get_PY0_GRID(density);
    int fdimz = Rmg_G->get_PZ0_GRID(density);
    int fixstart = Rmg_G->get_PX_OFFSET(density);
    int fiystart = Rmg_G->get_PY_OFFSET(density);
    int fizstart = Rmg_G->get_PZ_OFFSET(density);
    int NX_GRID = Rmg_G->get_NX_GRID(1);
    int NY_GRID = Rmg_G->get_NY_GRID(1);
    int NZ_GRID = Rmg_G->get_NZ_GRID(1);
    int FNX_GRID = Rmg_G->get_NX_GRID(density);
    int FNY_GRID = Rmg_G->get_NY_GRID(density);
    int FNZ_GRID = Rmg_G->get_NZ_GRID(density);
    int global_basis = NX_GRID * NY_GRID * NZ_GRID;
    int fglobal_basis = FNX_GRID * FNY_GRID * FNZ_GRID;
    int pbasis = Rmg_G->get_P0_BASIS(1);
    int fpbasis = Rmg_G->get_P0_BASIS(density);
    double hxx = get_hxgrid() * get_xside();
    double hyy = get_hygrid() * get_yside();
    double hzz = get_hzgrid() * get_zside();
    double vol = hxx * hyy * hzz * global_basis;
    hxx /= (double)density;
    hyy /= (double)density;
    hzz /= (double)density;

    std::complex<double> I_t(0.0, 1.0);
    std::complex<double> phase = PI * I_t;
    phase = std::exp(phase);


    std::complex<double> *weptr = new std::complex<double>[fpbasis];
    std::complex<double> *work1 = new std::complex<double>[(fdimx+2)*(fdimy+2)*(fdimz+2)];
    std::complex<double> *work2 = new std::complex<double>[(fdimx+2)*(fdimy+2)*(fdimz+2)];

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

    }

    RmgTimer *RT3= new RmgTimer("Weight: proj cal");
    for(int iproj = 0; iproj < tot_proj; iproj++)
    {
        proj = proj_iter[iproj];
        std::complex<double> IL = std::pow(-I_t, proj.l);
        SPECIES *sp = &ct.sp[proj.species];
        double gcut = PI / hxx + 1.0e-6;

        for(int kpt = 0; kpt <ct.num_kpts_pe; kpt++)
        {
//              phaseptr = (std::complex<double> *) &sp->phase[kpt * pbasis];
            std::complex<double> *betaptr = (std::complex<double> *)&sp->forward_beta[kpt *sp->num_projectors * pbasis + proj.proj_index * pbasis];
            for(int idx = 0;idx < fpbasis;idx++) weptr[idx] = std::complex<double>(0.0,0.0);
            for(int idx = 0;idx < fpbasis;idx++)
            {
                ax[0] = 2.0*PI*beta_pwaves->g[idx].a[0] / (hxx * FNX_GRID);
                ax[1] = 2.0*PI*beta_pwaves->g[idx].a[1] / (hyy * FNY_GRID);
                ax[2] = 2.0*PI*beta_pwaves->g[idx].a[2] / (hzz * FNZ_GRID);
                double gval = sqrt(ax[0]*ax[0] + ax[1]*ax[1] + ax[2]*ax[2]);
                if(gval >= 0.2*gcut) continue;
                double t1 = AtomicInterpolateInline_Ggrid(sp->beta_g[proj.ip], gval);
                weptr[idx] = IL * Ylm(proj.l, proj.m, beta_pwaves->g[idx].a) * t1;
            }

            // Shift atom to the center instead of corner.
            for(int ix = 0; ix < fdimx; ix++)
            {
                for(int iy = 0; iy < fdimy; iy++)
                {
                    for(int iz = 0; iz < fdimz; iz++)
                    {
                        int idx = ix * fdimy * fdimz + iy * fdimz + iz;
                        weptr[idx] *= std::pow(phase, ix + fixstart + iy + fiystart + iz + fizstart);
                    }
                }
            }

            // Transform to real space
            PfftInverse(weptr, weptr, *beta_pwaves);
            int ixoff, iyoff, izoff;
            int dx2 = MG.MG_SIZE (fdimx, 0, Rmg_G->get_NX_GRID(density), Rmg_G->get_PX_OFFSET(density), Rmg_G->get_PX0_GRID(density), &ixoff, ct.boundaryflag);
            int dy2 = MG.MG_SIZE (fdimy, 0, Rmg_G->get_NY_GRID(density), Rmg_G->get_PY_OFFSET(density), Rmg_G->get_PY0_GRID(density), &iyoff, ct.boundaryflag);
            int dz2 = MG.MG_SIZE (fdimz, 0, Rmg_G->get_NZ_GRID(density), Rmg_G->get_PZ_OFFSET(density), Rmg_G->get_PZ0_GRID(density), &izoff, ct.boundaryflag);

            CPP_pack_ptos (work1, weptr, fdimx, fdimy, fdimz);
            Rmg_T->trade_images (work1, fdimx, fdimy, fdimz, FULL_TRADE);
            MG.mg_restrict (work1, work2, fdimx, fdimy, fdimz, dx2, dy2, dz2, ixoff, iyoff, izoff);
            CPP_pack_stop (work2, betaptr, dimx, dimy, dimz);

            for(int idx=0;idx < pbasis;idx++) betaptr[idx] =  betaptr[idx]/vol/(double)global_basis;
            PfftForward(betaptr, betaptr, *coarse_pwaves);

        }

    }  // end for
    delete RT3;

    delete RT1;
    delete [] work2;
    delete [] work1;
    delete [] weptr;

} /* end InitDelocalizedWeight */
