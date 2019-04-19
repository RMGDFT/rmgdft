/************************** SVN Revision Information **************************
 **    $Id: UpdatePot.cpp 3198 2015-09-02 21:27:36Z luw $    **
 ******************************************************************************/


#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "rmgtypedefs.h"
#include "params.h"
#include "typedefs.h"
#include "RmgTimer.h"
#include "transition.h"
#include "vhartree.h"
#include "packfuncs.h"
#include "blas.h"
#include "RmgParallelFft.h"

float *vh_init;

double VhDriver(double *rho, double *rhoc, double *vh, double *vh_ext, double rms_target)
{
    int FP0_BASIS = Rmg_G->get_P0_BASIS(Rmg_G->default_FG_RATIO);
    int dimx = Rmg_G->get_PX0_GRID(Rmg_G->default_FG_RATIO);
    int dimy = Rmg_G->get_PY0_GRID(Rmg_G->default_FG_RATIO);
    int dimz = Rmg_G->get_PZ0_GRID(Rmg_G->default_FG_RATIO);


    double *rho_tot = new double[FP0_BASIS];
    double residual;

    for(int i = 0;i < FP0_BASIS;i++) rho_tot[i] = rho[i];
    if(ct.spin_flag) {
        // opposite spin array is stored sequentially
        for(int i = 0;i < FP0_BASIS;i++) rho_tot[i] += rho[i + FP0_BASIS];
    }

    /* Subtract off compensating charges from rho */
    for (int i = 0; i < FP0_BASIS; i++)
        rho_tot[i] = rho_tot[i] - rhoc[i];

    // Make sure it is completely neutral
    int global_basis = Rmg_G->get_GLOBAL_BASIS(Rmg_G->default_FG_RATIO);
    double sum = 0.0;
    for(int i = 0;i < FP0_BASIS;i++) sum += rho_tot[i];
    MPI_Allreduce(MPI_IN_PLACE, &sum, 1, MPI_DOUBLE, MPI_SUM, pct.grid_comm);
    sum = sum / (double)global_basis;
    for(int i = 0;i < FP0_BASIS;i++) rho_tot[i] -= sum;


    if(ct.poisson_solver == POISSON_PFFT_SOLVER)
    {
        RmgTimer *RT2 = new RmgTimer("VhPfft");
        VhPfft(rho_tot, NULL, vh);
        delete(RT2);
        residual = 0.0;  // FFT solver is exact
    }
    else
    {

        size_t coarse_size = FP0_BASIS;
//        if(ct.poi_parm.levels > 0) coarse_size /= 8;
        if(ct.poi_parm.levels > 1) coarse_size /= 8;
        if(ct.poi_parm.levels > 2) coarse_size /= 8;
        if(ct.poi_parm.levels > 3) coarse_size /= 8;
        if(!vh_init) 
        {
            vh_init = new float[coarse_size]();
        }

        RmgTimer *RT1 = new RmgTimer("VhMg");
        residual = vh_fmg (Rmg_G, &Rmg_L, Rmg_T, rho_tot, vh_ext,
                 ct.hartree_min_sweeps, ct.hartree_max_sweeps, ct.poi_parm.levels, ct.poi_parm.gl_pre,
                 ct.poi_parm.gl_pst, ct.poi_parm.mucycles, rms_target,
                 ct.poi_parm.gl_step, ct.poi_parm.sb_step, ct.boundaryflag, Rmg_G->get_default_FG_RATIO(), vh_init, ct.verbose);
        /* Pack the portion of the hartree potential used by the wavefunctions
         * back into the wavefunction hartree array. */
        CPP_pack_dtos (Rmg_G, vh, vh_ext, dimx, dimy, dimz, ct.boundaryflag);
//        if(ct.filter_dpot && (Rmg_G->default_FG_RATIO > 1)) 
//            FftFilter(vh, *fine_pwaves, sqrt(ct.filter_factor) / (double)ct.FG_RATIO, LOW_PASS);

        delete(RT1);
    }


    delete [] rho_tot;
    return residual;
}
