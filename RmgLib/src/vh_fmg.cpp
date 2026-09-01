/*
 *
 * Copyright (c) 1995, Emil Briggs
 * Copyright (C) 1998  Emil Briggs, Charles Brabec, Mark Wensell, 
 *                     Dan Sullivan, Chris Rapcewicz, Jerzy Bernholc
 * Copyright (C) 2001  Emil Briggs, Wenchang Lu,
 *                     Marco Buongiorno Nardelli,Charles Brabec, 
 *                     Mark Wensell,Dan Sullivan, Chris Rapcewicz,
 *                     Jerzy Bernholc
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the <organization> nor the
 *       names of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.

 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * 
*/


#include <iostream>
#include <omp.h>

#include "TradeImages.h"
#include "FiniteDiff.h"
#include "Mgrid.h"
#include "RmgSumAll.h"
#include "vhartree.h"
#include "rmg_error.h"
#include "packfuncs.h"
#include "boundary_conditions.h"
#include "RmgTimer.h"

#define         PI          3.14159265358979323


#define MAX_MG_LEVELS 8

/// Poisson solver that uses compact implicit (Mehrstellen) and multigrid techniques.
/// @param G Grid object that defines the layout of the 3-D grid and the MPI domains
/// @param rho Charge density. When using periodic boundary conditions the cell must be charge neutral.
/// @param vhartree Hartree potential corresponding to rho.
/// @param min_sweeps Minimum number of top level sweeps to perform.
/// @param max_sweeps Maximum number of top level sweeps to perform.
/// @param max_level Maximum number of multigrid levels to use.
/// @param global_presweeps Number of presweeps to use on the finest (0th) grid level.
/// @param global_postsweeps Number of postsweeps to use on the finest (0th) grid level.
/// @param mucycles Number of mu cycles (also known as W-cycles) to use in the multigrid solver.
/// @param rms_target Value for the root mean square residual at which to stop.
/// @param global_step Time step for the jacobi iteration on the finest (0th) grid level.
/// @param coarse_step Time step for the jacobi iteration on the coarse grid levels.
/// @param boundaryflag Type of boundary condition. Periodic is implemented internally.
/// @param density Density of the grid relative to the default grid
double vh_fmg (BaseGrid *G, Lattice *L, TradeImages *T, double * rho, double *vhartree,
                 int min_sweeps, int max_sweeps, int maxlevel, 
                 int global_presweeps, int global_postsweeps, int mucycles, 
                 double rms_target_in, double global_step, double coarse_step, int boundaryflag, int density,
                 float *vh_init, bool print_status)
{

    RmgTimer *RT0 = new RmgTimer("Hartree: init");
    double t1;
    double residual = 100.0;
    double rms_target = std::min(rms_target, 1.0e-6);
    rms_target = std::max(rms_target_in, 1.0e-10);
    Mgrid MG(L, T);

    if(maxlevel >= MAX_MG_LEVELS)
       rmg_error_handler(__FILE__, __LINE__, "Too many multigrid levels requested.");

    int dimx = G->get_PX0_GRID(density), dimy = G->get_PY0_GRID(density), dimz = G->get_PZ0_GRID(density);

    // Solve to a high degree of precision on the coarsest level
    int pbasis = dimx * dimy * dimz;
    int sbasis = (dimx + 2) * (dimy + 2) * (dimz + 2);

    /* Grab some memory for our multigrid structures */
    double *mgrhsarr = new double[2*sbasis];
    double *mglhsarr = new  double[sbasis];
    double *work = new  double[sbasis];
    double *sg_res = new  double[sbasis];

    float *mgrhsarr_f = new float[2*sbasis];
    float *mglhsarr_f = (float *)mglhsarr;
    float *sg_res_f = (float *)sg_res;

    int dx[MAX_MG_LEVELS];dx[0] = dimx;
    int dy[MAX_MG_LEVELS];dy[0] = dimy;
    int dz[MAX_MG_LEVELS];dz[0] = dimz;
    int ixoff, iyoff, izoff;
    double *mgrhsptr[MAX_MG_LEVELS];
    float *mgrhsptr_f[MAX_MG_LEVELS];


    // Set up RHS on finest level
    /* Multiply through by 4PI */
    t1 = -4.0 * PI;
    double rhs_avgcor = 0.0;
    for(int idx = 0;idx < pbasis;idx++) 
    {
        mgrhsarr[idx] = t1 * rho[idx];
        rhs_avgcor += mgrhsarr[idx];
    }
    rhs_avgcor = RmgSumAll(rhs_avgcor, T->get_MPI_comm());
    rhs_avgcor /= (double)G->get_GLOBAL_BASIS(density);
    for(int idx = 0;idx < pbasis;idx++) mgrhsarr[idx] -= rhs_avgcor;
    for(int idx = 0;idx < pbasis;idx++) mgrhsarr_f[idx] = (float)mgrhsarr[idx];


    // Restrict right hand side down to coarsest level
    size_t offset=pbasis;
    int dx2(0), dy2(0), dz2(0);
    mgrhsptr[0] = mgrhsarr;
    mgrhsptr_f[0] = mgrhsarr_f;

    for(int level=1;level <= maxlevel;level++)
    {
        dx2 = MG.MG_SIZE (dx[level-1], level-1, G->get_NX_GRID(density), G->get_PX_OFFSET(density), dimx, &ixoff, boundaryflag);
        dy2 = MG.MG_SIZE (dy[level-1], level-1, G->get_NY_GRID(density), G->get_PY_OFFSET(density), dimy, &iyoff, boundaryflag);
        dz2 = MG.MG_SIZE (dz[level-1], level-1, G->get_NZ_GRID(density), G->get_PZ_OFFSET(density), dimz, &izoff, boundaryflag);

        CPP_pack_ptos (work, mgrhsptr[level-1], dx[level-1], dy[level-1], dz[level-1]);
        T->trade_images (work, dx[level-1], dy[level-1], dz[level-1], FULL_TRADE);
        MG.mg_restrict (work, sg_res, dx[level-1], dy[level-1], dz[level-1], dx2, dy2, dz2, ixoff, iyoff, izoff);
      
        mgrhsptr[level] = &mgrhsarr[offset];
        mgrhsptr_f[level] = &mgrhsarr_f[offset];
        CPP_pack_stop (sg_res, mgrhsptr[level], dx2, dy2, dz2);

        // Make sure the restriction process didn't introduce a spurious source term in the RHS
        int global_basis = G->get_GLOBAL_BASIS(density) / std::round(pow(8.0, level));
        double *tptr=mgrhsptr[level];
        rhs_avgcor = 0.0;
        for(int idx=0;idx<dx2*dy2*dz2;idx++)rhs_avgcor += tptr[idx];
        rhs_avgcor = RmgSumAll(rhs_avgcor, T->get_MPI_comm());
        rhs_avgcor /= (double)global_basis;
        for(int idx=0;idx<dx2*dy2*dz2;idx++)tptr[idx] -= rhs_avgcor;
        //if((G->get_rank() == 0)) printf("FOR LEVEL=%d  SUM=%18.12e\n",level,rhs_avgcor);

        offset += dx2*dy2*dz2;

        dx[level] = dx2;
        dy[level] = dy2;
        dz[level] = dz2;

    } // end loop over levels

    for(int idx=0;idx < 2*sbasis;idx++)mgrhsarr_f[idx] = (float)mgrhsarr[idx];
    delete RT0;

    RmgTimer *RT1 = new RmgTimer("Hartree: SP solve");
    // Now solve from coarse grid to fine grid
    // If the calling routine has passed in an array to use for the coarse grid init
    // we use that and then save it for the next iteration. Otherwise just use 0.0 for
    // the coarse grid init.
    int pmaxlevel = dx[maxlevel]*dy[maxlevel]*dz[maxlevel];
    if(vh_init)
    {
        for(int ix=0;ix < pmaxlevel;ix++) mglhsarr_f[ix] = vh_init[ix];
    }
    else
    {
        for(int ix=0;ix < pmaxlevel;ix++) mglhsarr_f[ix] = 0.0;
    }

    for(int level=maxlevel;level >= 0;level--)
    {
        double lfactor = pow(2.0, (double)(level));
        coarse_vh (G, L, T, mgrhsptr_f[level], mglhsarr_f,
                 1, 2, maxlevel,
                 global_presweeps, global_postsweeps,
                 dx[level], dy[level], dz[level], level,
                 G->get_hxgrid(density)*lfactor, G->get_hygrid(density)*lfactor, G->get_hzgrid(density)*lfactor,
                 1.0e-8, global_step, coarse_step, boundaryflag, density, false, false);

        // Save coarse grid starting solution to use next time if vh_init is not null
        if((level == maxlevel) && vh_init) for(int ix=0;ix < dx2*dy2*dz2;ix++) vh_init[ix] = (float)mglhsarr_f[ix];
        if(level == 0) break;
        MG.mg_prolong_cubic (sg_res_f, mglhsarr_f, dx[level-1], dy[level-1], dz[level-1], dx[level], dy[level], dz[level], ixoff, iyoff, izoff);
        CPP_pack_stop (sg_res_f, mglhsarr_f, dx[level-1], dy[level-1], dz[level-1]);
    }
    delete RT1;

    RmgTimer *RT2 = new RmgTimer("Hartree: DP solve");
    for(int idx=0;idx < pbasis;idx++) work[idx] = (double)mglhsarr_f[idx];
    t1 = -4.0*PI;
    for(int idx=0;idx < pbasis;idx++) mgrhsarr[idx] = t1 * rho[idx];
    residual = coarse_vh (G, L, T, mgrhsarr, work,
             1, mucycles, maxlevel,
             global_presweeps, global_postsweeps,
             dimx, dimy, dimz, 0,
             G->get_hxgrid(density), G->get_hygrid(density), G->get_hzgrid(density),
             rms_target, global_step, coarse_step, boundaryflag, density, false, true);

    for(int idx=0;idx < pbasis;idx++) vhartree[idx] = work[idx];
    delete RT2;

    /* Release our memory */
    delete [] mgrhsarr_f;
    delete [] sg_res;
    delete [] work;
    delete [] mglhsarr;
    delete [] mgrhsarr;

    return residual;

} // vh_fmg



template <typename CalcType>
double coarse_vh (BaseGrid *G, Lattice *L, TradeImages *T, CalcType * rho, CalcType *vhartree,
                 int min_sweeps, int max_sweeps, int maxlevel, 
                 int global_presweeps, int global_postsweeps,
                 int dimx, int dimy, int dimz, int level,
                 double gridhx, double gridhy, double gridhz,
                 double rms_target, double global_step, double coarse_step, int boundaryflag, int density, bool print_status, bool setzero)
{

    int idx, its, cycles;
    double t1, vavgcor, diag=0.0, residual = 100.0, last_residual = 200.0;
    Mgrid MG(L, T);
    int global_basis = G->get_GLOBAL_BASIS(density) / pow(8.0, (double)level);

    /* Pre and post smoothings on each level */
    int poi_pre[MAX_MG_LEVELS] = { 4, 4, 4, 4, 4, 4, 4, 4};
    int poi_post[MAX_MG_LEVELS] = { 2, 2, 2, 2, 2, 2, 2, 2};

    int mu_cycles[MAX_MG_LEVELS] = {2, 2, 2, 2, 2, 2, 2, 2};
    if(maxlevel >= MAX_MG_LEVELS)
       rmg_error_handler(__FILE__, __LINE__, "Too many multigrid levels requested.");

    // Solve to a high degree of precision on the coarsest level
    int nits = global_presweeps + global_postsweeps;
    if(level == maxlevel) 
    {
        global_presweeps += global_postsweeps;
        nits = global_presweeps;
    }

    int factor =  (1 << maxlevel);
    int ixp = G->get_NX_GRID(density) / factor;
    int iyp = G->get_NY_GRID(density) / factor;
    int izp = G->get_NZ_GRID(density) / factor;
    factor = std::max(ixp, iyp);
    factor = std::max(factor, izp);
    poi_pre[maxlevel] = factor;

    int pbasis = dimx * dimy * dimz;
    int sbasis = (dimx + 2) * (dimy + 2) * (dimz + 2);

    /* Grab some memory for our multigrid structures */
    CalcType *mgrhsarr = new CalcType[std::max(sbasis, 512)];
    CalcType *mglhsarr = new CalcType[std::max(2*sbasis, 512)];
    CalcType *work = new CalcType[std::max(4*sbasis, 512)];
    CalcType *sg_res = new CalcType[std::max(2*sbasis, 512)];

    float *sg_res_f = (float *)sg_res; 
    float *mglhsarr_f = (float *)mglhsarr;
    float *work_f = (float *)work;

    CPP_app_cir_driver (L, T, rho, mgrhsarr, dimx, dimy, dimz, APP_CI_SIXTH);

    its = 0;
    while ( ((its < max_sweeps) && (residual > rms_target))  || (its < min_sweeps))
    {

        // Check if convergence is stalling at this level
        if(residual > 0.95*last_residual) break;
        last_residual = residual;

        /* Mehrstallen smoothings */
        for (cycles = 0; cycles < nits; cycles++)
        {

	    /*At the end of this loop, laplacian operator is reapplied to evalute the residual. Therefore,
	     * there is no need to reapply it, when this loop is called second, third, etc time. */
            if ( (cycles) || (!its))
            {
                /* Apply operator */
                diag = CPP_app_cil_driver (L, T, vhartree, mglhsarr, dimx, dimy, dimz, gridhx, gridhy, gridhz, APP_CI_SIXTH);
                diag = -1.0 / diag;
            }

            /* Pre and Post smoothings and multigrid */
            if (cycles == global_presweeps)
            {
                // Generate single precision residual vector and transfer into smoothing grid 
                for(int idx=0;idx < pbasis;idx++) work_f[idx] = (float)(mgrhsarr[idx] - mglhsarr[idx]);
                CPP_pack_ptos (sg_res_f, work_f, dimx, dimy, dimz);
                MG.mgrid_solv_pois (mglhsarr_f, sg_res_f, work_f,
                            dimx, dimy, dimz,
                            gridhx, gridhy, gridhz,
                            level, maxlevel, poi_pre,
                            poi_post, mu_cycles[level], coarse_step,
                            G->get_NX_GRID(density), G->get_NY_GRID(density), G->get_NZ_GRID(density),
                            G->get_PX_OFFSET(density), G->get_PY_OFFSET(density), G->get_PZ_OFFSET(density),
                            G->get_PX0_GRID(density), G->get_PY0_GRID(density), G->get_PZ0_GRID(density), boundaryflag);

                /* Transfer solution back to double array */
                for(int idx=0;idx < sbasis;idx++) work[idx] = (double)mglhsarr_f[idx];
                int t1 = 1.0;
                CPP_pack_stop_axpy (work, vhartree, t1, dimx, dimy, dimz);
            }
            else
            {
                /* Update vh */
                t1 = - global_step * diag;
                for(int i = 0;i < pbasis;i++) vhartree[i] = vhartree[i] + t1 * (mgrhsarr[i] - mglhsarr[i]);
            }                   /* end if */

        }                       /* end for */

        /*Get residual*/
        diag = CPP_app_cil_driver (L, T, vhartree, mglhsarr, dimx, dimy, dimz, gridhx, gridhy, gridhz, APP_CI_SIXTH);
        diag = -1.0 / diag;

        /* Compute residual */
        residual = 0.0;
        for (idx = 0; idx < pbasis; idx++)
        {
            residual += (double)(mgrhsarr[idx] - mglhsarr[idx])*(mgrhsarr[idx] - mglhsarr[idx]);
        } 

        residual = sqrt (RmgSumAll(residual, T->get_MPI_comm()) / (double)global_basis);
        //if(G->get_rank() == 0)printf("Hartree residual:   level=%d    sweep=%d    residual=%14.6e\n",level, its, residual);
        its ++;

    }   // end while


    t1 = - global_step * diag;
    for(int i = 0;i < pbasis;i++) vhartree[i] = vhartree[i] + t1 * (mgrhsarr[i] - mglhsarr[i]);


    if((G->get_rank() == 0) && print_status)
        std::cout << "\nget_vh: executed " << its << " sweeps, residual is " << residual << std::endl;


    if ((boundaryflag == PERIODIC) && setzero)
    {
        /* Evaluate the average potential */
        vavgcor = 0.0;
        for (idx = 0; idx < pbasis; idx++) vavgcor += (double)vhartree[idx];

        vavgcor =  RmgSumAll(vavgcor, T->get_MPI_comm());
        t1 = (double) global_basis;
        vavgcor = vavgcor / t1;

        /* Make sure that the average value is zero */
        for (idx = 0; idx < pbasis; idx++) vhartree[idx] -= (CalcType)vavgcor;
    }                   /* end if */


    /* Release our memory */
    delete [] sg_res;
    delete [] work;
    delete [] mglhsarr;
    delete [] mgrhsarr;

    return residual;

} // end coarse_vh


