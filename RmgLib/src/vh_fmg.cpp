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


double coarse_vh (BaseGrid *G, Lattice *L, TradeImages *T, double * rho, double *vhartree,
                 int min_sweeps, int max_sweeps, int maxlevel, 
                 int global_presweeps, int global_postsweeps, int mucycles, 
                 int dimx, int dimy, int dimz, int level,
                 double gridhx, double gridhy, double gridhz,
                 double rms_target, double global_step, double coarse_step, int boundaryflag, int density, bool print_status, bool setzero);

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
                 double rms_target, double global_step, double coarse_step, int boundaryflag, int density, bool print_status)
{

RmgTimer *RT0 = new RmgTimer("Hartree: init");
    int idx, its, cycles;
    double t1, vavgcor, diag=0.0;
    double residual = 100.0;
    Mgrid MG(L, T);
    int global_basis = G->get_GLOBAL_BASIS(density);

    if(maxlevel >= MAX_MG_LEVELS)
       rmg_error_handler(__FILE__, __LINE__, "Too many multigrid levels requested.");

    int dimx = G->get_PX0_GRID(density), dimy = G->get_PY0_GRID(density), dimz = G->get_PZ0_GRID(density);

    // Solve to a high degree of precision on the coarsest level
    int nits = global_presweeps + global_postsweeps;
    int pbasis = dimx * dimy * dimz;
    int sbasis = (dimx + 2) * (dimy + 2) * (dimz + 2);

    /* Grab some memory for our multigrid structures */
    double *mgrhsarr = new double[2*sbasis];
    double *mglhsarr = new  double[2*sbasis];
    double *mgresarr = new  double[sbasis];
    double *work = new  double[2*sbasis];
    double *sg_res = new  double[2*sbasis];

    int dx[MAX_MG_LEVELS];dx[0] = dimx;
    int dy[MAX_MG_LEVELS];dy[0] = dimy;
    int dz[MAX_MG_LEVELS];dz[0] = dimz;
    int ixoff, iyoff, izoff;
    double *mgrhsptr[MAX_MG_LEVELS];


    // Set up RHS on finest level
    /* Multiply through by 4PI */
    t1 = -4.0 * PI;
    for(idx = 0;idx < pbasis;idx++) mgrhsarr[idx] = t1 * rho[idx];

    // Restrict right hand side down to coarsest level
    size_t offset=pbasis;
    int dx2, dy2, dz2;
    mgrhsptr[0] = mgrhsarr;
    for(int level=1;level <= maxlevel;level++)
    {
        dx2 = MG.MG_SIZE (dx[level-1], level-1, G->get_NX_GRID(density), G->get_PX_OFFSET(density), dimx, &ixoff, boundaryflag);
        dy2 = MG.MG_SIZE (dy[level-1], level-1, G->get_NY_GRID(density), G->get_PY_OFFSET(density), dimy, &iyoff, boundaryflag);
        dz2 = MG.MG_SIZE (dz[level-1], level-1, G->get_NZ_GRID(density), G->get_PZ_OFFSET(density), dimz, &izoff, boundaryflag);

        CPP_pack_ptos (work, mgrhsptr[level-1], dx[level-1], dy[level-1], dz[level-1]);
        T->trade_images (work, dx[level-1], dy[level-1], dz[level-1], FULL_TRADE);

        MG.mg_restrict (work, sg_res, dx[level-1], dy[level-1], dz[level-1], dx2, dy2, dz2, ixoff, iyoff, izoff);
        mgrhsptr[level] = &mgrhsarr[offset];
        CPP_pack_stop (sg_res, mgrhsptr[level], dx2, dy2, dz2);
        offset += dx2*dy2*dz2;

        dx[level] = dx2;
        dy[level] = dy2;
        dz[level] = dz2;

    }
delete RT0;
RmgTimer *RT1 = new RmgTimer("Hartree: cg solve");
    // Now solve from coarse grid to fine grid
    for(int ix=0;ix < dx2*dy2*dz2;ix++) mglhsarr[ix] = 0.0;
    for(int level=maxlevel;level > 0;level--)
    {
        double lfactor = pow(2.0, (double)(level));
        coarse_vh (G, L, T, mgrhsptr[level], mglhsarr,
                 min_sweeps, max_sweeps, maxlevel,
                 global_presweeps, global_postsweeps, level,
                 dx[level], dy[level], dz[level], level,
                 G->get_hxgrid(density)*lfactor, G->get_hygrid(density)*lfactor, G->get_hzgrid(density)*lfactor,
                 1.0e-14, global_step, coarse_step, boundaryflag, density, false, false);
        CPP_pack_ptos (work, mglhsarr, dx[level], dy[level], dz[level]);
        T->trade_images (work, dx[level], dy[level], dz[level], FULL_TRADE);
        MG.mg_prolong_cubic (sg_res, work, dx[level], dy[level], dz[level], dx[level-1], dy[level-1], dz[level-1], ixoff, iyoff, izoff);
        CPP_pack_stop (sg_res, mglhsarr, dx[level-1], dy[level-1], dz[level-1]);
    }

delete RT1;
RmgTimer *RT2 = new RmgTimer("Hartree: fg solve");

    coarse_vh (G, L, T, mgrhsarr, mglhsarr,
             3, max_sweeps, maxlevel,
             global_presweeps, global_postsweeps, 3,
             dimx, dimy, dimz, 0,
             G->get_hxgrid(density), G->get_hygrid(density), G->get_hzgrid(density),
             1.0e-10, global_step, coarse_step, boundaryflag, density, false, true);

//printf("PPPP  %d  %d  %d\n",dx[0], dy[0], dz[0]);
    for(int idx=0;idx < pbasis;idx++) vhartree[idx] = mglhsarr[idx];

    /* Release our memory */
    delete [] sg_res;
    delete [] work;
    delete [] mgresarr;
    delete [] mglhsarr;
    delete [] mgrhsarr;

delete RT2;
    return residual;

} // vh_fmg





double coarse_vh (BaseGrid *G, Lattice *L, TradeImages *T, double * rho, double *vhartree,
                 int min_sweeps, int max_sweeps, int maxlevel, 
                 int global_presweeps, int global_postsweeps, int mucycles, 
                 int dimx, int dimy, int dimz, int level,
                 double gridhx, double gridhy, double gridhz,
                 double rms_target, double global_step, double coarse_step, int boundaryflag, int density, bool print_status, bool setzero)
{

    int idx, its, cycles;
    double t1, vavgcor, diag=0.0;
    double *mgrhsarr, *mglhsarr, *mgresarr, *work;
    double *sg_res, residual = 100.0;
    Mgrid MG(L, T);
    int global_basis = G->get_GLOBAL_BASIS(density) / pow(8.0, (double)level);

    /* Pre and post smoothings on each level */
    int poi_pre[MAX_MG_LEVELS] = { 0, 4, 4, 4, 4, 4, 4, 4};
    int poi_post[MAX_MG_LEVELS] = { 0, 2, 2, 2, 2, 2, 2, 2};

    int mu_cycles[MAX_MG_LEVELS] = {2, 2, 2, 2, 2, 2, 2, 2};
    if(maxlevel >= MAX_MG_LEVELS)
       rmg_error_handler(__FILE__, __LINE__, "Too many multigrid levels requested.");

    // Solve to a high degree of precision on the coarsest level
    poi_pre[maxlevel] = 10;
    int nits = global_presweeps + global_postsweeps;
    int pbasis = dimx * dimy * dimz;
    int sbasis = (dimx + 2) * (dimy + 2) * (dimz + 2);

    /* Grab some memory for our multigrid structures */
    mgrhsarr = new double[std::max(sbasis, 512)];
    mglhsarr = new  double[std::max(2*sbasis, 512)];
    mgresarr = new  double[std::max(sbasis, 512)];
    work = new  double[std::max(4*sbasis, 512)];
    sg_res = new  double[std::max(2*sbasis, 512)];

    float *sg_res_f = (float *)sg_res; 
    float *mglhsarr_f = (float *)mglhsarr;
    float *work_f = (float *)work;

    CPP_app_cir_driver<double> (L, T, rho, mgrhsarr, dimx, dimy, dimz, APP_CI_FOURTH);

    its = 0;
    while ( ((its < max_sweeps) && (residual > rms_target))  || (its < min_sweeps))
    {


        /* Mehrstallen smoothings */
        for (cycles = 0; cycles < nits; cycles++)
        {

	    /*At the end of this force loop, laplacian operator is reapplied to evalute the residual. Therefore,
	     * there is no need to reapply it, when this loop is called second, third, etc time. */
            if ( (cycles) || (!its))
            {

                /* Apply operator */
                diag = CPP_app_cil_driver (L, T, vhartree, mglhsarr, dimx, dimy, dimz,
                            gridhx, gridhy, gridhz, APP_CI_FOURTH);

                diag = -1.0 / diag;


            }

            /* Pre and Post smoothings and multigrid */
            if (cycles == global_presweeps)
            {

                /* Generate residual vector and transfer into smoothing grid */
                for(int idx=0;idx < pbasis;idx++) work_f[idx] = (float)(mgrhsarr[idx] - mglhsarr[idx]);
                CPP_pack_ptos<float> (sg_res_f, work_f, dimx, dimy, dimz);

                MG.mgrid_solv_pois<float> ((float *)mglhsarr, sg_res_f, (float *)work,
                            dimx, dimy, dimz,
                            gridhx, gridhy, gridhz,
                            level, maxlevel, poi_pre,
                            poi_post, mu_cycles[level], coarse_step,
                            G->get_NX_GRID(density), G->get_NY_GRID(density), G->get_NZ_GRID(density),
                            G->get_PX_OFFSET(density), G->get_PY_OFFSET(density), G->get_PZ_OFFSET(density),
                            G->get_PX0_GRID(density), G->get_PY0_GRID(density), G->get_PZ0_GRID(density), boundaryflag);


                /* Transfer solution back to mgresarr array */
                for(int idx=0;idx < sbasis;idx++) work[idx] = (double)mglhsarr_f[idx];
                int t1 = 1.0;
                CPP_pack_stop_axpy (work, vhartree, t1, dimx, dimy, dimz);


                /* Update vh */
                //for(int i = 0;i < pbasis;i++) vhartree[i] += mgresarr[i];

            }
            else
            {

                /* Update vh */
                t1 = - global_step * diag;
                for(int i = 0;i < pbasis;i++) vhartree[i] = vhartree[i] + t1 * (mgrhsarr[i] - mglhsarr[i]);

            }                   /* end if */

        }                       /* end for */

        if ((boundaryflag == PERIODIC) && setzero)
        {

            /* Evaluate the average potential */
            vavgcor = 0.0;
            for (idx = 0; idx < pbasis; idx++) vavgcor += vhartree[idx];

            vavgcor =  RmgSumAll(vavgcor, T->get_MPI_comm());
            t1 = (double) global_basis;
            vavgcor = vavgcor / t1;

            /* Make sure that the average value is zero */
            for (idx = 0; idx < pbasis; idx++) vhartree[idx] -= vavgcor;

        }                   /* end if */

        /*Get residual*/
        diag = CPP_app_cil_driver<double> (L, T, vhartree, mglhsarr, dimx, dimy, dimz,
                            gridhx, gridhy, gridhz, APP_CI_FOURTH);
        diag = -1.0 / diag;
        residual = 0.0;

        /* Generate residual vector */
        for (idx = 0; idx < pbasis; idx++)
        {

            mgresarr[idx] = mgrhsarr[idx] - mglhsarr[idx];
            residual += mgresarr[idx] * mgresarr[idx];

        }                   /* end for */


        residual = sqrt (RmgSumAll(residual, T->get_MPI_comm()) / global_basis);
//        if(G->get_rank() == 0) std::cout << "\n get_vh sweep " << its << " rms residual is " << residual << std::endl;
//if(G->get_rank() == 0)printf("LLLL  %d  %d  %14.6e\n",level, its, residual);
        its ++;
    }                           /* end for */


    if((G->get_rank() == 0) && print_status)
        std::cout << "\nget_vh: executed " << its << " sweeps, residual is " << residual << std::endl;


    /* Release our memory */
    delete [] sg_res;
    delete [] work;
    delete [] mgresarr;
    delete [] mglhsarr;
    delete [] mgrhsarr;

    return residual;

} // end coarse_vh


