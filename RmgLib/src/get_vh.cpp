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
double CPP_get_vh (BaseGrid *G, Lattice *L, TradeImages *T, double * rho, double *vhartree,
                 int min_sweeps, int max_sweeps, int maxlevel, 
                 int global_presweeps, int global_postsweeps, int mucycles, 
                 double rms_target, double global_step, double coarse_step, int boundaryflag, int density, bool print_status)
{

    int idx, its, cycles;
    double t1, vavgcor, diag=0.0;
    double *mgrhsarr, *mglhsarr, *mgresarr, *work;
    double *sg_res, residual = 100.0;
    Mgrid MG(L, T);

    int global_basis = G->get_GLOBAL_BASIS(density);

    /* Pre and post smoothings on each level */
    int poi_pre[MAX_MG_LEVELS] = { 0, 3, 3, 3, 3, 3, 3, 3 };
    int poi_post[MAX_MG_LEVELS] = { 0, 3, 3, 3, 3, 3, 3, 3 };

    if(maxlevel >= MAX_MG_LEVELS)
       rmg_error_handler(__FILE__, __LINE__, "Too many multigrid levels requested.");

    int dimx = G->get_PX0_GRID(density), dimy = G->get_PY0_GRID(density), dimz = G->get_PZ0_GRID(density);

    // Solve to a high degree of precision on the coarsest level
    poi_pre[maxlevel] = 20;
    int nits = global_presweeps + global_postsweeps;
    int pbasis = dimx * dimy * dimz;
    int sbasis = (dimx + 2) * (dimy + 2) * (dimz + 2);

    /* Grab some memory for our multigrid structures */
    mgrhsarr = new double[sbasis];
    mglhsarr = new  double[2*sbasis];
    mgresarr = new  double[sbasis];
    work = new  double[4*sbasis];
    sg_res = new  double[2*sbasis];


    CPP_app_cir_driver<double> (L, T, rho, mgrhsarr, dimx, dimy, dimz, APP_CI_FOURTH);

    /* Multiply through by 4PI */
    t1 = -4.0 * PI;
    for(idx = 0;idx < pbasis;idx++) mgrhsarr[idx] = t1 * mgrhsarr[idx];

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
                diag = CPP_app_cil_driver_threaded (L, T, vhartree, mglhsarr, dimx, dimy, dimz,
                            G->get_hxgrid(density), G->get_hygrid(density), G->get_hzgrid(density), APP_CI_FOURTH);

                diag = -1.0 / diag;

                /* Generate residual vector */
                for (idx = 0; idx < pbasis; idx++)
                {

                    mgresarr[idx] = mgrhsarr[idx] - mglhsarr[idx];

                }                   /* end for */

            }

            /* Pre and Post smoothings and multigrid */
            if (cycles == global_presweeps)
            {

                /* Transfer res into smoothing grid */
                CPP_pack_ptos<double> (sg_res, mgresarr, dimx, dimy, dimz);

                MG.mgrid_solv_pois<double> (mglhsarr, sg_res, work,
                            dimx, dimy, dimz,
                            G->get_hxgrid(density), G->get_hygrid(density), G->get_hzgrid(density),
                            0, G->get_neighbors(), maxlevel, poi_pre,
                            poi_post, mucycles, coarse_step,
                            G->get_NX_GRID(density), G->get_NY_GRID(density), G->get_NZ_GRID(density),
                            G->get_PX_OFFSET(density), G->get_PY_OFFSET(density), G->get_PZ_OFFSET(density),
                            G->get_PX0_GRID(density), G->get_PY0_GRID(density), G->get_PZ0_GRID(density), boundaryflag);


                /* Transfer solution back to mgresarr array */
                CPP_pack_stop<double> (mglhsarr, mgresarr, dimx, dimy, dimz);

                /* Update vh */
                for(int i = 0;i < pbasis;i++) vhartree[i] += mgresarr[i];

            }
            else
            {

                /* Update vh */
                t1 = - global_step * diag;
                for(int i = 0;i < pbasis;i++) vhartree[i] = vhartree[i] + t1 * mgresarr[i];

            }                   /* end if */

            if (boundaryflag == PERIODIC)
            {

                /* Evaluate the average potential */
                vavgcor = 0.0;
                for (idx = 0; idx < pbasis; idx++)
                    vavgcor += vhartree[idx];

                vavgcor =  RmgSumAll(vavgcor, T->get_MPI_comm());
                t1 = (double) global_basis;
                vavgcor = vavgcor / t1;


                /* Make sure that the average value is zero */
                for (idx = 0; idx < pbasis; idx++)
                    vhartree[idx] -= vavgcor;

            }                   /* end if */

        }                       /* end for */

        /*Get residual*/
        diag = CPP_app_cil_driver_threaded<double> (L, T, vhartree, mglhsarr, dimx, dimy, dimz,
                            G->get_hxgrid(density), G->get_hygrid(density), G->get_hzgrid(density), APP_CI_FOURTH);
        diag = -1.0 / diag;
        residual = 0.0;

        /* Generate residual vector */
        for (idx = 0; idx < pbasis; idx++)
        {

            mgresarr[idx] = mgrhsarr[idx] - mglhsarr[idx];
            residual += mgresarr[idx] * mgresarr[idx];

        }                   /* end for */

        residual = sqrt (RmgSumAll(residual, T->get_MPI_comm()) / global_basis);
        //if(G->get_rank() == 0) std::cout << "\n get_vh sweep " << its << " rms residual is " << residual << std::endl;

	    
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

} // end CPP_get_vh


