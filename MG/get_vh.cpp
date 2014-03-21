/************************** SVN Revision Information **************************
 **    $Id: get_vh.c 2174 2014-03-03 15:34:10Z ebriggs $    **
******************************************************************************/

/****f* QMD-MGDFT/get_vh.c *****
 * NAME
 *   Ab initio real space code with multigrid acceleration
 *   Quantum molecular dynamics package.
 *   Version: 2.1.5
 * COPYRIGHT
 *   Copyright (C) 1995  Emil Briggs
 *   Copyright (C) 1998  Emil Briggs, Charles Brabec, Mark Wensell, 
 *                       Dan Sullivan, Chris Rapcewicz, Jerzy Bernholc
 *   Copyright (C) 2001  Emil Briggs, Wenchang Lu,
 *                       Marco Buongiorno Nardelli,Charles Brabec, 
 *                       Mark Wensell,Dan Sullivan, Chris Rapcewicz,
 *                       Jerzy Bernholc
 * FUNCTION
 *   void get_vh(rmg_double_t *rho, rmg_double_t *rhoc, rmg_double_t *vh_eig, int sweeps, int maxlevel)
 *   Iterative procedure to obtain the hartree potential.
 *   Uses Mehrstallen finite differencing with multigrid accelerations.
 *   The multigrid scheme uses a standard W-cycle with Jacobi relaxation.
 *   (Underrelaxed on the coarser grids.)
 * INPUTS
 *   rho: total charge density
 *   rhoc: compensating charge density
 *   sweeps:  number of iterations 
 *   maxlevel:  maximum multigrid level
 * OUTPUT
 *   vh_eig: Hartree potential
 * PARENTS
 *   init.c scf.c
 * CHILDREN
 *   getpoi_bc.c pack_vhstod.c pack_ptos.c app_cir.c 
 *   set_bc.c app_cil.c mgrid_solv.c pack_vhdtos.c
 * SOURCE
 */



#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#include "TradeImages.h"
#include "FiniteDiff.h"
#include "Mgrid.h"
#include "BlasWrappers.h"
#include "FineGrid.h"
#include "auxiliary.h"
#include "const.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "common_prototypes.h"
#include "common_prototypes1.h"
#include "rmg_alloc.h"
#include "rmg_error.h"
#include "RmgTimer.h"


/* Pre and post smoothings on each level */
static int poi_pre[5] = { 0, 3, 3, 3, 3 };
static int poi_post[5] = { 0, 3, 3, 3, 3 };


void CPP_get_vh (rmg_double_t * rho, rmg_double_t * rhoc, rmg_double_t * vh_eig, int min_sweeps, int max_sweeps, int maxlevel, rmg_double_t rms_target, int boundaryflag)
{

    int idx, its;
    rmg_double_t t1, vavgcor, diag;
    rmg_double_t *mgrhsarr, *mglhsarr, *mgresarr, *work;
    rmg_double_t *sg_rho, *sg_vh, *sg_res, *nrho,  residual = 100.0;
    int incx = 1, cycles;
    double k_vh;
    BaseGrid G;
    FineGrid FG(2);
    Lattice L;
    Mgrid MG;
    TradeImages T;

    int dimx = FG.get_PE_GRIDX(), dimy = FG.get_PE_GRIDY(), dimz = FG.get_PE_GRIDZ();

    /*Taken from ON code, seems to help a lot with convergence*/
    poi_pre[maxlevel] = ct.poi_parm.coarsest_steps;

    k_vh = 0.0;

    int nits = ct.poi_parm.gl_pre + ct.poi_parm.gl_pst;
    int pbasis = dimx * dimy * dimz;
    int sbasis = (dimx + 2) * (dimy + 2) * (dimz + 2);


    /* Grab some memory for our multigrid structures */
    mgrhsarr = new rmg_double_t[sbasis];
    mglhsarr = new  rmg_double_t[sbasis];
    mgresarr = new  rmg_double_t[sbasis];
    work = new  rmg_double_t[4*sbasis];
    sg_rho = new  rmg_double_t[sbasis];
    sg_vh = new  rmg_double_t[sbasis];
    sg_res = new  rmg_double_t[sbasis];
    nrho = new  rmg_double_t[sbasis];

    /* Subtract off compensating charges from rho */
    for (idx = 0; idx < pbasis; idx++)
        work[idx] = rho[idx] - rhoc[idx];

    for (idx = 0; idx < sbasis; idx++)
        nrho[idx] = 0.0;
    CPP_pack_stod (work, nrho, G.get_FPX0_GRID(), G.get_FPY0_GRID(), G.get_FPZ0_GRID(), boundaryflag);

    CPP_app_cir_driver<double> (nrho, mgrhsarr, dimx, dimy, dimz, APP_CI_FOURTH);



    /* Multiply through by 4PI */
    t1 = -FOUR * PI;
    for(idx = 0;idx < pbasis;idx++) mgrhsarr[idx] = t1 * mgrhsarr[idx];

    its = 0;

    while ( ((its < max_sweeps) && (residual > rms_target))  || (its < min_sweeps))
    {


        /* Mehrstallen smoothings */
        for (cycles = 0; cycles < nits; cycles++)
        {


	    /*At the end of this force loop, laplacian operator is reapplied to evalute the residual. Therefore,
	     * when there is no need to apply it, when this loop is called second, third, etc time. */
            if ( (cycles) || (!its))
            {

                /* Apply operator */
                diag = CPP_app_cil_driver (ct.vh_ext, mglhsarr, dimx, dimy, dimz,
                            L.hxxgrid, L.hyygrid, L.hzzgrid, APP_CI_FOURTH);
                diag = -1.0 / diag;

                /* Generate residual vector */
                for (idx = 0; idx < pbasis; idx++)
                {

                    mgresarr[idx] = mgrhsarr[idx] - mglhsarr[idx];

                }                   /* end for */

            }

            /* Pre and Post smoothings and multigrid */
            if (cycles == ct.poi_parm.gl_pre)
            {

                /* Transfer res into smoothing grid */
                CPP_pack_ptos<double> (sg_res, mgresarr, dimx, dimy, dimz);

                MG.mgrid_solv<double> (mglhsarr, sg_res, work,
                            dimx, dimy, dimz, L.hxxgrid,
                            L.hyygrid, L.hzzgrid,
                            0, G.get_neighbors(), ct.poi_parm.levels, poi_pre,
                            poi_post, ct.poi_parm.mucycles, ct.poi_parm.sb_step, k_vh,
                            G.FG_NX*G.get_NX_GRID(), G.FG_NY*G.get_NY_GRID(), G.FG_NZ*G.get_NZ_GRID(),
                            G.get_FPX_OFFSET(), G.get_FPY_OFFSET(), G.get_FPZ_OFFSET(),
                            G.get_FPX0_GRID(), G.get_FPY0_GRID(), G.get_FPZ0_GRID(), ct.boundaryflag);


                /* Transfer solution back to mgresarr array */
                CPP_pack_stop<double> (mglhsarr, mgresarr, dimx, dimy, dimz);

                /* Update vh */
                t1 = ONE;
                QMD_axpy (pbasis, t1, mgresarr, incx, ct.vh_ext, incx);

            }
            else
            {

                /* Update vh */
                t1 = -ct.poi_parm.gl_step * diag;
                QMD_axpy (pbasis, t1, mgresarr, incx, ct.vh_ext, incx);

            }                   /* end if */

            if (ct.boundaryflag == PERIODIC)
            {

                /* Evaluate the average potential */
                vavgcor = 0.0;
                for (idx = 0; idx < pbasis; idx++)
                    vavgcor += ct.vh_ext[idx];

                vavgcor = real_sum_all (vavgcor, T.get_MPI_comm());
                t1 = (rmg_double_t) ct.psi_fnbasis;
                vavgcor = vavgcor / t1;


                /* Make sure that the average value is zero */
                for (idx = 0; idx < pbasis; idx++)
                    ct.vh_ext[idx] -= vavgcor;

            }                   /* end if */

        }                       /* end for */
            
        /*Get residual*/
        diag = CPP_app_cil_driver<double> (ct.vh_ext, mglhsarr, dimx, dimy, dimz,
                            L.hxxgrid, L.hyygrid, L.hzzgrid, APP_CI_FOURTH);
        diag = -1.0 / diag;
        residual = 0.0;

        /* Generate residual vector */
        for (idx = 0; idx < pbasis; idx++)
        {

            mgresarr[idx] = mgrhsarr[idx] - mglhsarr[idx];
            residual += mgresarr[idx] * mgresarr[idx];

        }                   /* end for */

        residual = sqrt (real_sum_all(residual, pct.grid_comm) / ct.psi_fnbasis);

       
	//printf("\n get_vh sweep %3d, rms residual is %10.5e", its, residual);

	    
        its ++;
    }                           /* end for */

//    printf ("\n");
//    progress_tag ();
//    printf ("Executed %3d sweeps, residual is %15.8e, rms is %15.8e\n", its, residual, ct.rms);


    /* Pack the portion of the hartree potential used by the wavefunctions
     * back into the wavefunction hartree array. */
    CPP_pack_dtos (vh_eig, ct.vh_ext, G.get_FPX0_GRID(), G.get_FPY0_GRID(), G.get_FPZ0_GRID(), boundaryflag);

    /* Release our memory */
    delete [] nrho;
    delete [] sg_res;
    delete [] sg_vh;
    delete [] sg_rho;
    delete [] work;
    delete [] mgresarr;
    delete [] mglhsarr;
    delete [] mgrhsarr;

} // end CPP_get_vh


extern "C" void get_vh (rmg_double_t * rho, rmg_double_t * rhoc, rmg_double_t * vh_eig, int min_sweeps, int max_sweeps, int maxlevel, rmg_double_t rms_target, int boundaryflag) 
{
    CPP_get_vh (rho, rhoc, vh_eig, min_sweeps, max_sweeps, maxlevel, rms_target, boundaryflag);
}

/******/
