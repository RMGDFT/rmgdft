/************************** SVN Revision Information **************************
 **    $Id$    **
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
 *   void get_vh(double *rho, double *rhoc, double *vh_eig, int sweeps, int maxlevel)
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



#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
//#include "main.h"

#include "negf_prototypes.h"
#include "const.h"
#include "params.h"
#include "rmg_alloc.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "common_prototypes.h"
#include "common_prototypes1.h"
#include "macros.h"
#include "FiniteDiff.h"
#include "transition.h"




//#include "init_var.h"
//#include "LCR.h"
#include "twoParts.h"


/* Pre and post smoothings on each level */
static int poi_pre[5] = { 0, 3, 3, 3, 3 };
static int poi_post[5] = { 0, 3, 3, 3, 3 };


void get_vh_negf (double * rho, double * rhoc, double * vh_eig, int min_sweeps, int max_sweeps, int maxlevel, double rms_target)
{

    int idx, its, nits, sbasis, pbasis;
    double t1, diag=1.0;
    double *mgrhsarr, *mglhsarr, *mgresarr, *work;
    double *sg_rho, *sg_vh, *sg_res, *nrho,  residual = 100.0;
    int incx = 1, cycles;
    double k_vh;

    /*Taken from ON code, seems to help a lot with convergence*/
    poi_pre[maxlevel] = ct.poi_parm.coarsest_steps;

    k_vh = 0.0;

    nits = ct.poi_parm.gl_pre + ct.poi_parm.gl_pst;
    sbasis = (ct.vh_pxgrid + 2) * (ct.vh_pygrid + 2) * (ct.vh_pzgrid + 2);
    pbasis = ct.vh_pxgrid * ct.vh_pygrid * ct.vh_pzgrid;


    /* Grab some memory for our multigrid structures */
    my_malloc (mgrhsarr, sbasis, double);
    my_malloc (mglhsarr, sbasis, double);
    my_malloc (mgresarr, sbasis, double);
    my_malloc (work, 4*sbasis, double);
    my_malloc (sg_rho, sbasis, double);
    my_malloc (sg_vh, sbasis, double);
    my_malloc (sg_res, sbasis, double);
    my_malloc (nrho, sbasis, double);

    /* Subtract off compensating charges from rho */
    for (idx = 0; idx < pbasis; idx++)
        work[idx] = rho[idx] - rhoc[idx];

    for (idx = 0; idx < sbasis; idx++)
        nrho[idx] = 0.0;
    pack_vhstod (work, nrho, get_FPX0_GRID(), get_FPY0_GRID(), get_FPZ0_GRID(), ct.boundaryflag);

    CPP_app_cir_driver (&Rmg_L, Rmg_T, nrho, mgrhsarr, ct.vh_pxgrid, ct.vh_pygrid, ct.vh_pzgrid, APP_CI_FOURTH);


    /* Multiply through by 4PI */
    t1 = -FOUR * PI;
    QMD_dscal (pbasis, t1, mgrhsarr, incx);


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
                diag = CPP_app_cil_driver (&Rmg_L, Rmg_T, ct.vh_ext, mglhsarr, ct.vh_pxgrid, ct.vh_pygrid, ct.vh_pzgrid,
                            get_hxxgrid(), get_hyygrid(), get_hzzgrid(), APP_CI_FOURTH);
                diag = -1.0 / diag;

                /* Generate residual vector */
#pragma omp parallel for schedule(static, 1024)
                for (idx = 0; idx < pbasis; idx++)
                {

                    mgresarr[idx] = mgrhsarr[idx] - mglhsarr[idx];

                }                   /* end for */

            }

         /*  Fix Hartree in some region  */
            confine (mgresarr, ct.vh_pxgrid, ct.vh_pygrid, ct.vh_pzgrid, potentialCompass, 0);

            /* Pre and Post smoothings and multigrid */
            if (cycles == ct.poi_parm.gl_pre)
            {

                /* Transfer res into smoothing grid */
                pack_ptos (sg_res, mgresarr, ct.vh_pxgrid, ct.vh_pygrid, ct.vh_pzgrid);

                 
                mgrid_solv_negf (mglhsarr, sg_res, work,
                            ct.vh_pxgrid, ct.vh_pygrid, ct.vh_pzgrid, get_hxxgrid(),
                            get_hyygrid(), get_hzzgrid(),
                            0, get_neighbors(), ct.poi_parm.levels, poi_pre,
                            poi_post, ct.poi_parm.mucycles, ct.poi_parm.sb_step, k_vh,
                            get_FG_RATIO()*get_NX_GRID(), get_FG_RATIO()*get_NY_GRID(), get_FG_RATIO()*get_NZ_GRID(),
                            get_FPX_OFFSET(), get_FPY_OFFSET(), get_FPZ_OFFSET(),
                            get_FPX0_GRID(), get_FPY0_GRID(), get_FPZ0_GRID());



                /* Transfer solution back to mgresarr array */
                pack_stop (mglhsarr, mgresarr, ct.vh_pxgrid, ct.vh_pygrid, ct.vh_pzgrid);

                /*  Fix Hartree in some region  */

                confine (mgresarr, ct.vh_pxgrid, ct.vh_pygrid, ct.vh_pzgrid, potentialCompass, 0);

                /* Update vh */
                t1 = ONE;
                QMD_daxpy (pbasis, t1, mgresarr, incx, ct.vh_ext, incx);

            }
            else
            {

                /* Update vh */
                t1 = -ct.poi_parm.gl_step * diag;
                QMD_daxpy (pbasis, t1, mgresarr, incx, ct.vh_ext, incx);

            }                   /* end if */



        }                       /* end for */

        /*Get residual*/
        diag = CPP_app_cil_driver (&Rmg_L, Rmg_T, ct.vh_ext, mglhsarr, ct.vh_pxgrid, ct.vh_pygrid, ct.vh_pzgrid,
                            get_hxxgrid(), get_hyygrid(), get_hzzgrid(), APP_CI_FOURTH);
        diag = -1.0 / diag;

        /* Generate residual vector */
        for (idx = 0; idx < pbasis; idx++)
        {

            mgresarr[idx] = mgrhsarr[idx] - mglhsarr[idx];

        }                   /* end for */

        confine (mgresarr, ct.vh_pxgrid, ct.vh_pygrid, ct.vh_pzgrid, potentialCompass, 0);

        residual = 0.0;
        for (idx = 0; idx < pbasis; idx++)
           residual += mgresarr[idx] * mgresarr[idx];


        residual = sqrt (real_sum_all(residual, pct.grid_comm) / ct.psi_fnbasis);


        //printf("\n get_vh sweep %3d, rms residual is %10.5e", its, residual);


        its ++;
    }                           /* end for */

    printf ("\n");
    progress_tag ();
    printf ("Executed %3d sweeps, residual is %15.8e, rms is %15.8e\n", its, residual, ct.rms);


    /* Pack the portion of the hartree potential used by the wavefunctions
     * back into the wavefunction hartree array. */
    pack_vhdtos (vh_eig, ct.vh_ext, get_FPX0_GRID(), get_FPY0_GRID(), get_FPZ0_GRID(), ct.boundaryflag);

    /* Release our memory */
    my_free (nrho);
    my_free (sg_res);
    my_free (sg_vh);
    my_free (sg_rho);
    my_free (work);
    my_free (mgresarr);
    my_free (mglhsarr);
    my_free (mgrhsarr);


}                               /* end get_vh */

