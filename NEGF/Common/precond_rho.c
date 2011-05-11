/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include "md.h"



void precond_rho (double *res)
{
    int cycles, ione = 1, nits, sbasis;
    double diag, t1, time1, time2, d1, one = 1.;
    double *sg_rho, *rho_res, *work1, *work2;
    int idx;

    sbasis = (FPX0_GRID + 2) * (FPY0_GRID + 2) * (FPZ0_GRID + 2);
    my_malloc_init( sg_rho, 2 * sbasis, REAL );
    rho_res = sg_rho + sbasis;
    my_malloc_init( work1, sbasis, REAL );
    my_malloc_init( work2, 4 * sbasis, REAL );

    /* Pre and post smoothings on each level */
    int eig_pre[6] = { 0, 2, 2, 2, 2, 2 };
    int eig_post[6] = { 0, 2, 2, 2, 2, 2 };


    time1 = my_crtc ();

    for (idx = 0; idx < FP0_BASIS; idx++)
    {
        work1[idx] = 0.05 * res[idx];
    }

    nits = ct.eig_parm.gl_pre + ct.eig_parm.gl_pst;


    /* Smoothing cycles */
    for (cycles = 0; cycles <= nits; cycles++)
    {


        pack_ptos (rho_res, work1, FPX0_GRID, FPY0_GRID, FPZ0_GRID);
        trade_images (rho_res, FPX0_GRID, FPY0_GRID, FPZ0_GRID, &pct.neighbors[0]);

        diag =
            app_cil (rho_res, work2, FPX0_GRID, FPY0_GRID, FPZ0_GRID, ct.hxxgrid, ct.hyygrid,
                     ct.hzzgrid);
        diag = -1.0 / diag;

        saxpy (&FP0_BASIS, &one, res, &ione, work2, &ione);

        /* Now either smooth the wavefunction or do a multigrid cycle */
        if (cycles == ct.eig_parm.gl_pre)
        {

            /* Pack the residual data into multigrid array */
            pack_ptos (sg_rho, work2, FPX0_GRID, FPY0_GRID, FPZ0_GRID);

            trade_images (sg_rho, FPX0_GRID, FPY0_GRID, FPZ0_GRID, &pct.neighbors[0]);

            /* Do multigrid step with solution in sg_twovpsi */


            mgrid_solv (rho_res, sg_rho, work2, FPX0_GRID, FPY0_GRID, FPZ0_GRID,
                        ct.hxxgrid, ct.hyygrid, ct.hzzgrid,
                        0, pct.neighbors, ct.eig_parm.levels, eig_pre, eig_post, 1, 0, 0, 0);


            pack_stop (rho_res, work2, FPX0_GRID, FPY0_GRID, FPZ0_GRID);


            t1 = -1.;

        }
        else
        {

            t1 = diag;

        }                       /* end if cycles == ct.eig_parm.gl_pre */

        /* Update correction for wavefuntion */
        saxpy (&FP0_BASIS, &t1, work2, &ione, work1, &ione);


    }                           /* end for Smoothing cycles */


    for (idx = 0; idx < FP0_BASIS; idx++)
    {
        res[idx] = work1[idx];
    }

    my_free(sg_rho);
    my_free(work1);
    my_free(work2);
    time2 = my_crtc ();
    d1 = time2 - time1;
    rmg_timings (PRECOND_TIME, d1);


}
