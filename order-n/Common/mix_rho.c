/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
/****f* QMD-MGDFT/mix_rho.c *****
 * NAME
 *   Ab initio real space code with multigrid acceleration
 *   Quantum molecular dynamics package.
 *   Version: 2.1.5
 * COPYRIGHT
 *   Copyright (C) 2002  Wenchang Lu, Jerzy Bernholc
 * FUNCTION
 *   void mix_rho(REAL *rho, REAL *rho_old, int steps, int mode)
 *   This module is used to mix the charge density
 * INPUTS
 *   
 * OUTPUT
 *  
 * PARENTS
 *  
 * CHILDREN
 *  
 * SOURCE
 */


#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "md.h"

static REAL *rho1, *rho_res;

/* This function returns a pointer to a block of memory of size nelem. */
void mix_rho(REAL * rho_out, REAL * rho_in, REAL mix, int steps, int mode)
{

    REAL tem11, tem12, tem22;
    REAL alpha1, alpha2;
    int idx, ione = 1;

    double *work1, *work2;

    my_malloc_init( work1, P0_BASIS, REAL );
    my_malloc_init( work2, 4 * P0_BASIS, REAL );



    if (mode == 1)
    {                           /*   linear mixing   */

        for (idx = 0; idx < P0_BASIS; idx++)
        {
            rho_out[idx] = rho_out[idx] - rho_in[idx];
            work1[idx] = 0.1 * rho_out[idx];
        }

/*	precond_rho(rho_out, work1, work2); 
 */

        for (idx = 0; idx < P0_BASIS; idx++)
        {
            rho_out[idx] = rho_in[idx] + mix * work1[idx];
        }



    }                           /* endif steps ==0 */

    if (mode == 2)
    {                           /* Pulay mixing with one previous step  */

        if (steps <= 1)
        {
            my_malloc_init( rho1, P0_BASIS, REAL );
            my_malloc_init( rho_res, P0_BASIS, REAL );

            for (idx = 0; idx < P0_BASIS; idx++)
            {
                rho1[idx] = rho_in[idx];
                rho_res[idx] = rho_out[idx] - rho_in[idx];
                rho_out[idx] = mix * rho_out[idx] + (1.0 - mix) * rho_in[idx];
            }

        }                       /* if(steps ==0 ) */
        else
        {                       /* mode == 2, steps >0) */

            for (idx = 0; idx < P0_BASIS; idx++)
            {
                rho_out[idx] = rho_out[idx] - rho_in[idx];      /* new residual */
            }

            tem22 = sdot(&P0_BASIS, rho_out, &ione, rho_out, &ione);
            tem12 = sdot(&P0_BASIS, rho_out, &ione, rho_res, &ione);
            tem11 = sdot(&P0_BASIS, rho_res, &ione, rho_res, &ione);
            tem22 = real_sum_all(tem22);
            tem12 = real_sum_all(tem12);
            tem11 = real_sum_all(tem11);

            alpha1 = (tem11 - tem12) / (tem22 - 2.0 * tem12 + tem11);
            alpha2 = (tem22 - tem12) / (tem22 - 2.0 * tem12 + tem11);

            if (pct.thispe == 0)
                printf("\n mixing rho parameter %f %f", alpha1, alpha2);

            for (idx = 0; idx < P0_BASIS; idx++)
            {
                tem11 = (alpha1 * rho_in[idx] + alpha2 * rho1[idx])
                    + (alpha1 * rho_out[idx] + alpha2 * rho_res[idx]) * mix;
                rho1[idx] = rho_in[idx];
                rho_res[idx] = rho_out[idx];
                rho_out[idx] = tem11;

            }

        }                       /*  endif steps ==0 */
    }                           /* endif mode ==2 */

    if (mode == 3)
    {
        if (steps == 0)
        {
            my_malloc_init( rho1, P0_BASIS, REAL );
            my_malloc_init( rho_res, P0_BASIS, REAL );

            for (idx = 0; idx < P0_BASIS; idx++)
            {
                rho1[idx] = rho_in[idx];
                rho_res[idx] = rho_out[idx];
                rho_out[idx] = mix * rho_out[idx] + (1.0 - mix) * rho_in[idx];
            }

        }
        else
        {

            for (idx = 0; idx < P0_BASIS; idx++)
            {
                rho1[idx] = 0.9 * rho1[idx] + 0.1 * rho_in[idx];
                rho_res[idx] = 0.9 * rho_res[idx] + 0.1 * rho_out[idx];
                rho_out[idx] = mix * rho1[idx] + (1.0 - mix) * rho_res[idx];
            }
        }
    }

}
