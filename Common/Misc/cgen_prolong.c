/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

/****f* QMD-MGDFT/mg_restrict_6.c *****
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
 *   void cgen_prolong(REAL  coef[], REAL fraction, int order)
 * INPUTS
 *   int order,(total number of neighbour points we need)
 *   double fraction,(the relative position of the interpolated point inside two nearest neighbour)
 *   fraction in [0,1] 
 *   
 * OUTPUT
 *   c[order] the coefficients for 1-D interpolation.
 * PARENTS
 *   mg_prolong_MAX10
 * CHILDREN
 *   nothing 
 * SOURCE
 */



#include "main.h"
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>





void cgen_prolong (REAL coef[], REAL fraction, int order)
{

    int ix, iy;
    REAL A[order * order];
    REAL b[order];
    REAL d[order];
    int ipvt[order];
    int info;
    int ione = 1;


    /*initialize A and b to be zero */


    for (ix = 0; ix < order; ix++)
    {

        b[ix] = 0.0;

        for (iy = 0; iy < order; iy++)
        {

            A[ix * order + iy] = 0.0;

        }                       /* end for */

    }                           /* end for */


    /*filling A , b and d (d is distance from different coarse grids to the interpolated pt, 
       coarse grid spacing is normalized to be ONE) */

    b[0] = 1.0;

    for (iy = 0; iy < order; iy++)
    {


        d[iy] = iy + 1.0 - fraction - (double) order / 2;

        for (ix = 0; ix < order; ix++)

        {

            A[iy * order + ix] = pow (d[iy], ix);
            /*  printf("  A[%d][%d]= %f  ", ix, iy, A[ix][iy]); */

            /*  here we flip ix and iy, just to transpose matrix A since C code treat array in row fashion */

        }                       /* end for */


    }                           /* end for */


    /*  solving Ac=b for c using  b = A^(-1) * b  */
    sgesv (&order, &ione, A, &order, ipvt, b, &order, &info);



    for (ix = 0; ix < order; ix++)
    {

        coef[ix] = b[ix];

    }


}                               /* end cgen-prolong */


/******/
