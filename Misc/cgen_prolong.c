/*
 *
 * Copyright 2014 The RMG Project Developers. See the COPYRIGHT file 
 * at the top-level directory of this distribution or in the current
 * directory.
 * 
 * This file is part of RMG. 
 * RMG is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * any later version.
 *
 * RMG is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
*/

/****f* QMD-MGDFT/mg_restrict_6.c *****
 * NAME
 *   Ab initio real space code with multigrid acceleration
 *   Quantum molecular dynamics package.
 *   Version: 2.1.5

 * FUNCTION
 *   void cgen_prolong(double  coef[], double fraction, int order)
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





void cgen_prolong (double coef[], double fraction, int order)
{

    int ix, iy;
    double A[order * order];
    double b[order];
    double d[order];
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
    dgesv (&order, &ione, A, &order, ipvt, b, &order, &info);



    for (ix = 0; ix < order; ix++)
    {

        coef[ix] = b[ix];

    }


}                               /* end cgen-prolong */


/******/
