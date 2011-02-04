/************************** SVN Revision Information **************************
 **    $Id: radiff.c 1242 2011-02-02 18:55:23Z luw $    **
******************************************************************************/
 
/****f* QMD-MGDFT/radiff.c *****
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
 *   void radiff(double *f, double *df, double *r, int n, double al)
 *   Radial differentiation functions for data defined on a log grid.
 *   We assume that the derivative at the large r-end is zero.
 *   If the log mesh parameter is zero then the grid is assumed linear.
 * INPUTS
 *   f:  function being differentiated
 *   r:  array of r-values in log-mesh
 *   n:  number of points in radial mesh
 *   al: log mesh parameter
 * OUTPUT
 *   df:   array of derivative values
 * PARENTS
 *   init_kbr.c
 * CHILDREN
 *   nothing
 * SOURCE
 */



#include "md.h"
#include <float.h>
#include <math.h>


void radiff (double *f, double *df, double *r, int n, double al)
{

    int idx;
    double dr;

    /* Check for linear grid case */
    if (al == 0.0)
    {

        dr = r[2] - r[1];
        for (idx = 2; idx < n - 2; idx++)
        {

            df[idx] = f[idx - 2] - 8.0 * f[idx - 1] + 8.0 * f[idx + 1] - f[idx + 2];

            df[idx] = df[idx] / (12.0 * dr);

        }                       /* end for */

        df[n - 2] = df[n - 1] = 0.0;

        /* Take care of small r endpoint */
        df[1] = -f[0] / 3.0 - f[1] / 2.0 + f[2] - f[3] / 6.0;

        df[0] = -3.0 * f[0] + 4.0 * f[1] - f[2];


    }
    else
    {

        for (idx = 2; idx < n - 2; idx++)
        {

            df[idx] = f[idx - 2] - 8.0 * f[idx - 1] + 8.0 * f[idx + 1] - f[idx + 2];
            df[idx] = df[idx] / (12.0 * al * r[idx]);

        }                       /* end for */

        df[n - 2] = df[n - 1] = 0.0;

        /* Take care of small r endpoint */
        df[1] = -f[0] / 3.0 - f[1] / 2.0 + f[2] - f[3] / 6.0;

        df[1] = df[1] / (al * r[1]);

        df[0] = -3.0 * f[0] + 4.0 * f[1] - f[2];
        df[0] = df[0] / (2.0 * al * r[0]);

    }                           /* end if */


}                               /* end radiff */



/******/
