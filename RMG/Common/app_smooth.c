/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

/****f* QMD-MGDFT/app_smooth.c *****
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
 *   void app_smooth(S0_GRID *f, S0_GRID *work, REAL sfac)
 *   Smooths the function f by averaging the neighbors 
 * INPUTS
 *   f: array to be smoothed
 *   sfac: scale factor, it is useless now
 * OUTPUT
 *   work: smoothed results 
 * PARENTS
 *   mg_eig_state.c
 * CHILDREN
 *   nothing
 * SOURCE
 */



#include "main.h"
#include <float.h>
#include <math.h>


/**

 Smooths f and returns result in work
*/
void app_smooth (REAL * f, REAL * work, int dimx, int dimy, int dimz)
{

    int iz, ix, iy, incx, incy;
    int ixs, iys, ixms, ixps, iyms, iyps;

    REAL scale, ec, fc, crn, cc, temp;

    incy = dimz + 2;
    incx = (dimz + 2) * (dimy + 2);

    cc = 8.0;
    fc = 8.0;
    ec = 8.0;
    crn = 8.0;
    scale = 1.0 / 216.0;


    for (ix = 1; ix <= dimx; ix++)
    {

        ixs = ix * incx;
        ixms = (ix - 1) * incx;
        ixps = (ix + 1) * incx;

        for (iy = 1; iy <= dimy; iy++)
        {

            iys = iy * incy;
            iyms = (iy - 1) * incy;
            iyps = (iy + 1) * incy;

            for (iz = 1; iz <= dimz; iz++)
            {

                temp = cc * f[ixs + iys + iz] +
                    fc * (f[ixms + iys + iz] +
                          f[ixps + iys + iz] +
                          f[ixs + iyms + iz] +
                          f[ixs + iyps + iz] +
                          f[ixs + iys + (iz - 1)] + 
                          f[ixs + iys + (iz + 1)]);


                temp +=
                    ec * (f[ixms + iys + iz - 1] +
                          f[ixps + iys + iz - 1] +
                          f[ixs + iyms + iz - 1] +
                          f[ixs + iyps + iz - 1] +
                          f[ixms + iyms + iz] +
                          f[ixms + iyps + iz] +
                          f[ixps + iyms + iz] +
                          f[ixps + iyps + iz] +
                          f[ixms + iys + iz + 1] +
                          f[ixps + iys + iz + 1] +
                          f[ixs + iyms + iz + 1] + 
                          f[ixs + iyps + iz + 1]);


                temp +=
                    crn * (f[ixms + iyms + iz - 1] +
                           f[ixms + iyms + iz + 1] +
                           f[ixms + iyps + iz - 1] +
                           f[ixms + iyps + iz + 1] +
                           f[ixps + iyms + iz - 1] +
                           f[ixps + iyms + iz + 1] +
                           f[ixps + iyps + iz - 1] +
                           f[ixps + iyps + iz + 1]);

                work[ixs + iys + iz] = scale * temp;


            }                   /* end for */

        }                       /* end for */

    }                           /* end for */

}                               /* end app_smooth */


void app_smooth_f (rmg_float_t * f, rmg_float_t * work, int dimx, int dimy, int dimz)
{

    int iz, ix, iy, incx, incy;
    int ixs, iys, ixms, ixps, iyms, iyps;

    REAL scale, ec, fc, crn, cc, temp;

    incy = dimz + 2;
    incx = (dimz + 2) * (dimy + 2);

    cc = 8.0;
    fc = 8.0;
    ec = 8.0;
    crn = 8.0;
    scale = 1.0 / 216.0;


    for (ix = 1; ix <= dimx; ix++)
    {

        ixs = ix * incx;
        ixms = (ix - 1) * incx;
        ixps = (ix + 1) * incx;

        for (iy = 1; iy <= dimy; iy++)
        {

            iys = iy * incy;
            iyms = (iy - 1) * incy;
            iyps = (iy + 1) * incy;

            for (iz = 1; iz <= dimz; iz++)
            {

                temp = cc * f[ixs + iys + iz] +
                    fc * (f[ixms + iys + iz] +
                          f[ixps + iys + iz] +
                          f[ixs + iyms + iz] +
                          f[ixs + iyps + iz] +
                          f[ixs + iys + (iz - 1)] + 
                          f[ixs + iys + (iz + 1)]);


                temp +=
                    ec * (f[ixms + iys + iz - 1] +
                          f[ixps + iys + iz - 1] +
                          f[ixs + iyms + iz - 1] +
                          f[ixs + iyps + iz - 1] +
                          f[ixms + iyms + iz] +
                          f[ixms + iyps + iz] +
                          f[ixps + iyms + iz] +
                          f[ixps + iyps + iz] +
                          f[ixms + iys + iz + 1] +
                          f[ixps + iys + iz + 1] +
                          f[ixs + iyms + iz + 1] + 
                          f[ixs + iyps + iz + 1]);


                temp +=
                    crn * (f[ixms + iyms + iz - 1] +
                           f[ixms + iyms + iz + 1] +
                           f[ixms + iyps + iz - 1] +
                           f[ixms + iyps + iz + 1] +
                           f[ixps + iyms + iz - 1] +
                           f[ixps + iyms + iz + 1] +
                           f[ixps + iyps + iz - 1] +
                           f[ixps + iyps + iz + 1]);

                work[ixs + iys + iz] = scale * temp;


            }                   /* end for */

        }                       /* end for */

    }                           /* end for */

}


/******/
