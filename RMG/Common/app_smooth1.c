/************************** SVN Revision Information **************************
 **    $Id: app_smooth.c 1066 2009-08-31 18:41:09Z froze $    **
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
void app_smooth1 (REAL * f, REAL * work, int dimx, int dimy, int dimz)
{

    int iz, ix, iy, incx, incy;
    int ixs, iys, ixms, ixps, iyms, iyps;

    REAL scale, ec, fc, crn, cc, temp;

    incy = dimz + 2;
    incx = (dimz + 2) * (dimy + 2);

    cc = 6.0;
    fc = 1.0;
    ec = 0.0;
    crn = 0.0;
    scale = 1.0 / 12.0;


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

                work[ixs + iys + iz] = scale * temp;

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */



}                               /* end app_smooth1 */






/******/
