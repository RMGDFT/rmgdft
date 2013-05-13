/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

/****f* QMD-MGDFT/app6_del2.c *****
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
 *   void app6_del2(S0_GRID *f, P0_GRID *work)
 *   Apply the sixth order Laplacian operator on a orthorhombic grid
 * INPUTS
 *   S0_GRID *f:  real array of (PX0_GRID+2) * (PY0_GRID+2) * ((PZ0_GRID+2)
 *   see main.h for defination of unions S0_GRID and P0_GRID
 * OUTPUT
 *   P0_GRID *work: real array of pct.PX0_GRID * pct.PY0_GRID * pct.PZ0_GRID
 * PARENTS
 *   get_ke.c xcgga.c
 * CHILDREN
 *   trade_image2.c
 * SOURCE
 */



#include <float.h>
#include <math.h>
#include "main.h"




void app6_del2 (rmg_double_t *rho, rmg_double_t * work, int dimx, int dimy, int dimz,
                rmg_double_t gridhx, rmg_double_t gridhy, rmg_double_t gridhz)

{

    int iz, ix, iy, incx, incy, incxr, incyr;
    int ixs, iys, ixms, ixps, iyms, iyps, ixmms, ixpps, iymms, iypps;

    rmg_double_t h2, t0, t1x, t2x;
    rmg_double_t t1y, t2y;
    rmg_double_t t1z, t2z;
    rmg_double_t *f;

    incx = (dimz + 4) * (dimy + 4);
    incy = dimz + 4;
    incxr = dimz * dimy;
    incyr = dimz;

    my_malloc(f, (dimx+4) * (dimy+4 ) * (dimz+4), rmg_double_t);

    trade_imagesx (rho, f, dimx, dimy, dimz, 2, CENTRAL_FD);

    h2 = gridhx * gridhx * ct.xside * ct.xside;
    t0 = -30.0 / (12.0 * h2);
    t1x = 16.0 / (12.0 * h2);
    t2x = -1.0 / (12.0 * h2);

    h2 = gridhy * gridhy * ct.yside * ct.yside;
    t0 -= 30.0 / (12.0 * h2);
    t1y = 16.0 / (12.0 * h2);
    t2y = -1.0 / (12.0 * h2);

    h2 = gridhz * gridhz * ct.zside * ct.zside;
    t0 -= 30.0 / (12.0 * h2);
    t1z = 16.0 / (12.0 * h2);
    t2z = -1.0 / (12.0 * h2);



    for (ix = 2; ix < dimx + 2; ix++)
    {

        ixs = ix * incx;
        ixms = (ix - 1) * incx;
        ixps = (ix + 1) * incx;
        ixmms = (ix - 2) * incx;
        ixpps = (ix + 2) * incx;

        for (iy = 2; iy < dimy + 2; iy++)
        {

            iys = iy * incy;
            iyms = (iy - 1) * incy;
            iyps = (iy + 1) * incy;
            iymms = (iy - 2) * incy;
            iypps = (iy + 2) * incy;

            for (iz = 2; iz < dimz + 2; iz++)
            {

                work[(ix - 2) * incxr + (iy - 2) * incyr + (iz - 2)] = t0 * f[ixs + iys + iz] +
                    t1x * f[ixms + iys + iz] +
                    t1x * f[ixps + iys + iz] +
                    t1y * f[ixs + iyms + iz] +
                    t1y * f[ixs + iyps + iz] +
                    t1z * f[ixs + iys + iz - 1] +
                    t1z * f[ixs + iys + iz + 1] +
                    t2x * f[ixmms + iys + iz] +
                    t2x * f[ixpps + iys + iz] +
                    t2y * f[ixs + iymms + iz] +
                    t2y * f[ixs + iypps + iz] +
                    t2z * f[ixs + iys + iz - 2] +
                    t2z * f[ixs + iys + iz + 2];

            }                   /* end for */
        }                       /* end for */
    }                           /* end for */

    my_free(f);

}                               /* end app6_del2 */


/******/
