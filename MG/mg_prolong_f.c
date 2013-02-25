/************************** SVN Revision Information **************************
 **    $Id: mg_prolong.c 1848 2013-01-27 14:46:54Z ebriggs $    **
******************************************************************************/

/****f* QMD-MGDFT/mg_prolong.c *****
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
 *   void mg_prolong(rmg_float_t *full, rmg_float_t *half, int dimx, int dimy, int dimz)
 *   prolongation operator: when we go from corse grid to fine grid,
 *                          we need to get full array from the half array
 * INPUTS
 *   half[(dimx/2+2) * (dimy/2+2) * (dimz/2+2)]: array in corse grid
 *     its image value must be there
 *   dimx, dimy, dimz: dimensions of array in the fine grid
 * OUTPUT
 *   full[(dimx+2) * (dimy+2) * (dimz+2)]: array in fine grid
 *     the missing points in the corse grid are interpolated from neighbors
 * PARENTS
 *   mgrid_solv.c
 * CHILDREN
 *   nothing
 * SOURCE
 */


#include "main.h"
#include <float.h>
#include <math.h>




void mg_prolong_f (rmg_float_t * full, rmg_float_t * half, int dimx, int dimy, int dimz, int dx2, int dy2, int dz2, int xoffset, int yoffset, int zoffset)
{

    int ix, iy, iz;
    int incx, incy, incxr, incyr;

    incy = dimz + 2;
    incx = (dimz + 2) * (dimy + 2);

    incyr = dz2 + 2;
    incxr = (dz2 + 2) * (dy2 + 2);


    /* transfer coarse grid points to fine grid along with the
     * high side image point
     */

    for (ix = 1; ix <= dx2 + 1; ix++)
    {

        for (iy = 1; iy <= dy2 + 1; iy++)
        {

            for (iz = 1; iz <= dz2 + 1; iz++)
            {

                full[(2 * ix - 1) * incx + (2 * iy - 1) * incy + 2 * iz - 1] =
                    half[ix * incxr + iy * incyr + iz];

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */



    /* interior center points
     */
    for (ix = 2; ix <= dimx; ix += 2)
    {

        for (iy = 2; iy <= dimy; iy += 2)
        {

            for (iz = 2; iz <= dimz; iz += 2)
            {

                full[ix * incx + iy * incy + iz] =
                    0.125 * full[(ix - 1) * incx + (iy - 1) * incy + iz - 1] +
                    0.125 * full[(ix - 1) * incx + (iy - 1) * incy + iz + 1] +
                    0.125 * full[(ix - 1) * incx + (iy + 1) * incy + iz - 1] +
                    0.125 * full[(ix - 1) * incx + (iy + 1) * incy + iz + 1] +
                    0.125 * full[(ix + 1) * incx + (iy - 1) * incy + iz - 1] +
                    0.125 * full[(ix + 1) * incx + (iy - 1) * incy + iz + 1] +
                    0.125 * full[(ix + 1) * incx + (iy + 1) * incy + iz - 1] +
                    0.125 * full[(ix + 1) * incx + (iy + 1) * incy + iz + 1];

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */



    for (ix = 1; ix <= dimx; ix += 2)
    {

        for (iy = 1; iy <= dimy; iy += 2)
        {

            for (iz = 2; iz <= dimz; iz += 2)
            {

                full[ix * incx + iy * incy + iz] =
                    0.5 * full[ix * incx + iy * incy + iz - 1] +
                    0.5 * full[ix * incx + iy * incy + iz + 1];

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */



    for (ix = 1; ix <= dimx; ix += 2)
    {

        for (iy = 2; iy <= dimy; iy += 2)
        {

            for (iz = 1; iz <= dimz; iz += 2)
            {

                full[ix * incx + iy * incy + iz] =
                    0.5 * full[ix * incx + (iy - 1) * incy + iz] +
                    0.5 * full[ix * incx + (iy + 1) * incy + iz];

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */



    for (ix = 2; ix <= dimx; ix += 2)
    {

        for (iy = 1; iy <= dimy; iy += 2)
        {

            for (iz = 1; iz <= dimz; iz += 2)
            {

                full[ix * incx + iy * incy + iz] =
                    0.5 * full[(ix - 1) * incx + iy * incy + iz] +
                    0.5 * full[(ix + 1) * incx + iy * incy + iz];

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */



    for (ix = 1; ix <= dimx; ix += 2)
    {

        for (iy = 2; iy <= dimy; iy += 2)
        {

            for (iz = 2; iz <= dimz; iz += 2)
            {

                full[ix * incx + iy * incy + iz] =
                    0.25 * full[ix * incx + (iy - 1) * incy + iz - 1] +
                    0.25 * full[ix * incx + (iy - 1) * incy + iz + 1] +
                    0.25 * full[ix * incx + (iy + 1) * incy + iz - 1] +
                    0.25 * full[ix * incx + (iy + 1) * incy + iz + 1];

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */



    for (ix = 2; ix <= dimx; ix += 2)
    {

        for (iy = 1; iy <= dimy; iy += 2)
        {

            for (iz = 2; iz <= dimz; iz += 2)
            {

                full[ix * incx + iy * incy + iz] =
                    0.25 * full[(ix - 1) * incx + iy * incy + iz - 1] +
                    0.25 * full[(ix - 1) * incx + iy * incy + iz + 1] +
                    0.25 * full[(ix + 1) * incx + iy * incy + iz - 1] +
                    0.25 * full[(ix + 1) * incx + iy * incy + iz + 1];

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */



    for (ix = 2; ix <= dimx; ix += 2)
    {

        for (iy = 2; iy <= dimy; iy += 2)
        {

            for (iz = 1; iz <= dimz; iz += 2)
            {

                full[ix * incx + iy * incy + iz] =
                    0.25 * full[(ix - 1) * incx + (iy - 1) * incy + iz] +
                    0.25 * full[(ix + 1) * incx + (iy - 1) * incy + iz] +
                    0.25 * full[(ix - 1) * incx + (iy + 1) * incy + iz] +
                    0.25 * full[(ix + 1) * incx + (iy + 1) * incy + iz];

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */



}                               /* end mg_prolong */


/******/
