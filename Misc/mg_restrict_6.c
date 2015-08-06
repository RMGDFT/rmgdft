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
 *   void mg_restrict(double *full, double *half, int dimx, int dimy, int dimz)
 *   get half dimensioned array in the corse grid from full dimensioned 
 *   array in the fine grid. The returned values are smoothed.
 * INPUTS
 *   full[(dimx+10)*(dimy+10)*(dimz+10)]: array in the fine grid
 *   dimx, dimy, dimz: dimensions of array in the fine grid
 * OUTPUT
 *   half[(dimx/grid_ratio)*(dimy/grid_ratio)*(dimz/grid_ratio)] array in the corse grid
 * PARENTS
 *   mgrid_solv.c
 * CHILDREN
 *   nothing 
 * SOURCE
 */



#include "main.h"
#include "common_prototypes.h"
#include <float.h>
#include <math.h>
#include <TradeImages.h>





void mg_restrict_6 (double * full, double * half, int dimx, int dimy, int dimz, int grid_ratio)
{

    int ix, iy, iz;
    int incy, incx, incy2, incx2, incx3;
    double a0, a1, a2, a3, a4, a5;
    double *fulla;
    double *fullb;
    double *sg_full;

    my_malloc (sg_full, (dimx + 10) * (dimy + 10) * (dimz + 10), double);
    trade_imagesx (full, sg_full, get_FPX0_GRID(), get_FPY0_GRID(), get_FPZ0_GRID(), 5, FULL_TRADE);

    incy = dimz + 10;
    incx = (dimz + 10) * (dimy + 10);

    incy2 = dimz / grid_ratio;
    incx2 = (dimz / grid_ratio) * (dimy / grid_ratio);

    incx3 = (dimz + 10) * (dimy / grid_ratio);

    a0 = 3.0 / 7.0;
    a1 = (1.0 - a0) * 15.0 / 20.0;
    a2 = (1.0 - a0) * (-6.0) / 20.0;
    a3 = (1.0 - a0) * 1.0 / 20.0;
    a4 = 0.0;
    a5 = 0.0;

    my_malloc (fulla, (get_FPX0_GRID() / grid_ratio) * (get_FPY0_GRID() + 10) * (get_FPZ0_GRID() + 10), double);
    my_malloc (fullb, (get_FPX0_GRID() / grid_ratio) * (get_FPY0_GRID() / grid_ratio) * (get_FPZ0_GRID() + 10), double);


    for (ix = 0; ix < dimx / grid_ratio; ix++)
    {

        for (iy = 0; iy < dimy + 10; iy++)
        {


            for (iz = 0; iz < dimz + 10; iz++)
            {


                fulla[ix * incx + iy * incy + iz] =
                    a5 * sg_full[(grid_ratio * ix + 0) * incx + iy * incy + iz] +
                    a4 * sg_full[(grid_ratio * ix + 1) * incx + iy * incy + iz] +
                    a3 * sg_full[(grid_ratio * ix + 2) * incx + iy * incy + iz] +
                    a2 * sg_full[(grid_ratio * ix + 3) * incx + iy * incy + iz] +
                    a1 * sg_full[(grid_ratio * ix + 4) * incx + iy * incy + iz] +
                    a0 * sg_full[(grid_ratio * ix + 5) * incx + iy * incy + iz] +
                    a1 * sg_full[(grid_ratio * ix + 6) * incx + iy * incy + iz] +
                    a2 * sg_full[(grid_ratio * ix + 7) * incx + iy * incy + iz] +
                    a3 * sg_full[(grid_ratio * ix + 8) * incx + iy * incy + iz] +
                    a4 * sg_full[(grid_ratio * ix + 9) * incx + iy * incy + iz] +
                    a5 * sg_full[(grid_ratio * ix + 10) * incx + iy * incy + iz];


            }                   /* end for */

        }                       /* end for */

    }                           /* end for */



    for (ix = 0; ix < dimx / grid_ratio; ix++)
    {

        for (iy = 0; iy < dimy / grid_ratio; iy++)
        {


            for (iz = 0; iz < dimz + 10; iz++)
            {


                fullb[ix * incx3 + iy * incy + iz] =
                    a5 * fulla[ix * incx + (grid_ratio * iy + 0) * incy + iz] +
                    a4 * fulla[ix * incx + (grid_ratio * iy + 1) * incy + iz] +
                    a3 * fulla[ix * incx + (grid_ratio * iy + 2) * incy + iz] +
                    a2 * fulla[ix * incx + (grid_ratio * iy + 3) * incy + iz] +
                    a1 * fulla[ix * incx + (grid_ratio * iy + 4) * incy + iz] +
                    a0 * fulla[ix * incx + (grid_ratio * iy + 5) * incy + iz] +
                    a1 * fulla[ix * incx + (grid_ratio * iy + 6) * incy + iz] +
                    a2 * fulla[ix * incx + (grid_ratio * iy + 7) * incy + iz] +
                    a3 * fulla[ix * incx + (grid_ratio * iy + 8) * incy + iz] +
                    a4 * fulla[ix * incx + (grid_ratio * iy + 9) * incy + iz] +
                    a5 * fulla[ix * incx + (grid_ratio * iy + 10) * incy + iz];


            }                   /* end for */

        }                       /* end for */

    }                           /* end for */



    for (ix = 0; ix < dimx / grid_ratio; ix++)
    {

        for (iy = 0; iy < dimy / grid_ratio; iy++)
        {


            for (iz = 0; iz < dimz / grid_ratio; iz++)
            {


                half[ix * incx2 + iy * incy2 + iz] =
                    a5 * fullb[ix * incx3 + iy * incy + grid_ratio * iz + 0] +
                    a4 * fullb[ix * incx3 + iy * incy + grid_ratio * iz + 1] +
                    a3 * fullb[ix * incx3 + iy * incy + grid_ratio * iz + 2] +
                    a2 * fullb[ix * incx3 + iy * incy + grid_ratio * iz + 3] +
                    a1 * fullb[ix * incx3 + iy * incy + grid_ratio * iz + 4] +
                    a0 * fullb[ix * incx3 + iy * incy + grid_ratio * iz + 5] +
                    a1 * fullb[ix * incx3 + iy * incy + grid_ratio * iz + 6] +
                    a2 * fullb[ix * incx3 + iy * incy + grid_ratio * iz + 7] +
                    a3 * fullb[ix * incx3 + iy * incy + grid_ratio * iz + 8] +
                    a4 * fullb[ix * incx3 + iy * incy + grid_ratio * iz + 9] +
                    a5 * fullb[ix * incx3 + iy * incy + grid_ratio * iz + 10];


            }                   /* end for */

        }                       /* end for */

    }                           /* end for */


    my_free (fulla);
    my_free (fullb);
    my_free (sg_full);




}                               /* end mg_restrict */


/******/
