/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

/****f* QMD-MGDFT/pack_stop_axpy.c *****
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
 *   void pack_stop_axpy(rmg_double_t *sg, rmg_double_t *pg, rmg_double_t alpha, int dimx,
 *                       int dimy, int dimz)
 *   Adds the contents of the smoothing grid sg * alpha
 *           into the physical grid pg
 * INPUTS
 *   sg[(dimx+2)*(dimy+2)*(dimz+2)]: array in the smoothing grid
 *   pg[dimx*dimy*dimz]: array in the physical grid
 *   alpha: scaling factor
 *   dimx, dimy, dimz: dimensions of x,y,z direction
 * OUTPUT
 *   pg = pg + sg * alpha
 * PARENTS
 *   mg_eig_state.c
 * CHILDREN
 *   nothing
 * SOURCE
 */

#include "main.h"
#include <float.h>
#include <math.h>

void pack_stop_axpy (rmg_double_t * sg, rmg_double_t * pg, rmg_double_t alpha, int dimx, int dimy, int dimz)
{

    int ix, iy, iz, ixh, iyh;
    int incx, incy, incxs, incys;

    incy = dimz;
    incx = dimy * dimz;
    incys = dimz + 2;
    incxs = (dimy + 2) * (dimz + 2);


    /* Transfer pg into smoothing grid */
    for (ix = 0; ix < dimx; ix++)
    {

        ixh = ix + 1;
        for (iy = 0; iy < dimy; iy++)
        {

            iyh = iy + 1;
            for (iz = 0; iz < dimz; iz++)
            {

                pg[ix * incx + iy * incy + iz] += alpha * sg[ixh * incxs + iyh * incys + iz + 1];

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */


}                               /* end pack_stop_axpy */



void pack_stop_axpy_f (rmg_float_t * sg, rmg_float_t * pg, rmg_double_t alpha, int dimx, int dimy, int dimz)
{

    int ix, iy, iz, ixh, iyh;
    int incx, incy, incxs, incys;

    incy = dimz;
    incx = dimy * dimz;
    incys = dimz + 2;
    incxs = (dimy + 2) * (dimz + 2);


    /* Transfer pg into smoothing grid */
    for (ix = 0; ix < dimx; ix++)
    {

        ixh = ix + 1;
        for (iy = 0; iy < dimy; iy++)
        {

            iyh = iy + 1;
            for (iz = 0; iz < dimz; iz++)
            {

                pg[ix * incx + iy * incy + iz] += alpha * sg[ixh * incxs + iyh * incys + iz + 1];

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */


}                               /* end pack_stop_axpy */




/******/
