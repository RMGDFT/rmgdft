/************************** SVN Revision Information **************************
 **    $Id: app_cir_ortho_0.c 1242 2011-02-02 18:55:23Z luw $    **
******************************************************************************/
 
/****f* QMD-MGDFT/app_cir_ortho.c *****
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
 *    void app_cir_ortho(REAL *a, REAL *b, int dimx, int dimy, int dimz)
 *    Applies the Mehrstellen RHS operator to a matrix for orthorhombic lattice 
 * INPUTS
 *    a[(dimx+2) * (dimy+2) * (dimz+2)]: the matrix to be applied
 *    dimx, dimy, dimz: array dimentions in x, y, z direction
 * OUTPUT
 *    b[dimx * dimy * dimz]:  RHS(a)
 * PARENTS
 *    get_cir.c
 * CHILDREN
 *    trade_images.c
 * SOURCE
 */


#include "md.h"
#include <float.h>
#include <math.h>
#include <stdlib.h>

void app_cir_ortho_0 (REAL * a, REAL * b, int dimx, int dimy, int dimz)
{

    int ix, iy, iz;
    int ixs, iys, ixms, ixps, iyms, iyps, iyps1, iyps2;
    int incy, incx;
    int incyr, incxr;
    REAL Bc, Bf;

    incy = dimz + 2;
    incx = (dimz + 2) * (dimy + 2);
    incyr = dimz;
    incxr = dimz * dimy;



    Bc = 1.0;
    Bf = 0.0 / 12.0;
    for (ix = 1; ix <= dimx; ix++)
    {
        ixs = ix * incx;
        ixms = (ix - 1) * incx;
        ixps = (ix + 1) * incx;

        /* We optimize cache utilization here by unrolling the outer loops */
        if (dimy % 4)
        {

            for (iy = 1; iy <= dimy; iy++)
            {
                iys = iy * incy;
                iyms = (iy - 1) * incy;
                iyps = (iy + 1) * incy;
                for (iz = 1; iz <= dimz; iz++)
                {

                    b[(ix - 1) * incxr + (iy - 1) * incyr + (iz - 1)] =
                        Bf * (a[ixs + iys + (iz - 1)] +
                              a[ixs + iys + (iz + 1)] +
                              a[ixms + iys + iz] +
                              a[ixps + iys + iz] +
                              a[ixs + iyms + iz] + a[ixs + iyps + iz]) + Bc * a[ixs + iys + iz];

                }               /* end for */

            }                   /* end for */

        }
        else
        {

            for (iy = 1; iy <= dimy; iy += 4)
            {

                iys = iy * incy;
                iyms = (iy - 1) * incy;
                iyps = (iy + 1) * incy;
                iyps1 = (iy + 2) * incy;
                iyps2 = (iy + 3) * incy;
                for (iz = 1; iz <= dimz; iz++)
                {

                    b[(ix - 1) * incxr + (iy - 1) * incyr + (iz - 1)] =
                        Bf * (a[ixs + iys + (iz - 1)] +
                              a[ixs + iys + (iz + 1)] +
                              a[ixms + iys + iz] +
                              a[ixps + iys + iz] +
                              a[ixs + iyms + iz] + a[ixs + iyps + iz]) + Bc * a[ixs + iys + iz];

                    b[(ix - 1) * incxr + iy * incyr + (iz - 1)] =
                        Bf * (a[ixs + iyps + (iz - 1)] +
                              a[ixs + iyps + (iz + 1)] +
                              a[ixms + iyps + iz] +
                              a[ixps + iyps + iz] +
                              a[ixs + iys + iz] +
                              a[ixs + iyps + incy + iz]) + Bc * a[ixs + iyps + iz];

                    b[(ix - 1) * incxr + (iy + 1) * incyr + (iz - 1)] =
                        Bf * (a[ixs + iyps1 + (iz - 1)] +
                              a[ixs + iyps1 + (iz + 1)] +
                              a[ixms + iyps1 + iz] +
                              a[ixps + iyps1 + iz] +
                              a[ixs + iyps + iz] +
                              a[ixs + iyps1 + incy + iz]) + Bc * a[ixs + iyps1 + iz];

                    b[(ix - 1) * incxr + (iy + 2) * incyr + (iz - 1)] =
                        Bf * (a[ixs + iyps2 + (iz - 1)] +
                              a[ixs + iyps2 + (iz + 1)] +
                              a[ixms + iyps2 + iz] +
                              a[ixps + iyps2 + iz] +
                              a[ixs + iyps1 + iz] +
                              a[ixs + iyps2 + incy + iz]) + Bc * a[ixs + iyps2 + iz];

                }               /* end for */

            }                   /* end for */

        }                       /* end if */

    }                           /* end for */

}                               /* end app_cir_ortho */



/******/
