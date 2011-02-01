/************************** SVN Revision Information **************************
 **    $Id: app_cir_hex.c 587 2006-08-18 22:57:56Z miro $    **
******************************************************************************/
 
/****f* QMD-MGDFT/app_cir_hex.c *****
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
 *    void app_cir_hex(REAL *a, REAL *b, int dimx, int dimy, int dimz)
 *    Applies the Mehrstellen RHS operator to a matrix for HEX lattice 
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



/** 
 *
 *    Applies the Mehrstellen RHS operator to matrix a dimensioned
 *      a[(dimx+2)][(dimy+2)][(dimz+2)]
 *    
 *    Returns result in matrix b dimensioned
 *      b[dimx][dimy][dimz] 
 *      
 *    For hexagonal grids.
 *
 * @return nothing
 */


void app_cir_hex(REAL * a, REAL * b, int dimx, int dimy, int dimz)
{

    int ix, iy, iz;
    int ixs, iys;
    int incy, incx;
    int incyr, incxr;
    REAL Bc, Bf, Bz;

    incy = dimz + 2;
    incx = (dimz + 2) * (dimy + 2);
    incyr = dimz;
    incxr = dimz * dimy;


    trade_images(a, dimx, dimy, dimz, pct.neighbors);


    Bc = 7.0 / 12.0;
    Bf = 1.0 / 24.0;
    Bz = 1.0 / 12.0;
    for (ix = 1; ix <= dimx; ix++)
    {
        ixs = ix * incx;
        for (iy = 1; iy <= dimy; iy++)
        {
            iys = iy * incy;
            for (iz = 1; iz <= dimz; iz++)
            {

                b[(ix - 1) * incxr + (iy - 1) * incyr + (iz - 1)] =
                    Bc * a[ixs + iys + iz] +
                    Bz * a[ixs + iys + iz - 1] +
                    Bz * a[ixs + iys + iz + 1] +
                    Bf * a[(ix + 1) * incx + iys + iz] +
                    Bf * a[(ix + 1) * incx + (iy - 1) * incy + iz] +
                    Bf * a[ixs + (iy - 1) * incy + iz] +
                    Bf * a[(ix - 1) * incx + iys + iz] +
                    Bf * a[(ix - 1) * incx + (iy + 1) * incy + iz] +
                    Bf * a[ixs + (iy + 1) * incy + iz];


            }                   /* end for */

        }                       /* end for */

    }                           /* end for */

}                               /* end app_cir_hex */

/******/
