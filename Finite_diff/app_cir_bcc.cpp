
/****f* QMD-MGDFT/app_cir_bcc.c *****
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
 *    void app_cir_bcc(double *a, double *b, int dimx, int dimy, int dimz)
 *    Applies the Mehrstellen RHS operator to a matrix for BCC lattice 
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


#include "BaseGrid.h"
#include "Lattice.h"
#include "FiniteDiff.h"
#include <cmath>
#include <complex>



template <typename RmgType>
void FiniteDiff::app_cir_bcc (RmgType * a, RmgType * b, int dimx, int dimy, int dimz)
{

    int ix, iy, iz;
    int ixs, iys;
    int incy, incx;
    int incyr, incxr;
    double Bc, Bf;

    incy = dimz + 2;
    incx = (dimz + 2) * (dimy + 2);
    incyr = dimz;
    incxr = dimz * dimy;


    Bc = 2.0 / 3.0;
    Bf = 1.0 / 24.0;
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
                    Bf * a[(ix - 1) * incx + (iy - 1) * incy + iz - 1] +
                    Bf * a[(ix - 1) * incx + iy * incy + iz] +
                    Bf * a[ix * incx + (iy - 1) * incy + iz] +
                    Bf * a[ix * incx + iy * incy + iz - 1] +
                    Bf * a[ix * incx + iy * incy + iz + 1] +
                    Bf * a[ix * incx + (iy + 1) * incy + iz] +
                    Bf * a[(ix + 1) * incx + iy * incy + iz] +
                    Bf * a[(ix + 1) * incx + (iy + 1) * incy + iz + 1];


            }                   /* end for */

        }                       /* end for */

    }                           /* end for */


}                               /* end app_cir_bcc */

/******/
