/****f* QMD-MGDFT/app_cir_fcc.c *****
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
 *    void app_cir_fcc(double *a, double *b, int dimx, int dimy, int dimz)
 *    Applies the Mehrstellen RHS operator to a matrix for FCC lattice 
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
void FiniteDiff::app_cir_fcc (RmgType * a, RmgType * b, int dimx, int dimy, int dimz)
{

    int ix, iy, iz;
    int ixs, iys, ixms, ixps, iyms, iyps;
    int incy, incx;
    int incyr, incxr;

    double Bc, Bf;


    incy = dimz + 2;
    incx = (dimz + 2) * (dimy + 2);
    incyr = dimz;
    incxr = dimz * dimy;

    Bc = 2.0 / 3.0;
    Bf = 1.0 / 36.0;
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

                b[(ix - 1) * incxr + (iy - 1) * incyr + (iz - 1)] =
                    Bc * a[ixs + iys + iz] +
                    Bf * (a[ixms + iys + iz] +
                          a[ixms + iys + iz + 1] +
                          a[ixms + iyps + iz] +
                          a[ixs + iyms + iz] +
                          a[ixs + iyms + iz + 1] +
                          a[ixs + iys + iz - 1] +
                          a[ixs + iys + iz + 1] +
                          a[ixs + iyps + iz - 1] +
                          a[ixs + iyps + iz] +
                          a[ixps + iyms + iz] + a[ixps + iys + iz - 1] + a[ixps + iys + iz]);


            }                   /* end for */

        }                       /* end for */

    }                           /* end for */


}                               /* end app_cir_fcc */

