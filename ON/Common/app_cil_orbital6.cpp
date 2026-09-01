/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

/****f* QMD-MGDFT/app_cil.c *****
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
 *    double app_cil(double *a, double *b, int dimx, int dimy, int dimz,
 *                 double gridhx, double gridhy, double gridhz)
 *    Applies left hand side (LHS) Mehrstellen operator to function  
 * INPUTS
 *    a[(dimx+2) * (dimy+2) * (dimz+2)]: function to be applied
 *    dimx, dimy, dimz: array dimentions in x, y, z direction
 *    gridhx, gridhy, gridhz:  grid spacings in crystal coordinates
 * OUTPUT
 *    b[dimx * dimy * dimz]: LHS(a)
 * PARENTS
 *    get_vh.c
 * CHILDREN
 *    trade_images.c
 * SOURCE
 */

#include "main.h"
#include "prototypes_on.h"
#include <float.h>
#include <math.h>
#include <stdlib.h>


double app_cil_orbital6 (double * a, double * b, int dimx, int dimy, int dimz,
              double gridhx, double gridhy, double gridhz)
{

    int iz, ix, iy, incx, incy, incxr, incyr;
    int ixs, iys, ixms, ixps, iyms, iyps, ixmms, ixpps, iymms, iypps;
    double ecxy, ecxz, ecyz, cc, fcx, fcy, fcz, cor;
    double fc2x, fc2y, fc2z, tcx, tcy, tcz;
    double ihx, ihy, ihz;
    double *rptr;

    incx = (dimz + 4) * (dimy + 4);
    incy = dimz + 4;
    incxr = dimz * dimy;
    incyr = dimz;


    rptr = a;


    ihx = 1.0 / (gridhx * gridhx * get_xside() * get_xside());
    ihy = 1.0 / (gridhy * gridhy * get_yside() * get_yside());
    ihz = 1.0 / (gridhz * gridhz * get_zside() * get_zside());

    cc = (-116.0 / 90.0) * (ihx + ihy + ihz);

    fcx = (31.0 / 232.0) * cc + (49.0 / 60.0) * ihx;
    fcy = (31.0 / 232.0) * cc + (49.0 / 60.0) * ihy;
    fcz = (31.0 / 232.0) * cc + (49.0 / 60.0) * ihz;

    ecxy = (-31.0 / 464.0) * cc - (1.0 / 10.0) * ihz;
    ecxz = (-31.0 / 464.0) * cc - (1.0 / 10.0) * ihy;
    ecyz = (-31.0 / 464.0) * cc - (1.0 / 10.0) * ihx;

    cor = (1.0 / 144.0) * (ihx + ihy + ihz);

    fc2x = (1.0 / 120.0) * (ihy + ihz);
    fc2y = (1.0 / 120.0) * (ihx + ihz);
    fc2z = (1.0 / 120.0) * (ihx + ihy);

    tcx = (-1.0 / 240.0) * ihx;
    tcy = (-1.0 / 240.0) * ihy;
    tcz = (-1.0 / 240.0) * ihz;

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
                b[(ix - 2) * incxr + (iy - 2) * incyr + (iz - 2)] = cc * rptr[ixs + iys + iz] +
                    fcx * (rptr[ixms + iys + iz] + rptr[ixps + iys + iz]) +
                    fcy * (rptr[ixs + iyms + iz] + rptr[ixs + iyps + iz]) +
                    fcz * (rptr[ixs + iys + (iz - 1)] + rptr[ixs + iys + (iz + 1)]);

                b[(ix - 2) * incxr + (iy - 2) * incyr + (iz - 2)] +=
                    ecxz * (rptr[ixms + iys + (iz - 1)] + rptr[ixps + iys + (iz - 1)] +
                            rptr[ixms + iys + (iz + 1)] + rptr[ixps + iys + (iz + 1)]) +
                    ecyz * (rptr[ixs + iyms + (iz - 1)] + rptr[ixs + iyps + (iz - 1)] +
                            rptr[ixs + iyms + (iz + 1)] + rptr[ixs + iyps + (iz + 1)]) +
                    ecxy * (rptr[ixms + iyms + iz] + rptr[ixms + iyps + iz] +
                            rptr[ixps + iyms + iz] + rptr[ixps + iyps + iz]);

                b[(ix - 2) * incxr + (iy - 2) * incyr + (iz - 2)] +=
                    cor * (rptr[ixms + iyms + (iz - 1)] + rptr[ixps + iyms + (iz - 1)] +
                           rptr[ixms + iyps + (iz - 1)] + rptr[ixps + iyps + (iz - 1)] +
                           rptr[ixms + iyms + (iz + 1)] + rptr[ixps + iyms + (iz + 1)] +
                           rptr[ixms + iyps + (iz + 1)] + rptr[ixps + iyps + (iz + 1)]);

                b[(ix - 2) * incxr + (iy - 2) * incyr + (iz - 2)] +=
                    fc2x * (rptr[ixmms + iys + iz] + rptr[ixpps + iys + iz]) +
                    fc2y * (rptr[ixs + iymms + iz] + rptr[ixs + iypps + iz]) +
                    fc2z * (rptr[ixs + iys + (iz - 2)] + rptr[ixs + iys + (iz + 2)]);

                b[(ix - 2) * incxr + (iy - 2) * incyr + (iz - 2)] +=
                    tcx * (rptr[ixps + iypps + iz] + rptr[ixps + iymms + iz] +
                           rptr[ixms + iypps + iz] + rptr[ixms + iymms + iz] +
                           rptr[ixps + iys + (iz + 2)] + rptr[ixps + iys + (iz - 2)] +
                           rptr[ixms + iys + (iz + 2)] + rptr[ixms + iys + (iz - 2)]) +
                    tcy * (rptr[ixpps + iyps + iz] + rptr[ixmms + iyps + iz] +
                           rptr[ixpps + iyms + iz] + rptr[ixmms + iyms + iz] +
                           rptr[ixs + iyps + (iz + 2)] + rptr[ixs + iyps + (iz - 2)] +
                           rptr[ixs + iyms + (iz + 2)] + rptr[ixs + iyms + (iz - 2)]) +
                    tcz * (rptr[ixpps + iys + (iz + 1)] + rptr[ixmms + iys + (iz + 1)] +
                           rptr[ixpps + iys + (iz - 1)] + rptr[ixmms + iys + (iz - 1)] +
                           rptr[ixs + iypps + (iz + 1)] + rptr[ixs + iymms + (iz + 1)] +
                           rptr[ixs + iypps + (iz - 1)] + rptr[ixs + iymms + (iz - 1)]);

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */


    return cc;

}

