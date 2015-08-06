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

double app_cil_orbital (double * a, double * b, int dimx, int dimy, int dimz,
              double gridhx, double gridhy, double gridhz)
{

    int iz, ix, iy, incx, incy, incxr, incyr;
    int ixs, iys, ixms, ixps, iyms, iyps;
    double ecxy, ecxz, ecyz, cc = 0.0, fcx, fcy, fcz;
    double ihx, ihy, ihz, a1, a2, a3;



    switch (get_ibrav_type())
    {

    case CUBIC_PRIMITIVE:
    case ORTHORHOMBIC_PRIMITIVE:


        if (get_anisotropy() < 1.000001)
        {

            ihx = 1.0 / (gridhx * gridhx * get_xside() * get_xside());
            cc = (-4.0 / 3.0) * (ihx + ihx + ihx);
            fcx = (5.0 / 6.0) * ihx + (cc / 8.0);
            ecxy = (1.0 / 12.0) * (ihx + ihx);
            incy = dimz + 2;
            incx = (dimz + 2) * (dimy + 2);
            incyr = dimz;
            incxr = dimz * dimy;

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
                            cc * a[ixs + iys + iz] +
                            fcx * (a[ixms + iys + iz] +
                                   a[ixps + iys + iz] +
                                   a[ixs + iyms + iz] +
                                   a[ixs + iyps + iz] +
                                   a[ixs + iys + (iz - 1)] + a[ixs + iys + (iz + 1)]);

                        b[(ix - 1) * incxr + (iy - 1) * incyr + (iz - 1)] +=
                            ecxy * (a[ixms + iys + iz - 1] +
                                    a[ixps + iys + iz - 1] +
                                    a[ixs + iyms + iz - 1] +
                                    a[ixs + iyps + iz - 1] +
                                    a[ixms + iyms + iz] +
                                    a[ixms + iyps + iz] +
                                    a[ixps + iyms + iz] +
                                    a[ixps + iyps + iz] +
                                    a[ixms + iys + iz + 1] +
                                    a[ixps + iys + iz + 1] +
                                    a[ixs + iyms + iz + 1] + a[ixs + iyps + iz + 1]);


                    }           /* end for */

                }               /* end for */

            }                   /* end for */
        }
        else
        {

            /* Compute coefficients for this grid spacing */
            ihx = 1.0 / (gridhx * gridhx * get_xside() * get_xside());
            ihy = 1.0 / (gridhy * gridhy * get_yside() * get_yside());
            ihz = 1.0 / (gridhz * gridhz * get_zside() * get_zside());

            cc = (-4.0 / 3.0) * (ihx + ihy + ihz);

            fcx = (5.0 / 6.0) * ihx + (cc / 8.0);
            fcy = (5.0 / 6.0) * ihy + (cc / 8.0);
            fcz = (5.0 / 6.0) * ihz + (cc / 8.0);

            ecxy = (1.0 / 12.0) * (ihx + ihy);
            ecxz = (1.0 / 12.0) * (ihx + ihz);
            ecyz = (1.0 / 12.0) * (ihy + ihz);


            incy = dimz + 2;
            incx = (dimz + 2) * (dimy + 2);
            incyr = dimz;
            incxr = dimz * dimy;



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
                            cc * a[ixs + iys + iz] +
                            fcx * a[ixms + iys + iz] +
                            fcx * a[ixps + iys + iz] +
                            fcy * a[ixs + iyms + iz] +
                            fcy * a[ixs + iyps + iz] +
                            fcz * a[ixs + iys + (iz - 1)] + fcz * a[ixs + iys + (iz + 1)];

                        b[(ix - 1) * incxr + (iy - 1) * incyr + (iz - 1)] +=
                            ecxz * a[ixms + iys + iz - 1] +
                            ecxz * a[ixps + iys + iz - 1] +
                            ecyz * a[ixs + iyms + iz - 1] +
                            ecyz * a[ixs + iyps + iz - 1] +
                            ecxy * a[ixms + iyms + iz] +
                            ecxy * a[ixms + iyps + iz] +
                            ecxy * a[ixps + iyms + iz] +
                            ecxy * a[ixps + iyps + iz] +
                            ecxz * a[ixms + iys + iz + 1] +
                            ecxz * a[ixps + iys + iz + 1] +
                            ecyz * a[ixs + iyms + iz + 1] + ecyz * a[ixs + iyps + iz + 1];


                    }           /* end for */

                }               /* end for */

            }                   /* end for */

        }                       /* end if */
        break;

    case CUBIC_BC:

        /* Compute coefficients for this grid spacing */
        ihx = 1.0 / (gridhx * gridhx * get_xside() * get_xside());
        ihx = 3.0 * ihx / 4.0;

        cc = (-22.0 / 3.0) * ihx;
        a1 = (2.0 / 3.0) * ihx;
        a2 = (1.0 / 3.0) * ihx;

        incy = dimz + 2;
        incx = (dimz + 2) * (dimy + 2);
        incyr = dimz;
        incxr = dimz * dimy;


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
                        cc * a[ix * incx + iys + iz] +
                        a2 * (a[ixms + iyms + iz] +
                              a[ixms + iys + iz - 1] +
                              a[ixs + iyms + iz - 1] +
                              a[ixs + iyps + iz + 1] +
                              a[ixps + iys + (iz + 1)] + a[ixps + iyps + iz]);

                    b[(ix - 1) * incxr + (iy - 1) * incyr + (iz - 1)] +=
                        a1 * (a[(ix - 1) * incx + (iy - 1) * incy + iz - 1] +
                              a[(ix - 1) * incx + iy * incy + iz] +
                              a[ix * incx + (iy - 1) * incy + iz] +
                              a[ix * incx + iy * incy + iz - 1] +
                              a[ix * incx + iy * incy + iz + 1] +
                              a[ix * incx + (iy + 1) * incy + iz] +
                              a[(ix + 1) * incx + iy * incy + iz] +
                              a[(ix + 1) * incx + (iy + 1) * incy + iz + 1]);

                }               /* end for */

            }                   /* end for */

        }                       /* end for */
        break;

    case CUBIC_FC:
        /* Compute coefficients for this grid spacing */
        ihx = 1.0 / (gridhx * gridhx * get_xside() * get_xside());
        ihx = ihx / 2.0;

        cc = (-34.0 / 3.0) * ihx;
        a1 = (8.0 / 9.0) * ihx;
        a2 = (1.0 / 9.0) * ihx;

        incy = dimz + 2;
        incx = (dimz + 2) * (dimy + 2);
        incyr = dimz;
        incxr = dimz * dimy;


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
                        cc * a[ix * incx + iys + iz] +
                        a1 * (a[ixms + iys + iz] +
                              a[ixms + iys + iz + 1] +
                              a[ixms + iyps + iz] +
                              a[ixs + iyms + iz] +
                              a[ixs + iyms + iz + 1] +
                              a[ixs + iys + iz - 1] +
                              a[ixs + iys + iz + 1] +
                              a[ixs + iyps + iz - 1] +
                              a[ixs + iyps + iz] +
                              a[ixps + iyms + iz] +
                              a[ixps + iy * incy + iz - 1] + a[ixps + iy * incy + iz]);


                    b[(ix - 1) * incxr + (iy - 1) * incyr + (iz - 1)] +=
                        a2 * (a[ixms + iyms + iz + 1] +
                              a[ixms + iyps + iz - 1] +
                              a[ixms + iyps + iz + 1] +
                              a[ixps + iyms + iz - 1] +
                              a[ixps + iyms + iz + 1] + a[ixps + iyps + iz - 1]);

                }               /* end for */

            }                   /* end for */

        }                       /* end for */
        break;

    case HEXAGONAL:
        /* Compute coefficients for this grid spacing */
        ihx = 1.0 / (gridhx * gridhx * get_xside() * get_xside());
        ihz = 1.0 / (gridhz * gridhz * get_zside() * get_zside());

        cc = ((-3.0 / 4.0) * ihz) - ((5.0 / 3.0) * ihx);
        a1 = ((3.0 / 8.0) * ihz) - ((1.0 / 6.0) * ihx);
        a2 = ((5.0 / 18.0) * ihx) - ((1.0 / 24.0) * ihz);
        a3 = ((1.0 / 48.0) * ihz) + ((1.0 / 36.0) * ihx);
        cc *= 2.0;
        a1 *= 2.0;
        a2 *= 2.0;
        a3 *= 2.0;


        incy = dimz + 2;
        incx = (dimz + 2) * (dimy + 2);
        incyr = dimz;
        incxr = dimz * dimy;


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
                        cc * a[ixs + iys + iz] +
                        a3 * (a[ixps + iys + iz - 1] +
                              a[ixps + iyms + iz - 1] +
                              a[ixs + iyms + iz - 1] +
                              a[ixms + iys + iz - 1] +
                              a[ixms + iyps + iz - 1] +
                              a[ixs + iyps + iz - 1] +
                              a[ixps + iys + iz + 1] +
                              a[ixps + iyms + iz + 1] +
                              a[ixs + iyms + iz + 1] +
                              a[ixms + iys + iz + 1] +
                              a[ixms + iyps + iz + 1] + a[ixs + iyps + iz + 1]);


                    b[(ix - 1) * incxr + (iy - 1) * incyr + (iz - 1)] +=
                        a2 * (a[ixps + iys + iz] +
                              a[ixps + iyms + iz] +
                              a[ixs + iyms + iz] +
                              a[ixms + iys + iz] + a[ixms + iyps + iz] + a[ixs + iyps + iz]);

                    b[(ix - 1) * incxr + (iy - 1) * incyr + (iz - 1)] +=
                        a1 * (a[ixs + iys + iz - 1] + a[ixs + iys + iz + 1]);

                }               /* end for */

            }                   /* end for */

        }                       /* end for */

        break;
    default:
        error_handler ("Lattice type not implemented");

    }                           /* end switch */



    return cc;

}                               /* end app_cil */

/******/
