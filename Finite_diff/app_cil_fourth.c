/************************** SVN Revision Information **************************
 **    $Id: app_cil_fourth.c 1722 2012-05-07 19:38:37Z ebriggs $    **
******************************************************************************/

#include "main.h"
#include <float.h>
#include <math.h>
#include <stdlib.h>


static REAL app_cil_fourth_global (REAL * a, REAL * b, REAL gridhx, REAL gridhy, REAL gridhz);


REAL app_cil_fourth (REAL * a, REAL * b, int dimx, int dimy, int dimz, REAL gridhx, REAL gridhy, REAL gridhz)
{

    int  numgrid;
    REAL *rptr;
    int iz, ix, iy, incx, incy, incxr, incyr;
    int ixs, iys, ixms, ixps, iyms, iyps;
    REAL ecxy, ecxz, ecyz, cc = 0.0, fcx, fcy, fcz;
    REAL ihx, ihy, ihz, a1, a2, a3;

    if((ct.ibrav != CUBIC_PRIMITIVE) && (ct.ibrav != ORTHORHOMBIC_PRIMITIVE)) {
        error_handler("Grid symmetry not programmed yet in app_cil_fourth.\n");
    }

    numgrid = dimx * dimy * dimz;
    if(numgrid == P0_BASIS && ct.anisotropy < 1.000001)
    {
        return app_cil_fourth_global (a, b, gridhx, gridhy, gridhz);
    }

    ihx = 1.0 / (gridhx * gridhx * ct.xside * ct.xside);
    ihy = 1.0 / (gridhy * gridhy * ct.yside * ct.yside);
    ihz = 1.0 / (gridhz * gridhz * ct.zside * ct.zside);


    incx = (dimz + 2) * (dimy + 2);
    incy = dimz + 2;
    incxr = dimz * dimy;
    incyr = dimz;

    my_malloc (rptr, (dimx + 2) * (dimy + 2) * (dimz + 2), REAL);

    trade_imagesx (a, rptr, dimx, dimy, dimz, 1, FULL_FD);



    if (ct.anisotropy < 1.000001)
    {

        ihx = 1.0 / (gridhx * gridhx * ct.xside * ct.xside);
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
                        cc * rptr[ixs + iys + iz] +
                        fcx * (rptr[ixms + iys + iz] +
                                rptr[ixps + iys + iz] +
                                rptr[ixs + iyms + iz] +
                                rptr[ixs + iyps + iz] +
                                rptr[ixs + iys + (iz - 1)] + rptr[ixs + iys + (iz + 1)]);

                    b[(ix - 1) * incxr + (iy - 1) * incyr + (iz - 1)] +=
                        ecxy * (rptr[ixms + iys + iz - 1] +
                                rptr[ixps + iys + iz - 1] +
                                rptr[ixs + iyms + iz - 1] +
                                rptr[ixs + iyps + iz - 1] +
                                rptr[ixms + iyms + iz] +
                                rptr[ixms + iyps + iz] +
                                rptr[ixps + iyms + iz] +
                                rptr[ixps + iyps + iz] +
                                rptr[ixms + iys + iz + 1] +
                                rptr[ixps + iys + iz + 1] +
                                rptr[ixs + iyms + iz + 1] + rptr[ixs + iyps + iz + 1]);


                }           /* end for */

            }               /* end for */

        }                   /* end for */
    }
    else
    {

        /* Compute coefficients for this grid spacing */
        ihx = 1.0 / (gridhx * gridhx * ct.xside * ct.xside);
        ihy = 1.0 / (gridhy * gridhy * ct.yside * ct.yside);
        ihz = 1.0 / (gridhz * gridhz * ct.zside * ct.zside);

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
                        cc * rptr[ixs + iys + iz] +
                        fcx * rptr[ixms + iys + iz] +
                        fcx * rptr[ixps + iys + iz] +
                        fcy * rptr[ixs + iyms + iz] +
                        fcy * rptr[ixs + iyps + iz] +
                        fcz * rptr[ixs + iys + (iz - 1)] + fcz * rptr[ixs + iys + (iz + 1)];

                    b[(ix - 1) * incxr + (iy - 1) * incyr + (iz - 1)] +=
                        ecxz * rptr[ixms + iys + iz - 1] +
                        ecxz * rptr[ixps + iys + iz - 1] +
                        ecyz * rptr[ixs + iyms + iz - 1] +
                        ecyz * rptr[ixs + iyps + iz - 1] +
                        ecxy * rptr[ixms + iyms + iz] +
                        ecxy * rptr[ixms + iyps + iz] +
                        ecxy * rptr[ixps + iyms + iz] +
                        ecxy * rptr[ixps + iyps + iz] +
                        ecxz * rptr[ixms + iys + iz + 1] +
                        ecxz * rptr[ixps + iys + iz + 1] +
                        ecyz * rptr[ixs + iyms + iz + 1] + ecyz * a[ixs + iyps + iz + 1];


                }           /* end for */

            }               /* end for */

        }                   /* end for */

    }                       /* end if */



    my_free (rptr);
    return cc;

}                               /* end app_cil */





REAL app_cil_fourth_global (REAL * a, REAL * b, REAL gridhx, REAL gridhy, REAL gridhz)
{

    int ix, iy, iz;
    int ixs, iys, ixms, ixps, iyms, iyps;
    int incy, incx;
    int incyr, incxr;
    REAL *rptr, rz, rzps, rzms, rzpps;
    REAL c000, c100, c110;
    REAL ihx;

    incx = (PZ0_GRID + 2) * (PY0_GRID + 2);
    incy = PZ0_GRID + 2;
    incxr = PZ0_GRID * PY0_GRID;
    incyr = PZ0_GRID;

    my_malloc (rptr, (PX0_GRID + 2) * (PY0_GRID + 2) * (PZ0_GRID + 2), REAL);

    trade_imagesx (a, rptr, PX0_GRID, PY0_GRID, PZ0_GRID, 1, FULL_FD);

    ihx = 1.0 / (gridhx * gridhx * ct.xside * ct.xside);
    c000 = (-4.0 / 3.0) * (ihx + ihx + ihx);
    c100 = (5.0 / 6.0) * ihx + (c000 / 8.0);
    c110 = (1.0 / 12.0) * (ihx + ihx);


    // Handle the general case first
    for (ix = 1; ix < PX0_GRID + 1; ix++)
    {
        ixs = ix * incx;
        ixms = (ix - 1) * incx;
        ixps = (ix + 1) * incx;

        for (iy = 1; iy < PY0_GRID + 1; iy++)
        {
            iys = iy * incy;
            iyms = (iy - 1) * incy;
            iyps = (iy + 1) * incy;

            for (iz = 1; iz < PZ0_GRID + 1; iz++)
            {

                b[(ix - 1) * incxr + (iy - 1) * incyr + (iz - 1)] =
                    c100 * (rptr[ixs + iys + (iz - 1)] +
                            rptr[ixs + iys + (iz + 1)] +
                            rptr[ixms + iys + iz] +
                            rptr[ixps + iys + iz] +
                            rptr[ixs + iyms + iz] +
                            rptr[ixs + iyps + iz]) + c000 * rptr[ixs + iys + iz];

                b[(ix - 1) * incxr + (iy - 1) * incyr + (iz - 1)] +=
                    c110 * (rptr[ixps + iyps + iz] +
                            rptr[ixps + iyms + iz] +
                            rptr[ixms + iyps + iz] +
                            rptr[ixms + iyms + iz] +
                            rptr[ixps + iys + (iz + 1)] +
                            rptr[ixps + iys + (iz - 1)] +
                            rptr[ixms + iys + (iz + 1)] +
                            rptr[ixms + iys + (iz - 1)] +
                            rptr[ixs + iyps + (iz + 1)] +
                            rptr[ixs + iyps + (iz - 1)] +
                            rptr[ixs + iyms + (iz + 1)] + rptr[ixs + iyms + (iz - 1)]);

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */


    my_free (rptr);
    return c000;

}
