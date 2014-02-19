/************************** SVN Revision Information **************************
 **    $Id: app_cil_fourth.c 1722 2012-05-07 19:38:37Z ebriggs $    **
******************************************************************************/

#include "const.h"
#include "rmgtypes.h"
#include "rmgtypedefs.h"
#include "common_prototypes.h"
#include "rmg_alloc.h"
#include "rmg_error.h"
#include "fixed_dims.h"
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include "common_prototypes.h"
#include "hybrid.h"

static rmg_double_t app_cil_fourth_global_f (rmg_float_t * a, rmg_float_t * b, rmg_double_t gridhx, rmg_double_t gridhy, rmg_double_t gridhz);


rmg_double_t app_cil_fourth_f (rmg_float_t * a, rmg_float_t * b, int dimx, int dimy, int dimz, rmg_double_t gridhx, rmg_double_t gridhy, rmg_double_t gridhz)
{

    int  numgrid, tid, P0_BASIS, ibrav;
    rmg_float_t *rptr;
    int iz, ix, iy, incx, incy, incxr, incyr;
    int ixs, iys, ixms, ixps, iyms, iyps;
    rmg_double_t ecxy, ecxz, ecyz, cc = 0.0, fcx, fcy, fcz;
    rmg_double_t ihx, ihy, ihz;

    ibrav = get_ibrav_type();

    if((ibrav != CUBIC_PRIMITIVE) && (ibrav != ORTHORHOMBIC_PRIMITIVE)) {
        error_handler("Grid symmetry not programmed yet in app_cil_fourth_f.\n");
    }

    P0_BASIS = get_P0_BASIS();

#if HYBRID_MODEL
    tid = get_thread_tid();
    if(tid < 0) tid = 0;  // OK in this case
#else
    tid = 0;
#endif

    numgrid = dimx * dimy * dimz;
    if(numgrid == P0_BASIS && (get_anisotropy() < 1.000001))
    {
        return app_cil_fourth_global_f (a, b, gridhx, gridhy, gridhz);
    }

    ihx = 1.0 / (gridhx * gridhx * get_xside() * get_xside());
    ihy = 1.0 / (gridhy * gridhy * get_yside() * get_yside());
    ihz = 1.0 / (gridhz * gridhz * get_zside() * get_zside());


    incx = (dimz + 2) * (dimy + 2);
    incy = dimz + 2;
    incxr = dimz * dimy;
    incyr = dimz;

    my_malloc (rptr, (dimx + 2) * (dimy + 2) * (dimz + 2), rmg_float_t);

    trade_imagesx_f (a, rptr, dimx, dimy, dimz, 1, FULL_FD);



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





rmg_double_t app_cil_fourth_global_f (rmg_float_t * a, rmg_float_t * b, rmg_double_t gridhx, rmg_double_t gridhy, rmg_double_t gridhz)
{

    int ix, iy, iz;
    int ixs, iys, ixms, ixps, iyms, iyps;
    int incy, incx;
    int incyr, incxr;
    rmg_float_t *rptr;
    rmg_double_t c000, c100, c110;
    rmg_double_t ihx;

    incx = (FIXED_ZDIM + 2) * (FIXED_YDIM + 2);
    incy = FIXED_ZDIM + 2;
    incxr = FIXED_ZDIM * FIXED_YDIM;
    incyr = FIXED_ZDIM;

    my_malloc (rptr, (FIXED_XDIM + 2) * (FIXED_YDIM + 2) * (FIXED_ZDIM + 2), rmg_float_t);

    trade_imagesx_f (a, rptr, FIXED_XDIM, FIXED_YDIM, FIXED_ZDIM, 1, FULL_FD);

    ihx = 1.0 / (gridhx * gridhx * get_xside() * get_xside());
    c000 = (-4.0 / 3.0) * (ihx + ihx + ihx);
    c100 = (5.0 / 6.0) * ihx + (c000 / 8.0);
    c110 = (1.0 / 12.0) * (ihx + ihx);


    // Handle the general case first
    for (ix = 1; ix < FIXED_XDIM + 1; ix++)
    {
        ixs = ix * incx;
        ixms = (ix - 1) * incx;
        ixps = (ix + 1) * incx;

        for (iy = 1; iy < FIXED_YDIM + 1; iy++)
        {
            iys = iy * incy;
            iyms = (iy - 1) * incy;
            iyps = (iy + 1) * incy;

            for (iz = 1; iz < FIXED_ZDIM + 1; iz++)
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
