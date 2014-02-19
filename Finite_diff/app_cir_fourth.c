/************************** SVN Revision Information **************************
 **    $Id: app_cir_sixth.c 1722 2012-05-07 19:38:37Z ebriggs $    **
******************************************************************************/

#include "const.h"
#include "common_prototypes.h"
#include "fixed_dims.h"
#include "rmg_alloc.h"
#include <float.h>
#include <math.h>
#include <stdlib.h>
#if HYBRID_MODEL
#include "hybrid.h"
#endif

static void app_cir_fourth_global (rmg_double_t * a, rmg_double_t * b);

void app_cir_fourth (rmg_double_t * a, rmg_double_t * b, int dimx, int dimy, int dimz)
{

    int ix, iy, iz, numgrid, tid, ibrav, P0_BASIS;
    int ixs, iys, ixms, ixps, iyms, iyps;
    int incy, incx;
    int incyr, incxr;
    rmg_double_t *rptr;
    rmg_double_t c000, c100;

    ibrav = get_ibrav_type();
    P0_BASIS = get_P0_BASIS();

    if((ibrav != CUBIC_PRIMITIVE) && (ibrav != ORTHORHOMBIC_PRIMITIVE)) {
        error_handler("Grid symmetry not programmed yet in app_cir_fourth.\n");
    }

#if HYBRID_MODEL
    tid = get_thread_tid();
    if(tid < 0) tid = 0;  // OK in this case
#else
    tid = 0;
#endif

    numgrid = dimx * dimy * dimz;
    if(numgrid == P0_BASIS) {
        app_cir_fourth_global (a, b);
        return;
    }

    incx = (dimz + 2) * (dimy + 2);
    incy = dimz + 2;
    incxr = dimz * dimy;
    incyr = dimz;

    my_malloc (rptr, (dimx + 2) * (dimy + 2) * (dimz + 2), rmg_double_t);

    trade_imagesx (a, rptr, dimx, dimy, dimz, 1, CENTRAL_FD);

    c000 = 0.5;
    c100 = 1.0 / 12.0;
    for (ix = 1; ix < dimx + 1; ix++)
    {
        ixs = ix * incx;
        ixms = (ix - 1) * incx;
        ixps = (ix + 1) * incx;

        for (iy = 1; iy < dimy + 1; iy++)
        {
            iys = iy * incy;
            iyms = (iy - 1) * incy;
            iyps = (iy + 1) * incy;

            for (iz = 1; iz < dimz + 1; iz++)
            {

                b[(ix - 1) * incxr + (iy - 1) * incyr + (iz - 1)] =
                    c100 * (rptr[ixs + iys + (iz - 1)] +
                            rptr[ixs + iys + (iz + 1)] +
                            rptr[ixms + iys + iz] +
                            rptr[ixps + iys + iz] +
                            rptr[ixs + iyms + iz] +
                            rptr[ixs + iyps + iz]) + 
                    c000 *  rptr[ixs + iys + iz];

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */

    my_free (rptr);
}

void app_cir_fourth_global (rmg_double_t * a, rmg_double_t * b)
{

    int ix, iy, iz;
    int ixs, iys, ixms, ixps, iyms, iyps;
    int incy, incx;
    int incyr, incxr;
    rmg_double_t *rptr, rz, rzps, rzms, rzpps;
    rmg_double_t c000, c100;

    incx = (FIXED_ZDIM + 2) * (FIXED_YDIM + 2);
    incy = FIXED_ZDIM + 2;
    incxr = FIXED_ZDIM * FIXED_YDIM;
    incyr = FIXED_ZDIM;

    my_malloc (rptr, (FIXED_XDIM + 2) * (FIXED_YDIM + 2) * (FIXED_ZDIM + 2), rmg_double_t);

    trade_imagesx (a, rptr, FIXED_XDIM, FIXED_YDIM, FIXED_ZDIM, 1, CENTRAL_FD);

    c000 = 0.5;
    c100 = 1.0 / 12.0;


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
                            rptr[ixs + iyps + iz]) + 
                    c000 * rptr[ixs + iys + iz];

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */


    my_free (rptr);
    return;

}
