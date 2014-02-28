/************************** SVN Revision Information **************************
 **    $Id: app_cir_sixth.c 1722 2012-05-07 19:38:37Z ebriggs $    **
******************************************************************************/

#include "const.h"
#include "common_prototypes.h"
#include "common_prototypes.h"
#include "rmg_alloc.h"
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include "hybrid.h"


void app_cir_beta_fourth (rmg_double_t * a, rmg_double_t * b, int dimx, int dimy, int dimz)
{

    int ix, iy, iz, numgrid, tid, ibrav;
    int ixs, iys, ixms, ixps, iyms, iyps;
    int incy, incx;
    int incyr, incxr;
    rmg_double_t *rptr;
    rmg_double_t c000, c100;

    ibrav = get_ibrav_type();
    if((ibrav != CUBIC_PRIMITIVE) && (ibrav != ORTHORHOMBIC_PRIMITIVE)) {
        error_handler("Grid symmetry not programmed yet in app_cir_fourth.\n");
    }

    incx = (dimz + 2) * (dimy + 2);
    incy = dimz + 2;
    incxr = dimz * dimy;
    incyr = dimz;

    my_malloc (rptr, (dimx + 2) * (dimy + 2) * (dimz + 2), rmg_double_t);

    for(ix = 0; ix < (dimx + 2) * (dimy + 2) * (dimz + 2); ix++)
        rptr[ix] = 0.0;
    for (ix = 1; ix < dimx + 1; ix++)
    {
        ixs = ix * incx;

        for (iy = 1; iy < dimy + 1; iy++)
        {
            iys = iy * incy;

            for (iz = 1; iz < dimz + 1; iz++)
            {

                rptr[ixs + iys + iz] = a[(ix-1)* incxr + (iy-1) * incyr + iz -1];
            }
        }
    }



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

