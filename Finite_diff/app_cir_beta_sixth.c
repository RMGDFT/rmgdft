/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

#include "const.h"
#include "common_prototypes.h"
#include "rmg_alloc.h"


#include <float.h>
#include <math.h>
#include <stdlib.h>


// This is B operator applied on projector |beta>
// no trade_images is needed.

void app_cir_beta_sixth (double * a, double * b, int dimx, int dimy, int dimz)
{

    int used_alloc=FALSE;
    double *rptr=NULL;


    int ix, iy, iz;
    int ixs, iys, ixms, ixps, iyms, iyps;
    int incy, incx, ixmms, ixpps, iymms, iypps;
    int incyr, incxr;
    double c000, c100, c110, c200;

    incx = (dimz + 4) * (dimy + 4);
    incy = dimz + 4;
    incxr = dimz * dimy;
    incyr = dimz;

    int sbasis = (dimx + 4) * (dimy + 4) * (dimz + 4);


    // If rptr is null then we must allocate it here
    if(rptr == NULL) {
        my_malloc (rptr, sbasis + 64, double);
        used_alloc = TRUE;
    }


    for(ix = 0; ix < (dimx + 4) * (dimy + 4) * (dimz + 4); ix++)
        rptr[ix] = 0.0;
    for (ix = 2; ix < dimx + 2; ix++)
    {
        ixs = ix * incx;

        for (iy = 2; iy < dimy + 2; iy++)
        {
            iys = iy * incy;

            for (iz = 2; iz < dimz + 2; iz++)
            {

                rptr[ixs + iys + iz] = a[(ix-2)* incxr + (iy-2) * incyr + iz -2];
            }
        }
    }




    c000 = 61.0 / 120.0;
    c100 = 13.0 / 180.0;
    c110 = 1.0 / 144.0;
    c200 = -1.0 / 240.0;
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

                b[(ix - 2) * incxr + (iy - 2) * incyr + (iz - 2)] =
                    c100 * (rptr[ixs + iys + (iz - 1)] +
                            rptr[ixs + iys + (iz + 1)] +
                            rptr[ixms + iys + iz] +
                            rptr[ixps + iys + iz] +
                            rptr[ixs + iyms + iz] +
                            rptr[ixs + iyps + iz]) + c000 * rptr[ixs + iys + iz];

                b[(ix - 2) * incxr + (iy - 2) * incyr + (iz - 2)] +=
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

                b[(ix - 2) * incxr + (iy - 2) * incyr + (iz - 2)] +=
                    c200 * (rptr[ixs + iys + (iz - 2)] +
                            rptr[ixs + iys + (iz + 2)] +
                            rptr[ixmms + iys + iz] +
                            rptr[ixpps + iys + iz] +
                            rptr[ixs + iymms + iz] + rptr[ixs + iypps + iz]);

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */

    if(used_alloc)
        my_free(rptr);

}

