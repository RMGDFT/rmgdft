/************************** SVN Revision Information **************************
 **    $Id: app_cir_sixth.c 1722 2012-05-07 19:38:37Z ebriggs $    **
******************************************************************************/

#include "const.h"
#include "common_prototypes.h"
#include "main.h"
#include "rmg_alloc.h"
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include <TradeImages.h>


void app_cir_beta_fourth (double * a, double * b, int dimx, int dimy, int dimz)
{

    int ix, iy, iz, ibrav;
    int ixs, iys, ixms, ixps, iyms, iyps;
    int incy, incx;
    int incyr, incxr;
    double *rptr;
    double c000, c100, Bc, Bf, Bz;

    ibrav = get_ibrav_type();

    incx = (dimz + 2) * (dimy + 2);
    incy = dimz + 2;
    incxr = dimz * dimy;
    incyr = dimz;

    my_malloc (rptr, (dimx + 2) * (dimy + 2) * (dimz + 2), double);

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


    switch(ibrav) {
        case CUBIC_PRIMITIVE:
        case ORTHORHOMBIC_PRIMITIVE:

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

            break;

       case CUBIC_FC: 

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
                            Bc * rptr[ixs + iys + iz] +
                            Bf * (rptr[ixms + iys + iz] +
                                  rptr[ixms + iys + iz + 1] +
                                  rptr[ixms + iyps + iz] +
                                  rptr[ixs + iyms + iz] +
                                  rptr[ixs + iyms + iz + 1] +
                                  rptr[ixs + iys + iz - 1] +
                                  rptr[ixs + iys + iz + 1] +
                                  rptr[ixs + iyps + iz - 1] +
                                  rptr[ixs + iyps + iz] +
                                  rptr[ixps + iyms + iz] + 
                                  rptr[ixps + iys + iz - 1] + 
                                  rptr[ixps + iys + iz]);


                    }                   /* end for */

                }                       /* end for */

            }                           /* end for */

            break;

       case HEXAGONAL:
            Bc = 7.0 / 12.0;
            Bf = 1.0 / 24.0;
            Bz = 1.0 / 12.0;
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
                            Bc * rptr[ixs + iys + iz] +
                            Bz * rptr[ixs + iys + iz - 1] +
                            Bz * rptr[ixs + iys + iz + 1] +
                            Bf * rptr[ixps + iys + iz] +
                            Bf * rptr[ixps + iyms + iz] +
                            Bf * rptr[ixs + iyms + iz] +
                            Bf * rptr[ixms + iys + iz] +
                            Bf * rptr[ixms + iyps + iz] +
                            Bf * rptr[ixs + iyps + iz];


                    }                   /* end for */

                }                       /* end for */

            }                           /* end for */
            break;

       default:
           error_handler("Grid symmetry not programmed yet in app_cir_fourth.\n");

    } // end switch

    my_free (rptr);
}

