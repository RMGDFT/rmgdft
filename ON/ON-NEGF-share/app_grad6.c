/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/


/* calculate the gradient with 6-order */

#include <float.h>
#include <math.h>
#include <stdlib.h>

#include "main.h"
#include "prototypes_on.h"

void app_grad6(double * f, double * wx, double * wy, double * wz, int dimx, int dimy, int dimz)
{

    int iz, ix, iy;
    double t1, t2, t1x, t2x, t1y, t2y, t1z, t2z;
    double t3x,t3y,t3z;
    double *rptr;
    double *wxr, *wyr, *wzr;
    int ixs, iys;
    int ix1, iy1;

    ixs = (dimy + 6) * (dimz + 6);
    iys = (dimz + 6);
    ix1 = dimy * dimz;
    iy1 = dimz;

    wxr = (double *) wx;
    wyr = (double *) wy;
    wzr = (double *) wz;

    my_malloc_init( rptr, (dimx + 6) * (dimy + 6) * (dimz + 6), double );


    trade_imagesx (f, rptr, dimx, dimy, dimz, 3, CENTRAL_TRADE);

    switch (get_ibrav_type())
    {

    case CUBIC_PRIMITIVE:
    case ORTHORHOMBIC_PRIMITIVE:
        t1x = 3.0 / ( 4.0 * get_hxxgrid() * get_xside());
        t2x =-3.0 / (20.0 * get_hxxgrid() * get_xside());
        t3x = 1.0 / (60.0 * get_hxxgrid() * get_xside());

        t1y = 3.0 / ( 4.0 * get_hyygrid() * get_yside());
        t2y =-3.0 / (20.0 * get_hyygrid() * get_yside());
        t3y = 1.0 / (60.0 * get_hyygrid() * get_yside());

        t1z = 3.0 / ( 4.0 * get_hzzgrid() * get_zside());
        t2z =-3.0 / (20.0 * get_hzzgrid() * get_zside());
        t3z = 1.0 / (60.0 * get_hzzgrid() * get_zside());

        for (ix = 3; ix < dimx + 3; ix++)
        {

            for (iy = 3; iy < dimy + 3; iy++)
            {

                for (iz = 3; iz < dimz + 3; iz++)
                {

                    wxr[(ix - 3) * ix1 + (iy - 3) * iy1 + iz - 3] =
                        -t3x * rptr[(ix - 3) * ixs + iy * iys + iz] +
                        -t2x * rptr[(ix - 2) * ixs + iy * iys + iz] +
                        -t1x * rptr[(ix - 1) * ixs + iy * iys + iz] +
                         t1x * rptr[(ix + 1) * ixs + iy * iys + iz] +
                         t2x * rptr[(ix + 2) * ixs + iy * iys + iz] +
                         t3x * rptr[(ix + 3) * ixs + iy * iys + iz];

                    wyr[(ix - 3) * ix1 + (iy - 3) * iy1 + iz - 3] =
                        -t3y * rptr[ix * ixs + (iy - 3) * iys + iz] +
                        -t2y * rptr[ix * ixs + (iy - 2) * iys + iz] +
                        -t1y * rptr[ix * ixs + (iy - 1) * iys + iz] +
                         t1y * rptr[ix * ixs + (iy + 1) * iys + iz] +
                         t2y * rptr[ix * ixs + (iy + 2) * iys + iz] +
                         t3y * rptr[ix * ixs + (iy + 3) * iys + iz];

                    wzr[(ix - 3) * ix1 + (iy - 3) * iy1 + iz - 3] =
                        -t3z * rptr[ix * ixs + iy * iys + iz - 3] +
                        -t2z * rptr[ix * ixs + iy * iys + iz - 2] +
                        -t1z * rptr[ix * ixs + iy * iys + iz - 1] +
                         t1z * rptr[ix * ixs + iy * iys + iz + 1] +
                         t2z * rptr[ix * ixs + iy * iys + iz + 2] +
                         t3z * rptr[ix * ixs + iy * iys + iz + 3];



                }               /* end for */
            }                   /* end for */
        }                       /* end for */

        break;

    default:
        error_handler("Lattice type not implemented");
    }                           /* end switch */

    my_free(rptr);



}   


/******/
