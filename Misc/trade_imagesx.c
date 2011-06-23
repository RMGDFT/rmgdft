/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "main.h"

#if HYBRID_MODEL
#include "hybrid.h"
#endif

/*
 * INPUTS
 * f[dimx*dimy*dimz] - raw array without images. pack_ptos should not be called, this is 
 *                     handled inside the function now
 * OUTPUT
 * w[(dimx+2*images) * (dimy+2*images) * (dimz+2*images)]
 *    This array is completely filled, i.e. the original data is filled in and then 
 *    the image data are added*/




void trade_imagesx (REAL * f, REAL * w, int dimx, int dimy, int dimz, int images)
{
    int ix, iy, iz, incx, incy, incx0, incy0, index, tim, ione = 1;
    int ixs, iys, ixs2, iys2, c1, c2;
    int xlen, ylen, zlen;
    int *nb_ids, basetag=0;
    MPI_Status mrstatus;
    REAL *frdx1, *frdx2, *frdy1, *frdy2, *frdz1, *frdz2;
    REAL *frdx1n, *frdx2n, *frdy1n, *frdy2n, *frdz1n, *frdz2n;

#if MD_TIMERS
    REAL time1, time2;
    time1 = my_crtc ();
#endif

#if HYBRID_MODEL
    basetag = get_thread_basetag();
#endif

    tim = 2 * images;

    incx = (dimy + tim) * (dimz + tim);
    incy = dimz + tim;
    incx0 = dimy * dimz;
    incy0 = dimz;

    zlen = dimx * dimy * images;
    ylen = dimx * images * (dimz + tim);
    xlen = images * (dimy + tim) * (dimz + tim);



    my_malloc (frdx1, 2 * (xlen + ylen + zlen), REAL);
    frdx2 = frdx1 + xlen;
    frdy1 = frdx2 + xlen;
    frdy2 = frdy1 + ylen;
    frdz1 = frdy2 + ylen;
    frdz2 = frdz1 + zlen;

    my_malloc (frdx1n, 2 * (xlen + ylen + zlen), REAL);
    frdx2n = frdx1n + xlen;
    frdy1n = frdx2n + xlen;
    frdy2n = frdy1n + ylen;
    frdz1n = frdy2n + ylen;
    frdz2n = frdz1n + zlen;


    nb_ids = &pct.neighbors[0];

    /* Load up w with the basic stuff */
    for (ix = 0; ix < dimx; ix++)
    {
        ixs = ix * incx0;
        ixs2 = (ix + images) * incx;

        for (iy = 0; iy < dimy; iy++)
        {
            iys = ixs + iy * incy0;
            iys2 = ixs2 + (iy + images) * incy;

            scopy (&dimz, &f[iys], &ione, &w[iys2 + images], &ione);

        }                       /* end for */

    }                           /* end for */


    /* Collect the positive z-plane and negative z-planes */
    for (ix = 0; ix < dimx; ix++)
    {
        ixs = ix * dimy * images;
        ixs2 = (ix + images) * incx;
        for (iy = 0; iy < dimy; iy++)
        {
            iys = ixs + iy * images;
            iys2 = ixs2 + (iy + images) * incy;
            for (iz = 0; iz < images; iz++)
            {

                index = iys2 + iz;

                frdz1[iys + iz] = w[index + images];
                frdz2[iys + iz] = w[index + dimz];

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */


    if(!is_loop_over_states() || !HYBRID_MODEL) {

        MPI_Sendrecv (frdz1, zlen, MPI_DOUBLE, nb_ids[NB_D], basetag + (1>>16),
                  frdz2n, zlen, MPI_DOUBLE, nb_ids[NB_U], basetag + (1>>16), pct.grid_comm, &mrstatus);

        MPI_Sendrecv (frdz2, zlen, MPI_DOUBLE, nb_ids[NB_U], basetag + (2>>16),
                  frdz1n, zlen, MPI_DOUBLE, nb_ids[NB_D], basetag + (2>>16), pct.grid_comm, &mrstatus);

    }
    else {

        RMG_MPI_thread_order_lock();
        MPI_Sendrecv (frdz1, zlen, MPI_DOUBLE, nb_ids[NB_D], basetag + (1>>16),
                  frdz2n, zlen, MPI_DOUBLE, nb_ids[NB_U], basetag + (1>>16), pct.grid_comm, &mrstatus);

        MPI_Sendrecv (frdz2, zlen, MPI_DOUBLE, nb_ids[NB_U], basetag + (2>>16),
                  frdz1n, zlen, MPI_DOUBLE, nb_ids[NB_D], basetag + (2>>16), pct.grid_comm, &mrstatus);
        RMG_MPI_thread_order_unlock();

    }
 

    /* Unpack them */
    c1 = dimz + images;
    for (ix = 0; ix < dimx; ix++)
    {

        ixs = ix * dimy * images;
        ixs2 = (ix + images) * incx;
        for (iy = 0; iy < dimy; iy++)
        {
            iys = ixs + iy * images;
            iys2 = ixs2 + (iy + images) * incy;
            for (iz = 0; iz < images; iz++)
            {

                index = iys2 + iz;

                w[index] = frdz1n[iys + iz];
                w[index + c1] = frdz2n[iys + iz];

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */



    /* Collect the north and south planes */
    c1 = images * incy;
    c2 = dimy * incy;
    for (ix = 0; ix < dimx; ix++)
    {
        ixs = ix * images * (dimz + tim);
        ixs2 = (ix + images) * incx;
        for (iy = 0; iy < images; iy++)
        {

            iys = ixs + iy * (dimz + tim);
            iys2 = ixs2 + iy * incy;
            for (iz = 0; iz < dimz + tim; iz++)
            {

                index = iys2 + iz;

                frdy1[iys + iz] = w[index + c1];
                frdy2[iys + iz] = w[index + c2];

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */


    if(!is_loop_over_states() || !HYBRID_MODEL) {

        MPI_Sendrecv (frdy1, ylen, MPI_DOUBLE, nb_ids[NB_S], basetag + (3>>16),
                  frdy2n, ylen, MPI_DOUBLE, nb_ids[NB_N], basetag + (3>>16), pct.grid_comm, &mrstatus);

        MPI_Sendrecv (frdy2, ylen, MPI_DOUBLE, nb_ids[NB_N], basetag + (4>>16),
                  frdy1n, ylen, MPI_DOUBLE, nb_ids[NB_S], basetag + (4>>16), pct.grid_comm, &mrstatus);

     }
     else {

        RMG_MPI_thread_order_lock();
        MPI_Sendrecv (frdy1, ylen, MPI_DOUBLE, nb_ids[NB_S], basetag + (3>>16),
                  frdy2n, ylen, MPI_DOUBLE, nb_ids[NB_N], basetag + (3>>16), pct.grid_comm, &mrstatus);

        MPI_Sendrecv (frdy2, ylen, MPI_DOUBLE, nb_ids[NB_N], basetag + (4>>16),
                  frdy1n, ylen, MPI_DOUBLE, nb_ids[NB_S], basetag + (4>>16), pct.grid_comm, &mrstatus);
        RMG_MPI_thread_order_unlock();

     }


    /* Unpack them */
    c1 = (dimy + images) * incy;
    for (ix = 0; ix < dimx; ix++)
    {
        ixs = ix * images * (dimz + tim);
        ixs2 = (ix + images) * incx;
        for (iy = 0; iy < images; iy++)
        {
            iys = ixs + iy * (dimz + tim);
            iys2 = ixs2 + iy * incy;
            for (iz = 0; iz < dimz + tim; iz++)
            {

                index = iys2 + iz;

                w[index] = frdy1n[iys + iz];
                w[index + c1] = frdy2n[iys + iz];

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */


    /* Collect the east and west planes */
    c1 = images * incx;
    c2 = dimx * incx;
    for (ix = 0; ix < images; ix++)
    {
        ixs = ix * (dimy + tim) * (dimz + tim);
        ixs2 = ix * incx;
        for (iy = 0; iy < dimy + tim; iy++)
        {
            iys = ixs + iy * (dimz + tim);
            iys2 = ixs2 + iy * incy;
            for (iz = 0; iz < dimz + tim; iz++)
            {

                index = iys2 + iz;

                frdx1[iys + iz] = w[index + c1];
                frdx2[iys + iz] = w[index + c2];

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */


    if(!is_loop_over_states() || !HYBRID_MODEL) {

        MPI_Sendrecv (frdx1, xlen, MPI_DOUBLE, nb_ids[NB_W], basetag + (5>>16),
                  frdx2n, xlen, MPI_DOUBLE, nb_ids[NB_E], basetag + (5>>16), pct.grid_comm, &mrstatus);

        MPI_Sendrecv (frdx2, xlen, MPI_DOUBLE, nb_ids[NB_E], basetag + (6>>16),
                  frdx1n, xlen, MPI_DOUBLE, nb_ids[NB_W], basetag + (6>>16), pct.grid_comm, &mrstatus);

    }
    else {

        RMG_MPI_thread_order_lock();
        MPI_Sendrecv (frdx1, xlen, MPI_DOUBLE, nb_ids[NB_W], basetag + (5>>16),
                  frdx2n, xlen, MPI_DOUBLE, nb_ids[NB_E], basetag + (5>>16), pct.grid_comm, &mrstatus);

        MPI_Sendrecv (frdx2, xlen, MPI_DOUBLE, nb_ids[NB_E], basetag + (6>>16),
                  frdx1n, xlen, MPI_DOUBLE, nb_ids[NB_W], basetag + (6>>16), pct.grid_comm, &mrstatus);
        RMG_MPI_thread_order_unlock();
    }

    /* Unpack them */
    c1 = (dimx + images) * incx;
    for (ix = 0; ix < images; ix++)
    {
        ixs = ix * (dimy + tim) * (dimz + tim);
        ixs2 = ix * incx;
        for (iy = 0; iy < dimy + tim; iy++)
        {
            iys = ixs + iy * (dimz + tim);
            iys2 = ixs2 + iy * incy;
            for (iz = 0; iz < dimz + tim; iz++)
            {

                index = iys2 + iz;

                w[index] = frdx1n[iys + iz];
                w[index + c1] = frdx2n[iys + iz];

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */


    my_free (frdx1n);
    my_free (frdx1);

#if MD_TIMERS
    time2 = my_crtc ();
    rmg_timings (IMAGE_TIME, (time2 - time1));
#endif

}                               /* end trade_images2 */
