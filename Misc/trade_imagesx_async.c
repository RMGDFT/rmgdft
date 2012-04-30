/************************** SVN Revision Information **************************
 **    $Id: trade_imagesx.c 1694 2012-04-09 19:46:45Z ebriggs $    **
******************************************************************************/

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "main.h"

#if HYBRID_MODEL
#include "hybrid.h"
#endif

#include <pthread.h>


/*
 * INPUTS
 * f[dimx*dimy*dimz] - raw array without images. pack_ptos should not be called, this is 
 *                     handled inside the function now
 * OUTPUT
 * w[(dimx+2*images) * (dimy+2*images) * (dimz+2*images)]
 *    This array is completely filled, i.e. the original data is filled in and then 
 *    the image data are added*/




void trade_imagesx_async (REAL * f, REAL * w, int dimx, int dimy, int dimz, int images)
{
    int ix, iy, iz, incx, incy, incx0, incy0, index, tim, ione = 1;
    int ixs, iys, ixs2, iys2, c1, c2, alloc, idx;
    int xlen, ylen, zlen, stop, yzlen;
    int basetag=0, tid;
    MPI_Status mrstatus;
    REAL *frdx1, *frdx2, *frdy1, *frdy2, *frdz1, *frdz2;
    REAL *frdx1n, *frdx2n, *frdy1n, *frdy2n, *frdz1n, *frdz2n;
    REAL *cpsms_r, *cpsps_r, *cmsms_r, *cmsps_r;
    REAL *cpsms_s, *cpsps_s, *cmsms_s, *cmsps_s;
 
#if MD_TIMERS
    REAL time1, time2, time3;
    time1 = my_crtc ();
#endif

#if HYBRID_MODEL
    basetag = get_thread_basetag();
    tid = get_thread_tid();
#endif

    tim = 2 * images;

    incx = (dimy + tim) * (dimz + tim);
    incy = dimz + tim;
    incx0 = dimy * dimz;
    incy0 = dimz;

    zlen = dimx * dimy * images;
    ylen = dimx * images * (dimz + tim);
    xlen = images * (dimy + tim) * (dimz + tim);

    yzlen = images * images * dimx;

    alloc = 2 * (xlen + ylen + zlen);

#if HYBRID_MODEL
    frdx1 = get_trade_mem(alloc, 0);
    frdx2 = get_trade_mem(alloc, 1);
    frdy1 = get_trade_mem(alloc, 2);
    frdy2 = get_trade_mem(alloc, 3);
    frdz1 = get_trade_mem(alloc, 4);
    frdz2 = get_trade_mem(alloc, 5);

    frdx1n = get_trade_mem(alloc, 6);
    frdx2n = get_trade_mem(alloc, 7);
    frdy1n = get_trade_mem(alloc, 8);
    frdy2n = get_trade_mem(alloc, 9);
    frdz1n = get_trade_mem(alloc, 10);
    frdz2n = get_trade_mem(alloc, 11);
#else
    my_malloc(frdx1, alloc, REAL);
    frdx2 = frdx1 + xlen;
    frdy1 = frdx2 + xlen;
    frdy2 = frdy1 + ylen;
    frdz1 = frdy2 + ylen;
    frdz2 = frdz1 + zlen;

    my_malloc(frdx1n, alloc, REAL);
    frdx2n = frdx1n + xlen;
    frdy1n = frdx2n + xlen;
    frdy2n = frdy1n + ylen;
    frdz1n = frdy2n + ylen;
    frdz2n = frdz1n + zlen;
#endif

    my_malloc(cpsms_r, 8 * yzlen, REAL);
    cpsps_r = cpsms_r + yzlen;
    cmsms_r = cpsps_r + yzlen;
    cmsps_r = cmsms_r + yzlen;
    cpsms_s = cmsps_r + yzlen;
    cpsps_s = cpsms_s + yzlen;
    cmsms_s = cpsps_s  + yzlen;
    cmsps_s = cmsms_s  + yzlen;


    // Post the receives for the z planes
    RMG_MPI_trade(frdz2n, dimx*dimy*images, RMG_MPI_IRECV, 0, 0, 1, pct.grid_comm, 0);
    RMG_MPI_trade(frdz1n, dimx*dimy*images, RMG_MPI_IRECV, 0, 0, -1, pct.grid_comm, 1);

    // The y planes
    RMG_MPI_trade(frdy2n, dimx*dimz*images, RMG_MPI_IRECV, 0, 1, 0, pct.grid_comm, 2);
    RMG_MPI_trade(frdy1n, dimx*dimz*images, RMG_MPI_IRECV, 0, -1, 0, pct.grid_comm, 3);

    // And the yz-plane edges
    RMG_MPI_trade(cpsms_r, yzlen, RMG_MPI_IRECV, 0, 1, -1, pct.grid_comm, 6);
    RMG_MPI_trade(cpsps_r, yzlen, RMG_MPI_IRECV, 0, 1, 1, pct.grid_comm, 7);
    RMG_MPI_trade(cmsms_r, yzlen, RMG_MPI_IRECV, 0, -1, -1, pct.grid_comm, 8);
    RMG_MPI_trade(cmsps_r, yzlen, RMG_MPI_IRECV, 0, -1, 1, pct.grid_comm, 9);


    /* Collect the positive z-plane and negative z-planes */
    idx = 0;
    for (ix = 0; ix < dimx; ix++)
    {
        ixs2 = ix * incx0;
        for (iy = 0; iy < dimy; iy++)
        {
            iys2 = ixs2 + iy * incy0;
            for (iz = 0; iz < images; iz++)
            {

                index = iys2 + iz;

                frdz1[idx] = f[index];
                frdz2[idx] = f[index+dimz-images];
                idx++;

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */

    // Send z planes
    RMG_MPI_trade(frdz1, dimx*dimy*images, RMG_MPI_ISEND, 0, 0, -1, pct.grid_comm, 0);
    RMG_MPI_trade(frdz2, dimx*dimy*images, RMG_MPI_ISEND, 0, 0, 1, pct.grid_comm, 1);


    RMG_MPI_wake_manager();

    /* Collect the north and south planes */
    c1 = (dimy - images)*incy0;
    idx = 0;
    for (ix = 0; ix < dimx; ix++)
    {
        ixs2 = ix * incx0;
        for (iy = 0; iy < images; iy++)
        {

            iys2 = ixs2 + iy * incy0;
            for (iz = 0; iz < dimz; iz++)
            {

                index = iys2 + iz;

                frdy1[idx] = f[index];
                frdy2[idx] = f[index + c1];
                idx++;

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */


    // Send the north and south planes
    RMG_MPI_trade(frdy1, dimx*dimz*images, RMG_MPI_ISEND, 0, -1, 0, pct.grid_comm, 2);
    RMG_MPI_trade(frdy2, dimx*dimz*images, RMG_MPI_ISEND, 0, 1, 0, pct.grid_comm, 3);

    /* Load up w with the basic stuff */
    for (ix = 0; ix < dimx; ix++)
    {
        ixs = ix * incx0;
        ixs2 = (ix + images) * incx;

        for (iy = 0; iy < dimy; iy++)
        {
            iys = ixs + iy * incy0;
            iys2 = ixs2 + (iy + images) * incy;

            QMD_scopy (dimz, &f[iys], ione, &w[iys2 + images], ione);

        }                       /* end for */

    }                           /* end for */



    // Pack yz-plane edges
    idx = 0;
    c1 = (dimy-images)*incy0;
    for (ix = 0; ix < dimx; ix++)
    {
        ixs2 = ix * incx0;
        for(iy = 0;iy < images;iy++) {
            iys2 = ixs2 + iy*incy0;
            for(iz = 0;iz < images;iz++) {
                index = iys2 + iz;
                cpsms_s[idx] = f[index + (dimz - images)]; 
                cpsps_s[idx] = f[index];
                cmsms_s[idx] = f[index + (dimz - images) + c1];
                cmsps_s[idx] = f[index + c1];
                idx++;
            }
        }
    }

    // Send yz-plane edges
    RMG_MPI_trade(cpsms_s, yzlen, RMG_MPI_ISEND, 0, -1, 1, pct.grid_comm, 6);
    RMG_MPI_trade(cpsps_s, yzlen, RMG_MPI_ISEND, 0, -1, -1, pct.grid_comm, 7);
    RMG_MPI_trade(cmsms_s, yzlen, RMG_MPI_ISEND, 0, 1, 1, pct.grid_comm, 8);
    RMG_MPI_trade(cmsps_s, yzlen, RMG_MPI_ISEND, 0, 1, -1, pct.grid_comm, 9);



    // Wait for all the sends to finish
    RMG_MPI_trade_waitall();


    // Post receives for x planes
    RMG_MPI_trade(frdx2n, xlen, RMG_MPI_IRECV, 1, 0, 0, pct.grid_comm, 4);
    RMG_MPI_trade(frdx1n, xlen, RMG_MPI_IRECV, -1, 0, 0, pct.grid_comm, 5);


    // Unpack yz-plane edges
    idx = 0;
    c1 = (dimy + images) * incy;
    for (ix = 0; ix < dimx; ix++)
    {
        ixs2 = (ix + images) * incx;
        for(iy = 0;iy < images;iy++) {
            iys2 = ixs2 + iy*incy;
            for(iz = 0;iz < images;iz++) {
                index = iys2 + iz;
                w[index + c1] = cpsms_r[idx];
                w[index + dimz+images + c1] = cpsps_r[idx];
                w[index] =  cmsms_r[idx];
                w[index + dimz+images] = cmsps_r[idx];
                idx++;
            }
        }
    }

    RMG_MPI_wake_manager();

    /* Unpack z-planes */
    c1 = dimz + images;
    idx = 0;
    for (ix = 0; ix < dimx; ix++)
    {
        ixs2 = (ix + images) * incx;
        for (iy = 0; iy < dimy; iy++)
        {
            iys2 = ixs2 + (iy + images) * incy;
            for (iz = 0; iz < images; iz++)
            {

                index = iys2 + iz;
                w[index] = frdz1n[idx];
                w[index + c1] = frdz2n[idx];
                idx++;

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */


    RMG_MPI_wake_manager();

    /* Unpack the north and south planes */
    c1 = (dimy + images) * incy;
    idx = 0;
    for (ix = 0; ix < dimx; ix++)
    {
        ixs2 = (ix + images) * incx;
        for (iy = 0; iy < images; iy++)
        {
            iys2 = ixs2 + iy * incy;
            for (iz = 0; iz < dimz; iz++)
            {
                index = iys2 + iz + images;
                w[index] = frdy1n[idx];
                w[index + c1] = frdy2n[idx];
                idx++;

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */

    RMG_MPI_wake_manager();

    /* Collect the east and west planes */
    c1 = images * incx;
    c2 = dimx * incx;
    idx = 0;
    for (ix = 0; ix < images; ix++)
    {
        ixs2 = ix * incx;
        for (iy = 0; iy < dimy + tim; iy++)
        {
            iys2 = ixs2 + iy * incy;
            for (iz = 0; iz < dimz + tim; iz++)
            {

                index = iys2 + iz;

                frdx1[idx] = w[index + c1];
                frdx2[idx] = w[index + c2];
                idx++;

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */


    RMG_MPI_trade(frdx1, xlen, RMG_MPI_ISEND, -1, 0, 0, pct.grid_comm, 4);
    RMG_MPI_trade(frdx2, xlen, RMG_MPI_ISEND, 1, 0, 0, pct.grid_comm, 5);
    RMG_MPI_trade_waitall();


    /* Unpack them */
    c1 = (dimx + images) * incx;
    idx = 0;
    for (ix = 0; ix < images; ix++)
    {
        ixs2 = ix * incx;
        for (iy = 0; iy < dimy + tim; iy++)
        {
            iys2 = ixs2 + iy * incy;
            for (iz = 0; iz < dimz + tim; iz++)
            {

                index = iys2 + iz;

                w[index] = frdx1n[idx];
                w[index + c1] = frdx2n[idx];
                idx++;

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */

my_free(cpsms_r);

#if HYBRID_MODEL
    free_trade_mem (frdx1n);
    free_trade_mem (frdx1);
    free_trade_mem (frdy1n);
    free_trade_mem (frdy1);
    free_trade_mem (frdz1n);
    free_trade_mem (frdz1);

    free_trade_mem (frdx2n);
    free_trade_mem (frdx2);
    free_trade_mem (frdy2n);
    free_trade_mem (frdy2);
    free_trade_mem (frdz2n);
    free_trade_mem (frdz2);
#else
    my_free(frdx1n);
    my_free(frdx1);
#endif

#if MD_TIMERS
    time2 = my_crtc ();
    rmg_timings (IMAGE_TIME, (time2 - time1));
#endif

} // end trade_imagesx
