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
    int ix, iy, iz, ix1, iy1, iz1, incx, incy, incx0, incy0, index, tim, ione = 1;
    int ixs, iys, ixs2, iys2, c1, c2, c3, alloc, idx, idx1, img3;
    int xlen, ylen, zlen, stop, yzlen, xylen, xzlen;
    int basetag=0, tid;
    MPI_Status mrstatus;
    REAL *frdx1, *frdx2, *frdy1, *frdy2, *frdz1, *frdz2;
    REAL *frdx1n, *frdx2n, *frdy1n, *frdy2n, *frdz1n, *frdz2n;
    REAL *yzpsms_r, *yzpsps_r, *yzmsms_r, *yzmsps_r;
    REAL *yzpsms_s, *yzpsps_s, *yzmsms_s, *yzmsps_s;
    REAL *xzpsms_r, *xzpsps_r, *xzmsms_r, *xzmsps_r;
    REAL *xzpsms_s, *xzpsps_s, *xzmsms_s, *xzmsps_s;
    REAL *xypsms_r, *xypsps_r, *xymsms_r, *xymsps_r;
    REAL *xypsms_s, *xypsps_s, *xymsms_s, *xymsps_s;

//  Need to be dimensioned [8][images*images*images] so 64 is good for 4th order
    REAL m0_s[8][64];
    REAL m0_r[8][64];

#if MD_TIMERS
    REAL time1, time2, time3;
    time1 = my_crtc ();
#endif

#if HYBRID_MODEL
    basetag = get_thread_basetag();
    tid = get_thread_tid();
#endif

    tim = 2 * images;
    img3 = images*images*images;

    incx = (dimy + tim) * (dimz + tim);
    incy = dimz + tim;
    incx0 = dimy * dimz;
    incy0 = dimz;

    zlen = dimx * dimy * images;
    ylen = dimx * images * (dimz + tim);
    xlen = images * (dimy + tim) * (dimz + tim);

    yzlen = images * images * dimx;
    xylen = images * images * dimz;
    xzlen = images * images * dimy;

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

    my_malloc(yzpsms_r, 8 * yzlen, REAL);
    yzpsps_r = yzpsms_r + yzlen;
    yzmsms_r = yzpsps_r + yzlen;
    yzmsps_r = yzmsms_r + yzlen;
    yzpsms_s = yzmsps_r + yzlen;
    yzpsps_s = yzpsms_s + yzlen;
    yzmsms_s = yzpsps_s + yzlen;
    yzmsps_s = yzmsms_s + yzlen;

    my_malloc(xypsms_r, 8 * xylen, REAL);
    xypsps_r = xypsms_r + xylen;
    xymsms_r = xypsps_r + xylen;
    xymsps_r = xymsms_r + xylen;
    xypsms_s = xymsps_r + xylen;
    xypsps_s = xypsms_s + xylen;
    xymsms_s = xypsps_s + xylen;
    xymsps_s = xymsms_s + xylen;

    my_malloc(xzpsms_r, 8 * xzlen, REAL);
    xzpsps_r = xzpsms_r + xzlen;
    xzmsms_r = xzpsps_r + xzlen;
    xzmsps_r = xzmsms_r + xzlen;
    xzpsms_s = xzmsps_r + xzlen;
    xzpsps_s = xzpsms_s + xzlen;
    xzmsms_s = xzpsps_s + xzlen;
    xzmsps_s = xzmsms_s + xzlen;


    // Post receives for the corner nodes
    RMG_MPI_trade( &m0_r[0][0], img3, RMG_MPI_IRECV, -1, -1, -1, pct.grid_comm, 18);
    RMG_MPI_trade( &m0_r[1][0], img3, RMG_MPI_IRECV, -1, -1, 1, pct.grid_comm, 19);
    RMG_MPI_trade( &m0_r[2][0], img3, RMG_MPI_IRECV, -1, 1, -1, pct.grid_comm, 20);
    RMG_MPI_trade( &m0_r[3][0], img3, RMG_MPI_IRECV, -1, 1, 1, pct.grid_comm, 21);
    RMG_MPI_trade( &m0_r[4][0], img3, RMG_MPI_IRECV, 1, -1, -1, pct.grid_comm, 22);
    RMG_MPI_trade( &m0_r[5][0], img3, RMG_MPI_IRECV, 1, -1, 1, pct.grid_comm, 23);
    RMG_MPI_trade( &m0_r[6][0], img3, RMG_MPI_IRECV, 1, 1, -1, pct.grid_comm, 24);
    RMG_MPI_trade( &m0_r[7][0], img3, RMG_MPI_IRECV, 1, 1, 1, pct.grid_comm, 25);

    // Post the receives for the z planes
    RMG_MPI_trade(frdz2n, dimx*dimy*images, RMG_MPI_IRECV, 0, 0, 1, pct.grid_comm, 0);
    RMG_MPI_trade(frdz1n, dimx*dimy*images, RMG_MPI_IRECV, 0, 0, -1, pct.grid_comm, 1);

    // The y planes
    RMG_MPI_trade(frdy2n, dimx*dimz*images, RMG_MPI_IRECV, 0, 1, 0, pct.grid_comm, 2);
    RMG_MPI_trade(frdy1n, dimx*dimz*images, RMG_MPI_IRECV, 0, -1, 0, pct.grid_comm, 3);

    // The x planes
    RMG_MPI_trade(frdx2n, dimy*dimz*images, RMG_MPI_IRECV, 1, 0, 0, pct.grid_comm, 4);
    RMG_MPI_trade(frdx1n, dimy*dimz*images, RMG_MPI_IRECV, -1, 0, 0, pct.grid_comm, 5);

    // And the yz-plane edges
    RMG_MPI_trade(yzpsms_r, yzlen, RMG_MPI_IRECV, 0, 1, -1, pct.grid_comm, 6);
    RMG_MPI_trade(yzpsps_r, yzlen, RMG_MPI_IRECV, 0, 1, 1, pct.grid_comm, 7);
    RMG_MPI_trade(yzmsms_r, yzlen, RMG_MPI_IRECV, 0, -1, -1, pct.grid_comm, 8);
    RMG_MPI_trade(yzmsps_r, yzlen, RMG_MPI_IRECV, 0, -1, 1, pct.grid_comm, 9);

    // And the xy-plane edges
    RMG_MPI_trade(xypsms_r, xylen, RMG_MPI_IRECV, 1, -1, 0, pct.grid_comm, 10);
    RMG_MPI_trade(xypsps_r, xylen, RMG_MPI_IRECV, 1, 1, 0, pct.grid_comm, 11);
    RMG_MPI_trade(xymsms_r, xylen, RMG_MPI_IRECV, -1, -1, 0, pct.grid_comm, 12);
    RMG_MPI_trade(xymsps_r, xylen, RMG_MPI_IRECV, -1, 1, 0, pct.grid_comm, 13);

    // And the xz-plane edges
    RMG_MPI_trade(xzpsms_r, xzlen, RMG_MPI_IRECV, 1, 0, -1, pct.grid_comm, 14);
    RMG_MPI_trade(xzpsps_r, xzlen, RMG_MPI_IRECV, 1, 0, 1, pct.grid_comm, 15);
    RMG_MPI_trade(xzmsms_r, xzlen, RMG_MPI_IRECV, -1, 0, -1, pct.grid_comm, 16);
    RMG_MPI_trade(xzmsps_r, xzlen, RMG_MPI_IRECV, -1, 0, 1, pct.grid_comm, 17);


    // Collect data for the corner nodes
    idx = 0;
    for(ix = -1;ix <= 1;ix+=2) 
    {
        c1 = (ix == -1) ? 0 : (dimx-images)*incx0;

        for(iy = -1;iy <= 1;iy+=2) 
        {
            c2 = (iy == -1) ? 0 : (dimy-images)*incy0;

            for(iz = -1;iz <= 1;iz+=2) 
            {
                c3 = (iz == -1) ? 0 : (dimz-images);

                // Pack the send arrays
                idx1 = 0;
                for(ix1 = 0;ix1 < images;ix1++) 
                {
                    ixs2 = ix1*incx0;

                    for(iy1 = 0;iy1 < images;iy1++) 
                    {
                        iys2 = ixs2 + iy1*incy0;

                        for(iz1 = 0;iz1 < images;iz1++) 
                        {
                            index = iys2 + iz1;
                            m0_s[idx][idx1] = f[index + c1 + c2 + c3];
                            idx1++;
                        }
                    }
                }

                idx++;
            }
        }
    }

    /* Collect the positive z-plane and negative z-planes */
    c1 = (dimz-images);
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
                frdz2[idx] = f[index+c1];
                idx++;

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */


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



    /* Collect the east and west planes */
    c1 = (dimx - images)*incx0;
    idx = 0;
    for (ix = 0; ix < images; ix++)
    {
        ixs2 = ix * incx0;
        for (iy = 0; iy < dimy; iy++)
        {
            iys2 = ixs2 + iy * incy0;
            for (iz = 0; iz < dimz; iz++)
            {

                index = iys2 + iz;

                frdx1[idx] = f[index];
                frdx2[idx] = f[index + c1];
                idx++;

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */



    // Pack yz-plane edges
    idx = 0;
    c1 = (dimy-images)*incy0;
    c2 = (dimz - images);
    for (ix = 0; ix < dimx; ix++)
    {
        ixs2 = ix * incx0;
        for(iy = 0;iy < images;iy++) {
            iys2 = ixs2 + iy*incy0;
            for(iz = 0;iz < images;iz++) {
                index = iys2 + iz;
                yzpsms_s[idx] = f[index + c2]; 
                yzpsps_s[idx] = f[index];
                yzmsms_s[idx] = f[index + c2 + c1];
                yzmsps_s[idx] = f[index + c1];
                idx++;
            }
        }
    }


    // Pack xy plane edges
    idx = 0;
    c1 = (dimy-images)*incy0;
    c2 = (dimx-images)*incx0;
    for (ix = 0; ix < images; ix++)
    {
        ixs2 = ix * incx0;
        for(iy = 0;iy < images;iy++) {
            iys2 = ixs2 + iy*incy0;
            for(iz = 0;iz < dimz;iz++) {
                index = iys2 + iz;
                xypsms_s[idx] = f[index + c2]; 
                xypsps_s[idx] = f[index];
                xymsms_s[idx] = f[index + c2 + c1];
                xymsps_s[idx] = f[index + c1];
                idx++;
            }
        }
    }


    // Pack xz plane edges
    idx = 0;
    c1 = (dimz-images);
    c2 = (dimx-images)*incx0;
    for (ix = 0; ix < images; ix++)
    {
        ixs2 = ix * incx0;
        for(iy = 0;iy < dimy;iy++) {
            iys2 = ixs2 + iy*incy0;
            for(iz = 0;iz < images;iz++) {
                index = iys2 + iz;
                xzpsms_s[idx] = f[index + c2]; 
                xzpsps_s[idx] = f[index];
                xzmsms_s[idx] = f[index + c2 + c1];
                xzmsps_s[idx] = f[index + c1];
                idx++;
            }
        }
    }


    // Send the corners
    RMG_MPI_trade( &m0_s[0][0], img3, RMG_MPI_ISEND, 1, 1, 1, pct.grid_comm, 18);
    RMG_MPI_trade( &m0_s[1][0], img3, RMG_MPI_ISEND, 1, 1, -1, pct.grid_comm, 19);
    RMG_MPI_trade( &m0_s[2][0], img3, RMG_MPI_ISEND, 1, -1, 1, pct.grid_comm, 20);
    RMG_MPI_trade( &m0_s[3][0], img3, RMG_MPI_ISEND, 1, -1, -1, pct.grid_comm, 21);
    RMG_MPI_trade( &m0_s[4][0], img3, RMG_MPI_ISEND, -1, 1, 1, pct.grid_comm, 22);
    RMG_MPI_trade( &m0_s[5][0], img3, RMG_MPI_ISEND, -1, 1, -1, pct.grid_comm, 23);
    RMG_MPI_trade( &m0_s[6][0], img3, RMG_MPI_ISEND, -1, -1, 1, pct.grid_comm, 24);
    RMG_MPI_trade( &m0_s[7][0], img3, RMG_MPI_ISEND, -1, -1, -1, pct.grid_comm, 25);

    // Send z planes
    RMG_MPI_trade(frdz1, dimx*dimy*images, RMG_MPI_ISEND, 0, 0, -1, pct.grid_comm, 0);
    RMG_MPI_trade(frdz2, dimx*dimy*images, RMG_MPI_ISEND, 0, 0, 1, pct.grid_comm, 1);

    // Send the north and south planes
    RMG_MPI_trade(frdy1, dimx*dimz*images, RMG_MPI_ISEND, 0, -1, 0, pct.grid_comm, 2);
    RMG_MPI_trade(frdy2, dimx*dimz*images, RMG_MPI_ISEND, 0, 1, 0, pct.grid_comm, 3);

    // Send the east and west planes
    RMG_MPI_trade(frdx1, dimy*dimz*images, RMG_MPI_ISEND, -1, 0, 0, pct.grid_comm, 4);
    RMG_MPI_trade(frdx2, dimy*dimz*images, RMG_MPI_ISEND, 1, 0, 0, pct.grid_comm, 5);

    // Send yz-plane edges
    RMG_MPI_trade(yzpsms_s, yzlen, RMG_MPI_ISEND, 0, -1, 1, pct.grid_comm, 6);
    RMG_MPI_trade(yzpsps_s, yzlen, RMG_MPI_ISEND, 0, -1, -1, pct.grid_comm, 7);
    RMG_MPI_trade(yzmsms_s, yzlen, RMG_MPI_ISEND, 0, 1, 1, pct.grid_comm, 8);
    RMG_MPI_trade(yzmsps_s, yzlen, RMG_MPI_ISEND, 0, 1, -1, pct.grid_comm, 9);

    // Send xy-plane edges
    RMG_MPI_trade(xypsms_s, xylen, RMG_MPI_ISEND, -1, 1, 0, pct.grid_comm, 10);
    RMG_MPI_trade(xypsps_s, xylen, RMG_MPI_ISEND, -1, -1, 0, pct.grid_comm, 11);
    RMG_MPI_trade(xymsms_s, xylen, RMG_MPI_ISEND, 1, 1, 0, pct.grid_comm, 12);
    RMG_MPI_trade(xymsps_s, xylen, RMG_MPI_ISEND, 1, -1, 0, pct.grid_comm, 13);

    // Send xz-plane edges
    RMG_MPI_trade(xzpsms_s, xzlen, RMG_MPI_ISEND, -1, 0, 1, pct.grid_comm, 14);
    RMG_MPI_trade(xzpsps_s, xzlen, RMG_MPI_ISEND, -1, 0, -1, pct.grid_comm, 15);
    RMG_MPI_trade(xzmsms_s, xzlen, RMG_MPI_ISEND, 1, 0, 1, pct.grid_comm, 16);
    RMG_MPI_trade(xzmsps_s, xzlen, RMG_MPI_ISEND, 1, 0, -1, pct.grid_comm, 17);


    // Load up w with the interior points
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


    // Wait for all the recvs to finish
    RMG_MPI_recv_waitall();

    // Unpack the corners
    idx = 0;
    for(ix = -1;ix <= 1;ix+=2)
    {
        c1 = (ix == -1) ? (dimx+images)*incx : 0;

        for(iy = -1;iy <= 1;iy+=2)
        {
            c2 = (iy == -1) ? (dimy+images)*incy : 0;

            for(iz = -1;iz <= 1;iz+=2)
            {
                c3 = (iz == -1) ? (dimz+images) : 0;

                // unpack the recv arrays
                idx1 = 0;
                for(ix1 = 0;ix1 < images;ix1++)
                {
                    ixs2 = ix1*incx;

                    for(iy1 = 0;iy1 < images;iy1++)
                    {
                        iys2 = ixs2 + iy1*incy;

                        for(iz1 = 0;iz1 < images;iz1++)
                        {
                            index = iys2 + iz1;
                            w[index + c1 + c2 + c3] = m0_r[idx][idx1];
                            idx1++;
                        }
                    }
                }
                idx++;

            }
        }
    }


    // Unpack yz-plane edges
    idx = 0;
    c1 = (dimy + images) * incy;
    c2 = (dimz + images);
    for (ix = 0; ix < dimx; ix++)
    {
        ixs2 = (ix + images) * incx;
        for(iy = 0;iy < images;iy++) {
            iys2 = ixs2 + iy*incy;
            for(iz = 0;iz < images;iz++) {
                index = iys2 + iz;
                w[index + c1] = yzpsms_r[idx];
                w[index + c2 + c1] = yzpsps_r[idx];
                w[index] =  yzmsms_r[idx];
                w[index + c2] = yzmsps_r[idx];
                idx++;
            }
        }
    }


    // Unpack xy-plane edges
    idx = 0;
    c1 = (dimy + images) * incy;
    c2 = (dimx + images) * incx;
    for (ix = 0; ix < images; ix++)
    {
        ixs2 = ix * incx;
        for(iy = 0;iy < images;iy++) {
            iys2 = ixs2 + iy*incy;
            for(iz = 0;iz < dimz;iz++) {
                index = iys2 + iz + images;
                w[index + c1] = xypsms_r[idx];
                w[index + c2 + c1] = xypsps_r[idx];
                w[index] =  xymsms_r[idx];
                w[index + c2] = xymsps_r[idx];
                idx++;
            }
        }
    }


    // Unpack xz-plane edges
    idx = 0;
    c1 = (dimz + images);
    c2 = (dimx + images) * incx;
    for (ix = 0; ix < images; ix++)
    {
        ixs2 = ix * incx;
        for(iy = 0;iy < dimy;iy++) {
            iys2 = ixs2 + (iy + images)*incy;
            for(iz = 0;iz < images;iz++) {
                index = iys2 + iz;
                w[index + c1] = xzpsms_r[idx];
                w[index + c2 + c1] = xzpsps_r[idx];
                w[index] =  xzmsms_r[idx];
                w[index + c2] = xzmsps_r[idx];
                idx++;
            }
        }
    }


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


    /* Unpack the east and west planes */
    c1 = (dimx + images) * incx;
    idx = 0;
    for (ix = 0; ix < images; ix++)
    {
        ixs2 = ix * incx;
        for (iy = 0; iy < dimy; iy++)
        {
            iys2 = ixs2 + (iy + images) * incy;
            for (iz = 0; iz < dimz; iz++)
            {

                index = iys2 + iz + images;
                w[index] = frdx1n[idx];
                w[index + c1] = frdx2n[idx];
                idx++;

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */

    // Finally wait for all the sends to finish
    RMG_MPI_send_waitall();

    my_free(xzpsms_r);
    my_free(xypsms_r);
    my_free(yzpsms_r);

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
