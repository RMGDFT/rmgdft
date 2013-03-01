/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <ieee754.h>
#include <stdint.h>
#include "main.h"

#if HYBRID_MODEL
#include "hybrid.h"
#endif

#include <pthread.h>

static REAL *swbuf1x = NULL;
static REAL *swbuf2x = NULL;
static int max_alloc;

/*
 * INPUTS
 * f[dimx*dimy*dimz] - raw array without images. pack_ptos should not be called, this is 
 *                     handled inside the function now
 
 * w[(dimx+2*images) * (dimy+2*images) * (dimz+2*images)]
 *    This array is completely filled, i.e. the original data is filled in and then 
 *    the image data are added*/




void trade_imagesx (REAL * f, REAL * w, int dimx, int dimy, int dimz, int images, int type)
{
    int ix, iy, iz, incx, incy, incx0, incy0, index, tim, ione = 1, retval;
    int ixs, iys, ixs2, iys2, c1, c2, alloc;
    int xlen, ylen, zlen, stop, tid, grid_max1, grid_max2;
    int *nb_ids, basetag=0;
    MPI_Status mrstatus;
    REAL *frdx1, *frdx2, *frdy1, *frdy2, *frdz1, *frdz2;
    REAL *frdx1n, *frdx2n, *frdy1n, *frdy2n, *frdz1n, *frdz2n;
    int ACTIVE_THREADS = 1;

    // To initialize call from init.c with NULL args
    if(swbuf1x == NULL) {
        int grid_xp, grid_yp, grid_zp;

        grid_xp = pct.FPX0_GRID + 2*MAX_TRADE_IMAGES;
        grid_yp = pct.FPY0_GRID + 2*MAX_TRADE_IMAGES;
        grid_zp = pct.FPZ0_GRID + 2*MAX_TRADE_IMAGES;
        if(grid_xp > grid_yp) {
            grid_max1 = grid_xp;
            if(grid_yp > grid_zp) {
                grid_max2 = grid_yp;
            }
            else {
                grid_max2 = grid_zp;
            }
         }
         else {
            grid_max1 = grid_yp;
              if(grid_xp > grid_zp) {
                  grid_max2 = grid_xp;
              }
              else {
                  grid_max2 = grid_zp;
              }
         }
         my_malloc(swbuf1x, 6 * MAX_TRADE_IMAGES * grid_max1 * grid_max2 * ct.THREADS_PER_NODE, REAL);
         my_malloc(swbuf2x, 6 * MAX_TRADE_IMAGES * grid_max1 * grid_max2 * ct.THREADS_PER_NODE, REAL);
         max_alloc = 6 * MAX_TRADE_IMAGES * grid_max1 * grid_max2 * ct.THREADS_PER_NODE;
         return;
    }
#if ASYNC_TRADES
    if(type == CENTRAL_FD) {
        trade_imagesx_central_async (f, w, dimx, dimy, dimz, images);
        return;
    }
    else {
        trade_imagesx_async (f, w, dimx, dimy, dimz, images);
        return;
    }
#endif
    //if(ct.trade_compression_level != 0) {
        //trade_imagesx_compress4 (f, w, dimx, dimy, dimz, images);
        //return;
    //}

#if MD_TIMERS
    REAL time1, time2, time3;
    time1 = my_crtc ();
#endif

#if HYBRID_MODEL
    basetag = get_thread_basetag();
    tid = get_thread_tid();
    if(tid < 0) tid = 0;
    if(is_loop_over_states()) ACTIVE_THREADS = ct.THREADS_PER_NODE;
#endif

    tim = 2 * images;

    incx = (dimy + tim) * (dimz + tim);
    incy = dimz + tim;
    incx0 = dimy * dimz;
    incy0 = dimz;

    zlen = dimx * dimy * images;
    ylen = dimx * images * (dimz + tim);
    xlen = images * (dimy + tim) * (dimz + tim);

    alloc = 2 * (xlen + ylen + zlen) * ct.THREADS_PER_NODE;
    // Verify that the required memory will fit into the statically allocated storage
    if(alloc > max_alloc)
        error_handler("Not enough memory. This should never happen.");


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

            QMD_dcopy (dimz, &f[iys], ione, &w[iys2 + images], ione);

        }                       /* end for */

    }                           /* end for */


    /* Collect the positive z-plane and negative z-planes */
    frdz1 = &swbuf1x[tid * zlen];
    frdz2 = &swbuf1x[tid * zlen + ACTIVE_THREADS * zlen];
    frdz2n = &swbuf2x[tid * zlen];
    frdz1n = &swbuf2x[tid * zlen + ACTIVE_THREADS * zlen];
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


#if MD_TIMERS
    time3 = my_crtc ();
#endif

    stop = zlen * ACTIVE_THREADS;
    scf_barrier_wait();
    if(tid == 0) {
        MPI_Sendrecv (&swbuf1x[0], stop, MPI_DOUBLE, nb_ids[NB_D], (1>>16),
                  &swbuf2x[0], stop, MPI_DOUBLE, nb_ids[NB_U], (1>>16), pct.grid_comm, &mrstatus);

        MPI_Sendrecv (&swbuf1x[stop], stop, MPI_DOUBLE, nb_ids[NB_U], (2>>16),
                  &swbuf2x[stop], stop, MPI_DOUBLE, nb_ids[NB_D], (2>>16), pct.grid_comm, &mrstatus);
    }
    scf_barrier_wait();

#if MD_TIMERS
    rmg_timings (TRADE_MPI_TIME, (my_crtc () - time3));
#endif

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
    frdy1 = &swbuf1x[tid * ylen];
    frdy2 = &swbuf1x[tid * ylen + ACTIVE_THREADS * ylen];
    frdy2n = &swbuf2x[tid * ylen];
    frdy1n = &swbuf2x[tid * ylen + ACTIVE_THREADS * ylen];
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


#if MD_TIMERS
    time3 = my_crtc ();
#endif
    stop = ylen * ACTIVE_THREADS;
    scf_barrier_wait();
    if(tid == 0) {
        MPI_Sendrecv (&swbuf1x[0], stop, MPI_DOUBLE, nb_ids[NB_S], (3>>16),
                  &swbuf2x[0], stop, MPI_DOUBLE, nb_ids[NB_N], (3>>16), pct.grid_comm, &mrstatus);

        MPI_Sendrecv (&swbuf1x[stop], stop, MPI_DOUBLE, nb_ids[NB_N], (4>>16),
                  &swbuf2x[stop], stop, MPI_DOUBLE, nb_ids[NB_S], (4>>16), pct.grid_comm, &mrstatus);
    }
    scf_barrier_wait();

#if MD_TIMERS
    rmg_timings (TRADE_MPI_TIME, (my_crtc () - time3));
#endif

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
    frdx1 = &swbuf1x[tid * xlen];
    frdx2 = &swbuf1x[tid * xlen + ACTIVE_THREADS * xlen];
    frdx2n = &swbuf2x[tid * xlen];
    frdx1n = &swbuf2x[tid * xlen + ACTIVE_THREADS * xlen];
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


#if MD_TIMERS
    time3 = my_crtc ();
#endif
    stop = xlen * ACTIVE_THREADS;
    scf_barrier_wait();
    if(tid == 0) {

        MPI_Sendrecv (&swbuf1x[0], stop, MPI_DOUBLE, nb_ids[NB_W], (5>>16),
                   &swbuf2x[0], stop, MPI_DOUBLE, nb_ids[NB_E], (5>>16), pct.grid_comm, &mrstatus);

        MPI_Sendrecv (&swbuf1x[stop], stop, MPI_DOUBLE, nb_ids[NB_E], (6>>16),
                  &swbuf2x[stop], stop, MPI_DOUBLE, nb_ids[NB_W], (6>>16), pct.grid_comm, &mrstatus);

    }
    scf_barrier_wait();

#if MD_TIMERS
    rmg_timings (TRADE_MPI_TIME, (my_crtc () - time3));
#endif


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


#if MD_TIMERS
    time2 = my_crtc ();
    rmg_timings (IMAGE_TIME, (time2 - time1));
#endif

}                               /* end trade_imagesx */


// Compresss/decompress image data inline. Uses a packed 5 byte format
// 
void trade_imagesx_compress5 (REAL * f, REAL * w, int dimx, int dimy, int dimz, int images)
{
    int ix, iy, iz, incx, incy, incx0, incy0, index, tim, ione = 1, retval;
    int ixs, iys, ixs2, iys2, c1, c2, alloc;
    int xlen, ylen, zlen, stop, tid, xlen_c, ylen_c, zlen_c;
    int *nb_ids, basetag=0;
    MPI_Status mrstatus;
    unsigned char *swbuf1x_c, *swbuf2x_c;
    unsigned char *frdx1, *frdx2, *frdy1, *frdy2, *frdz1, *frdz2;
    unsigned char *frdx1n, *frdx2n, *frdy1n, *frdy2n, *frdz1n, *frdz2n;
    int ACTIVE_THREADS = 1;
    union ieee754_double x;
    uint32_t x_32, *x_32p;
    uint8_t x_8, *xp1, *xp2;

#if MD_TIMERS
    REAL time1, time2, time3;
    time1 = my_crtc ();
#endif

#if HYBRID_MODEL
    basetag = get_thread_basetag();
    tid = get_thread_tid();
    if(tid < 0) tid = 0;
    if(is_loop_over_states()) ACTIVE_THREADS = ct.THREADS_PER_NODE;
#endif

    tim = 2 * images;

    incx = (dimy + tim) * (dimz + tim);
    incy = dimz + tim;
    incx0 = dimy * dimz;
    incy0 = dimz;

    zlen = dimx * dimy * images;
    ylen = dimx * images * (dimz + tim);
    xlen = images * (dimy + tim) * (dimz + tim);

    // Base byte pointers
    swbuf1x_c = (unsigned char *)swbuf1x;
    swbuf2x_c = (unsigned char *)swbuf2x;

    // Size in bytes required for each thread
    zlen_c = 5 * zlen;
    iz = 8 - zlen_c % 8; 
    zlen_c += iz;

    ylen_c = 5 * ylen;
    iy = 8 - ylen_c % 8; 
    ylen_c += iy;

    xlen_c = 5 * xlen;
    ix = 8 - xlen_c % 8;
    xlen_c += ix;


    alloc = 2 * (xlen + ylen + zlen) * ct.THREADS_PER_NODE;
    // Verify that the required memory will fit into the statically allocated storage
    if(alloc > max_alloc)
        error_handler("Not enough memory. This should never happen.");


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

            QMD_dcopy (dimz, &f[iys], ione, &w[iys2 + images], ione);

        }                       /* end for */

    }                           /* end for */


    /* Collect the positive z-plane and negative z-planes */
    frdz1 = &swbuf1x_c[tid * zlen_c];
    frdz2 = &swbuf1x_c[tid * zlen_c + ACTIVE_THREADS * zlen_c];
    xp1 = frdz1;
    xp2 = frdz2;
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

                x.d = w[index + images];
                x_32 = (x.ieee.negative << 31) + (x.ieee.exponent << 20) + x.ieee.mantissa0;
                x_8 = (unsigned char)(x.ieee.mantissa1 >> 24);
                x_32p = (uint32_t *)xp1;
                *x_32p = x_32;
                xp1[4] = x_8;
                xp1 += 5;
//                frdz1[iys + iz] = w[index + images];
                
                x.d = w[index + dimz];
                x_32 = (x.ieee.negative << 31) + (x.ieee.exponent << 20) + x.ieee.mantissa0;
                x_8 = (unsigned char)(x.ieee.mantissa1 >> 24);
                x_32p = (uint32_t *)xp2;
                *x_32p = x_32;
                xp2[4] = x_8;
                xp2 += 5;
//                frdz2[iys + iz] = w[index + dimz];


            }                   /* end for */

        }                       /* end for */

    }                           /* end for */


#if MD_TIMERS
    time3 = my_crtc ();
#endif

    stop = zlen_c * ACTIVE_THREADS;
    scf_barrier_wait();
    if(tid == 0) {
        MPI_Sendrecv (&swbuf1x_c[0], stop, MPI_BYTE, nb_ids[NB_D], (1>>16),
                  &swbuf2x_c[0], stop, MPI_BYTE, nb_ids[NB_U], (1>>16), pct.grid_comm, &mrstatus);

        MPI_Sendrecv (&swbuf1x_c[stop], stop, MPI_BYTE, nb_ids[NB_U], (2>>16),
                  &swbuf2x_c[stop], stop, MPI_BYTE, nb_ids[NB_D], (2>>16), pct.grid_comm, &mrstatus);
    }
    scf_barrier_wait();

#if MD_TIMERS
    rmg_timings (TRADE_MPI_TIME, (my_crtc () - time3));
#endif

    /* Unpack them */
    frdz1n = &swbuf2x_c[tid * zlen_c + ACTIVE_THREADS * zlen_c];
    frdz2n = &swbuf2x_c[tid * zlen_c];
    xp1 = frdz1n;
    xp2 = frdz2n;
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

                x_32p = (uint32_t *)xp1;
                x_32 = *x_32p;
                x_8 = xp1[4];
                x.ieee.negative = x_32 >> 31;
                x.ieee.exponent = (x_32 >> 20) & 0b11111111111;
                x.ieee.mantissa0 = x_32 & 0xfffff;
                x.ieee.mantissa1 = (x_8 << 24);
                w[index] = x.d;
                xp1 += 5;
//                w[index] = frdz1n[iys + iz];

                x_32p = (uint32_t *)xp2;
                x_32 = *x_32p;
                x_8 = xp2[4];
                x.ieee.negative = x_32 >> 31;
                x.ieee.exponent = (x_32 >> 20) & 0b11111111111;
                x.ieee.mantissa0 = x_32 & 0xfffff;
                x.ieee.mantissa1 = (x_8 << 24);
                w[index + c1] = x.d;
                xp2 += 5;
//                w[index + c1] = frdz2n[iys + iz];

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */

    /* Collect the north and south planes */
    frdy1 = &swbuf1x_c[tid * ylen_c];
    frdy2 = &swbuf1x_c[tid * ylen_c + ACTIVE_THREADS * ylen_c];
    xp1 = frdy1;
    xp2 = frdy2;
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

                x.d = w[index + c1];
                x_32 = (x.ieee.negative << 31) + (x.ieee.exponent << 20) + x.ieee.mantissa0;
                x_8 = (unsigned char)(x.ieee.mantissa1 >> 24);
                x_32p = (uint32_t *)xp1;
                *x_32p = x_32;
                xp1[4] = x_8;
                xp1 += 5;
//                frdy1[iys + iz] = w[index + c1];

                x.d = w[index + c2];
                x_32 = (x.ieee.negative << 31) + (x.ieee.exponent << 20) + x.ieee.mantissa0;
                x_8 = (unsigned char)(x.ieee.mantissa1 >> 24);
                x_32p = (uint32_t *)xp2;
                *x_32p = x_32;
                xp2[4] = x_8;
                xp2 += 5;
//                frdy2[iys + iz] = w[index + c2];

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */


#if MD_TIMERS
    time3 = my_crtc ();
#endif
    stop = ylen_c * ACTIVE_THREADS;
    scf_barrier_wait();
    if(tid == 0) {
        MPI_Sendrecv (&swbuf1x_c[0], stop, MPI_BYTE, nb_ids[NB_S], (3>>16),
                  &swbuf2x_c[0], stop, MPI_BYTE, nb_ids[NB_N], (3>>16), pct.grid_comm, &mrstatus);

        MPI_Sendrecv (&swbuf1x_c[stop], stop, MPI_BYTE, nb_ids[NB_N], (4>>16),
                  &swbuf2x_c[stop], stop, MPI_BYTE, nb_ids[NB_S], (4>>16), pct.grid_comm, &mrstatus);
    }
    scf_barrier_wait();

#if MD_TIMERS
    rmg_timings (TRADE_MPI_TIME, (my_crtc () - time3));
#endif

    /* Unpack them */
    frdy1n = &swbuf2x_c[tid * ylen_c + ACTIVE_THREADS * ylen_c];
    frdy2n = &swbuf2x_c[tid * ylen_c];
    xp1 = frdy1n;
    xp2 = frdy2n;
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

                x_32p = (uint32_t *)xp1;
                x_32 = *x_32p;
                x_8 = xp1[4];
                x.ieee.negative = x_32 >> 31;
                x.ieee.exponent = (x_32 >> 20) & 0b11111111111;
                x.ieee.mantissa0 = x_32 & 0xfffff;
                x.ieee.mantissa1 = (x_8 << 24);
                w[index] = x.d;
                xp1 += 5;
//                w[index] = frdy1n[iys + iz];

                x_32p = (uint32_t *)xp2;
                x_32 = *x_32p;
                x_8 = xp2[4];
                x.ieee.negative = x_32 >> 31;
                x.ieee.exponent = (x_32 >> 20) & 0b11111111111;
                x.ieee.mantissa0 = x_32 & 0xfffff;
                x.ieee.mantissa1 = (x_8 << 24);
                w[index + c1] = x.d;
                xp2 += 5;
//                w[index + c1] = frdy2n[iys + iz];

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */


    /* Collect the east and west planes */
    frdx1 = &swbuf1x_c[tid * xlen_c];
    frdx2 = &swbuf1x_c[tid * xlen_c + ACTIVE_THREADS * xlen_c];
    c1 = images * incx;
    c2 = dimx * incx;
    xp1 = frdx1;
    xp2 = frdx2;
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

                x.d = w[index + c1];
                x_32 = (x.ieee.negative << 31) + (x.ieee.exponent << 20) + x.ieee.mantissa0;
                x_8 = (unsigned char)(x.ieee.mantissa1 >> 24);
                x_32p = (uint32_t *)xp1;
                *x_32p = x_32;
                xp1[4] = x_8;
                xp1 += 5;
//                frdx1[iys + iz] = w[index + c1];

                x.d = w[index + c2];
                x_32 = (x.ieee.negative << 31) + (x.ieee.exponent << 20) + x.ieee.mantissa0;
                x_8 = (unsigned char)(x.ieee.mantissa1 >> 24);
                x_32p = (uint32_t *)xp2;
                *x_32p = x_32;
                xp2[4] = x_8;
                xp2 += 5;
//                frdx2[iys + iz] = w[index + c2];

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */


#if MD_TIMERS
    time3 = my_crtc ();
#endif
    stop = xlen_c * ACTIVE_THREADS;
    scf_barrier_wait();
    if(tid == 0) {

        MPI_Sendrecv (&swbuf1x_c[0], stop, MPI_BYTE, nb_ids[NB_W], (5>>16),
                   &swbuf2x_c[0], stop, MPI_BYTE, nb_ids[NB_E], (5>>16), pct.grid_comm, &mrstatus);

        MPI_Sendrecv (&swbuf1x_c[stop], stop, MPI_BYTE, nb_ids[NB_E], (6>>16),
                  &swbuf2x_c[stop], stop, MPI_BYTE, nb_ids[NB_W], (6>>16), pct.grid_comm, &mrstatus);

    }
    scf_barrier_wait();

#if MD_TIMERS
    rmg_timings (TRADE_MPI_TIME, (my_crtc () - time3));
#endif


    /* Unpack them */
    frdx1n = &swbuf2x_c[tid * xlen_c + ACTIVE_THREADS * xlen_c];
    frdx2n = &swbuf2x_c[tid * xlen_c];
    xp1 = frdx1n;
    xp2 = frdx2n;
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

                x_32p = (uint32_t *)xp1;
                x_32 = *x_32p;
                x_8 = xp1[4];
                x.ieee.negative = x_32 >> 31;
                x.ieee.exponent = (x_32 >> 20) & 0b11111111111;
                x.ieee.mantissa0 = x_32 & 0xfffff;
                x.ieee.mantissa1 = (x_8 << 24);
                w[index] = x.d;
                xp1 += 5;
//                w[index] = frdx1n[iys + iz];

                x_32p = (uint32_t *)xp2;
                x_32 = *x_32p;
                x_8 = xp2[4];
                x.ieee.negative = x_32 >> 31;
                x.ieee.exponent = (x_32 >> 20) & 0b11111111111;
                x.ieee.mantissa0 = x_32 & 0xfffff;
                x.ieee.mantissa1 = (x_8 << 24);
                w[index + c1] = x.d;
                xp2 += 5;
//                w[index + c1] = frdx2n[iys + iz];

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */


#if MD_TIMERS
    time2 = my_crtc ();
    rmg_timings (IMAGE_TIME, (time2 - time1));
#endif

}                               /* end trade_imagesx_compress */

// Compresss/decompress image data inline. Works by converting double to float
// 
void trade_imagesx_compress4 (REAL * f, REAL * w, int dimx, int dimy, int dimz, int images)
{
    int ix, iy, iz, incx, incy, incx0, incy0, index, tim, ione = 1, retval;
    int ixs, iys, ixs2, iys2, c1, c2, alloc;
    int xlen, ylen, zlen, stop, tid, xlen_c, ylen_c, zlen_c;
    int *nb_ids, basetag=0;
    MPI_Status mrstatus;
    unsigned char *frdx1, *frdx2, *frdy1, *frdy2, *frdz1, *frdz2;
    unsigned char *frdx1n, *frdx2n, *frdy1n, *frdy2n, *frdz1n, *frdz2n;
    unsigned char *swbuf1x_c, *swbuf2x_c;
    int ACTIVE_THREADS = 1;
    float *fxp1, *fxp2;

#if MD_TIMERS
    REAL time1, time2, time3;
    time1 = my_crtc ();
#endif

#if HYBRID_MODEL
    basetag = get_thread_basetag();
    tid = get_thread_tid();
    if(tid < 0) tid = 0;
    if(is_loop_over_states()) ACTIVE_THREADS = ct.THREADS_PER_NODE;
#endif

    tim = 2 * images;

    incx = (dimy + tim) * (dimz + tim);
    incy = dimz + tim;
    incx0 = dimy * dimz;
    incy0 = dimz;

    zlen = dimx * dimy * images;
    ylen = dimx * images * (dimz + tim);
    xlen = images * (dimy + tim) * (dimz + tim);

    // Base byte pointers
    swbuf1x_c = (unsigned char *)swbuf1x;
    swbuf2x_c = (unsigned char *)swbuf2x;

    // Size in bytes required for each thread
    zlen_c = sizeof(float) * zlen;
    ylen_c = sizeof(float) * ylen;
    xlen_c = sizeof(float) * xlen;

    alloc = 2 * (xlen + ylen + zlen) * ct.THREADS_PER_NODE;
    // Verify that the required memory will fit into the statically allocated storage
    if(alloc > max_alloc)
        error_handler("Not enough memory. This should never happen.");


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

            QMD_dcopy (dimz, &f[iys], ione, &w[iys2 + images], ione);

        }                       /* end for */

    }                           /* end for */


    /* Collect the positive z-plane and negative z-planes */
    fxp1 = (float *)&swbuf1x_c[tid * zlen_c];
    fxp2 = (float *)&swbuf1x_c[tid * zlen_c + ACTIVE_THREADS * zlen_c];
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

                *fxp1 = (float)w[index + images];
                fxp1++;
//              frdz1[iys + iz] = w[index + images];

                *fxp2 = (float)w[index + dimz];
                fxp2++;
//              frdz2[iys + iz] = w[index + dimz];

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */


#if MD_TIMERS
    time3 = my_crtc ();
#endif

    stop = zlen_c * ACTIVE_THREADS;
    scf_barrier_wait();
    if(tid == 0) {
        MPI_Sendrecv (&swbuf1x_c[0], stop, MPI_BYTE, nb_ids[NB_D], (1>>16),
                  &swbuf2x_c[0], stop, MPI_BYTE, nb_ids[NB_U], (1>>16), pct.grid_comm, &mrstatus);

        MPI_Sendrecv (&swbuf1x_c[stop], stop, MPI_BYTE, nb_ids[NB_U], (2>>16),
                  &swbuf2x_c[stop], stop, MPI_BYTE, nb_ids[NB_D], (2>>16), pct.grid_comm, &mrstatus);
    }
    scf_barrier_wait();

#if MD_TIMERS
    rmg_timings (TRADE_MPI_TIME, (my_crtc () - time3));
#endif

    /* Unpack them */
    fxp1 = (float *)&swbuf2x_c[tid * zlen_c + ACTIVE_THREADS * zlen_c];
    fxp2 = (float *)&swbuf2x_c[tid * zlen_c];
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

                w[index] = (double)*fxp1;
                fxp1++;
//              w[index] = frdz1n[iys + iz];

                w[index + c1] = (double)*fxp2;
                fxp2++;
//              w[index + c1] = frdz2n[iys + iz];

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */

    /* Collect the north and south planes */
    fxp1 = (float *)&swbuf1x_c[tid * ylen_c];
    fxp2 = (float *)&swbuf1x_c[tid * ylen_c + ACTIVE_THREADS * ylen_c];
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

                *fxp1 = (float)w[index + c1];
                fxp1++;
//              frdy1[iys + iz] = w[index + c1];

                *fxp2 = (float)w[index + c2];
                fxp2++;
//              frdy2[iys + iz] = w[index + c2];

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */


#if MD_TIMERS
    time3 = my_crtc ();
#endif
    stop = ylen_c * ACTIVE_THREADS;
    scf_barrier_wait();
    if(tid == 0) {
        MPI_Sendrecv (&swbuf1x_c[0], stop, MPI_BYTE, nb_ids[NB_S], (3>>16),
                  &swbuf2x_c[0], stop, MPI_BYTE, nb_ids[NB_N], (3>>16), pct.grid_comm, &mrstatus);

        MPI_Sendrecv (&swbuf1x_c[stop], stop, MPI_BYTE, nb_ids[NB_N], (4>>16),
                  &swbuf2x_c[stop], stop, MPI_BYTE, nb_ids[NB_S], (4>>16), pct.grid_comm, &mrstatus);
    }
    scf_barrier_wait();

#if MD_TIMERS
    rmg_timings (TRADE_MPI_TIME, (my_crtc () - time3));
#endif

    /* Unpack them */
    fxp1 = (float *)&swbuf2x_c[tid * ylen_c + ACTIVE_THREADS * ylen_c];
    fxp2 = (float *)&swbuf2x_c[tid * ylen_c];
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

                w[index] = (double)*fxp1;
                fxp1++;
//              w[index] = frdy1n[iys + iz];

                w[index + c1] = (double)*fxp2;
                fxp2++;
//              w[index + c1] = frdy2n[iys + iz];

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */


    /* Collect the east and west planes */
    fxp1 = (float *)&swbuf1x_c[tid * xlen_c];
    fxp2 = (float *)&swbuf1x_c[tid * xlen_c + ACTIVE_THREADS * xlen_c];
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

                *fxp1 = (float)w[index + c1];
                fxp1++;
//              frdx1[iys + iz] = w[index + c1];

                *fxp2 = (float)w[index + c2];
                fxp2++;
//              frdx2[iys + iz] = w[index + c2];

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */


#if MD_TIMERS
    time3 = my_crtc ();
#endif
    stop = xlen_c * ACTIVE_THREADS;
    scf_barrier_wait();
    if(tid == 0) {

        MPI_Sendrecv (&swbuf1x_c[0], stop, MPI_BYTE, nb_ids[NB_W], (5>>16),
                   &swbuf2x_c[0], stop, MPI_BYTE, nb_ids[NB_E], (5>>16), pct.grid_comm, &mrstatus);

        MPI_Sendrecv (&swbuf1x_c[stop], stop, MPI_BYTE, nb_ids[NB_E], (6>>16),
                  &swbuf2x_c[stop], stop, MPI_BYTE, nb_ids[NB_W], (6>>16), pct.grid_comm, &mrstatus);

    }
    scf_barrier_wait();

#if MD_TIMERS
    rmg_timings (TRADE_MPI_TIME, (my_crtc () - time3));
#endif


    /* Unpack them */
    fxp1 = (float *)&swbuf2x_c[tid * xlen_c + ACTIVE_THREADS * xlen_c];
    fxp2 = (float *)&swbuf2x_c[tid * xlen_c];
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

                w[index] = (double)*fxp1;
                fxp1++;
//              w[index] = frdx1n[iys + iz];

                w[index + c1] = (double)*fxp2;
                fxp2++;
//              w[index + c1] = frdx2n[iys + iz];

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */


#if MD_TIMERS
    time2 = my_crtc ();
    rmg_timings (IMAGE_TIME, (time2 - time1));
#endif

}                               /* end trade_imagesx_compress */


// Single precision version of trade_imagesx
void trade_imagesx_f (rmg_float_t *f, rmg_float_t *w, int dimx, int dimy, int dimz, int images, int type)
{
    int ix, iy, iz, incx, incy, incx0, incy0, index, tim, ione = 1, retval;
    int ixs, iys, ixs2, iys2, c1, c2, alloc;
    int xlen, ylen, zlen, stop, tid, grid_max1, grid_max2;
    int *nb_ids, basetag=0;
    MPI_Status mrstatus;
    rmg_float_t *frdx1, *frdx2, *frdy1, *frdy2, *frdz1, *frdz2;
    rmg_float_t *frdx1n, *frdx2n, *frdy1n, *frdy2n, *frdz1n, *frdz2n;
    rmg_float_t *swbuf1x_f, *swbuf2x_f;
    int ACTIVE_THREADS = 1;

#if MD_TIMERS
    REAL time1, time2, time3;
    time1 = my_crtc ();
#endif

#if HYBRID_MODEL
    basetag = get_thread_basetag();
    tid = get_thread_tid();
    if(tid < 0) tid = 0;
    if(is_loop_over_states()) ACTIVE_THREADS = ct.THREADS_PER_NODE;
#endif

    swbuf1x_f = (rmg_float_t *)swbuf1x;
    swbuf2x_f = (rmg_float_t *)swbuf2x;

    tim = 2 * images;

    incx = (dimy + tim) * (dimz + tim);
    incy = dimz + tim;
    incx0 = dimy * dimz;
    incy0 = dimz;

    zlen = dimx * dimy * images;
    ylen = dimx * images * (dimz + tim);
    xlen = images * (dimy + tim) * (dimz + tim);

    alloc = 2 * (xlen + ylen + zlen) * ct.THREADS_PER_NODE;
    // Verify that the required memory will fit into the statically allocated storage
    if(alloc > max_alloc)
        error_handler("Not enough memory. This should never happen.");


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

            QMD_scopy (dimz, &f[iys], ione, &w[iys2 + images], ione);

        }                       /* end for */

    }                           /* end for */


    /* Collect the positive z-plane and negative z-planes */
    frdz1 = (rmg_float_t *)&swbuf1x_f[tid * zlen];
    frdz2 = (rmg_float_t *)&swbuf1x_f[tid * zlen + ACTIVE_THREADS * zlen];
    frdz2n = (rmg_float_t *)&swbuf2x_f[tid * zlen];
    frdz1n = (rmg_float_t *)&swbuf2x_f[tid * zlen + ACTIVE_THREADS * zlen];
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


#if MD_TIMERS
    time3 = my_crtc ();
#endif

    stop = zlen * ACTIVE_THREADS;
    scf_barrier_wait();
    if(tid == 0) {
        MPI_Sendrecv (&swbuf1x_f[0], stop, MPI_FLOAT, nb_ids[NB_D], (1>>16),
                  &swbuf2x_f[0], stop, MPI_FLOAT, nb_ids[NB_U], (1>>16), pct.grid_comm, &mrstatus);

        MPI_Sendrecv (&swbuf1x_f[stop], stop, MPI_FLOAT, nb_ids[NB_U], (2>>16),
                  &swbuf2x_f[stop], stop, MPI_FLOAT, nb_ids[NB_D], (2>>16), pct.grid_comm, &mrstatus);
    }
    scf_barrier_wait();

#if MD_TIMERS
    rmg_timings (TRADE_MPI_TIME, (my_crtc () - time3));
#endif

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
    frdy1 = (rmg_float_t *)&swbuf1x_f[tid * ylen];
    frdy2 = (rmg_float_t *)&swbuf1x_f[tid * ylen + ACTIVE_THREADS * ylen];
    frdy2n = (rmg_float_t *)&swbuf2x_f[tid * ylen];
    frdy1n = (rmg_float_t *)&swbuf2x_f[tid * ylen + ACTIVE_THREADS * ylen];
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


#if MD_TIMERS
    time3 = my_crtc ();
#endif
    stop = ylen * ACTIVE_THREADS;
    scf_barrier_wait();
    if(tid == 0) {
        MPI_Sendrecv (&swbuf1x_f[0], stop, MPI_FLOAT, nb_ids[NB_S], (3>>16),
                  &swbuf2x_f[0], stop, MPI_FLOAT, nb_ids[NB_N], (3>>16), pct.grid_comm, &mrstatus);

        MPI_Sendrecv (&swbuf1x_f[stop], stop, MPI_FLOAT, nb_ids[NB_N], (4>>16),
                  &swbuf2x_f[stop], stop, MPI_FLOAT, nb_ids[NB_S], (4>>16), pct.grid_comm, &mrstatus);
    }
    scf_barrier_wait();

#if MD_TIMERS
    rmg_timings (TRADE_MPI_TIME, (my_crtc () - time3));
#endif

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
    frdx1 = (rmg_float_t *)&swbuf1x_f[tid * xlen];
    frdx2 = (rmg_float_t *)&swbuf1x_f[tid * xlen + ACTIVE_THREADS * xlen];
    frdx2n = (rmg_float_t *)&swbuf2x_f[tid * xlen];
    frdx1n = (rmg_float_t *)&swbuf2x_f[tid * xlen + ACTIVE_THREADS * xlen];
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


#if MD_TIMERS
    time3 = my_crtc ();
#endif
    stop = xlen * ACTIVE_THREADS;
    scf_barrier_wait();
    if(tid == 0) {

        MPI_Sendrecv (&swbuf1x_f[0], stop, MPI_FLOAT, nb_ids[NB_W], (5>>16),
                   &swbuf2x_f[0], stop, MPI_FLOAT, nb_ids[NB_E], (5>>16), pct.grid_comm, &mrstatus);

        MPI_Sendrecv (&swbuf1x_f[stop], stop, MPI_FLOAT, nb_ids[NB_E], (6>>16),
                  &swbuf2x_f[stop], stop, MPI_FLOAT, nb_ids[NB_W], (6>>16), pct.grid_comm, &mrstatus);

    }
    scf_barrier_wait();

#if MD_TIMERS
    rmg_timings (TRADE_MPI_TIME, (my_crtc () - time3));
#endif


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


#if MD_TIMERS
    time2 = my_crtc ();
    rmg_timings (IMAGE_TIME, (time2 - time1));
#endif

}                               /* end trade_imagesx */
