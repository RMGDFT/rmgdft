/************************** SVN Revision Information **************************
 **    $Id: trade_imagesx.c 1694 2012-04-09 19:46:45Z ebriggs $    **
******************************************************************************/

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "grid.h"
#include "main.h"
#include "common_prototypes.h"



int TRADE_GRID_EDGES;
int GRID_MAX1;
int GRID_MAX2;


#include "hybrid.h"
#include <pthread.h>


/*
 * INPUTS
 * f[dimx*dimy*dimz] - raw array without images. pack_ptos should not be called, this is 
 *                     handled inside the function now
 * OUTPUT
 * w[(dimx+2*images) * (dimy+2*images) * (dimz+2*images)]
 *    This array is completely filled, i.e. the original data is filled in and then 
 *    the image data are added*/


#define MAX_IMG2 (MAX_TRADE_IMAGES*MAX_TRADE_IMAGES)
#define MAX_IMG3 (MAX_TRADE_IMAGES*MAX_TRADE_IMAGES*MAX_TRADE_IMAGES)


#if ASYNC_TRADES

// Rank of target node based on offsets from current node
static int target_node[3][3][3];

static rmg_double_t *frdx1, *frdx2, *frdy1, *frdy2, *frdz1, *frdz2;
static rmg_double_t *frdx1n, *frdx2n, *frdy1n, *frdy2n, *frdz1n, *frdz2n;
static rmg_double_t *yzpsms_r, *yzpsps_r, *yzmsms_r, *yzmsps_r;
static rmg_double_t *yzpsms_s, *yzpsps_s, *yzmsms_s, *yzmsps_s;
static rmg_double_t *xzpsms_r, *xzpsps_r, *xzmsms_r, *xzmsps_r;
static rmg_double_t *xzpsms_s, *xzpsps_s, *xzmsms_s, *xzmsps_s;
static rmg_double_t *xypsms_r, *xypsps_r, *xymsms_r, *xymsps_r;
static rmg_double_t *xypsms_s, *xypsps_s, *xymsms_s, *xymsps_s;
static rmg_double_t *m0_s, *m0_r;

static rmg_float_t *frdx1_f, *frdx2_f, *frdy1_f, *frdy2_f, *frdz1_f, *frdz2_f;
static rmg_float_t *frdx1n_f, *frdx2n_f, *frdy1n_f, *frdy2n_f, *frdz1n_f, *frdz2n_f;
static rmg_float_t *yzpsms_r_f, *yzpsps_r_f, *yzmsms_r_f, *yzmsps_r_f;
static rmg_float_t *yzpsms_s_f, *yzpsps_s_f, *yzmsms_s_f, *yzmsps_s_f;
static rmg_float_t *xzpsms_r_f, *xzpsps_r_f, *xzmsms_r_f, *xzmsps_r_f;
static rmg_float_t *xzpsms_s_f, *xzpsps_s_f, *xzmsms_s_f, *xzmsps_s_f;
static rmg_float_t *xypsms_r_f, *xypsps_r_f, *xymsms_r_f, *xymsps_r_f;
static rmg_float_t *xypsms_s_f, *xypsps_s_f, *xymsms_s_f, *xymsps_s_f;
static rmg_float_t *m0_s_f, *m0_r_f;

static MPI_Request sreqs[26];
static MPI_Request rreqs[26];


void trade_imagesx_async (rmg_double_t * f, rmg_double_t * w, int dimx, int dimy, int dimz, int images)
{
    int ix, iy, iz, ix1, iy1, iz1, incx, incy, incx0, incy0, index, tim, ione = 1;
    int ixs, iys, ixs2, iys2, c1, c2, c3, idx, idx1, img3;
    int xlen, ylen, zlen, stop, yzlen, xylen, xzlen;
    int tid=0, corner_node_stride, node_idx, retval;
    int ACTIVE_THREADS = 1;

#if MD_TIMERS
    rmg_double_t time1, time2;
    time1 = my_crtc ();
#endif

    if(images > MAX_TRADE_IMAGES) {
        error_handler ("Images count too high in trade_imagesx_async. Modify and recompile may be required.\n");
    }

    tid = get_thread_tid();
    if(tid < 0) tid = 0;
    if(is_loop_over_states()) ACTIVE_THREADS = ct.THREADS_PER_NODE;

    corner_node_stride = ct.THREADS_PER_NODE * MAX_IMG3; 
    tim = 2 * images;
    img3 = images*images*images;

    incx = (dimy + tim) * (dimz + tim);
    incy = dimz + tim;
    incx0 = dimy * dimz;
    incy0 = dimz;

    zlen = dimx * dimy * images;
    ylen = dimx * dimz * images;
    xlen = dimy * dimz * images;

    yzlen = images * images * dimx;
    xylen = images * images * dimz;
    xzlen = images * images * dimy;

    // Thread 0 posts all of the receives
    if(tid == 0) {

        // Corner nodes
        RMG_MPI_trade( &m0_r[0*corner_node_stride], ACTIVE_THREADS * img3, RMG_MPI_IRECV, -1, -1, -1, pct.grid_comm, 18, &rreqs[18]);
        RMG_MPI_trade( &m0_r[1*corner_node_stride], ACTIVE_THREADS * img3, RMG_MPI_IRECV, -1, -1, 1, pct.grid_comm, 19, &rreqs[19]);
        RMG_MPI_trade( &m0_r[2*corner_node_stride], ACTIVE_THREADS * img3, RMG_MPI_IRECV, -1, 1, -1, pct.grid_comm, 20, &rreqs[20]);
        RMG_MPI_trade( &m0_r[3*corner_node_stride], ACTIVE_THREADS * img3, RMG_MPI_IRECV, -1, 1, 1, pct.grid_comm, 21, &rreqs[21]);
        RMG_MPI_trade( &m0_r[4*corner_node_stride], ACTIVE_THREADS * img3, RMG_MPI_IRECV, 1, -1, -1, pct.grid_comm, 22, &rreqs[22]);
        RMG_MPI_trade( &m0_r[5*corner_node_stride], ACTIVE_THREADS * img3, RMG_MPI_IRECV, 1, -1, 1, pct.grid_comm, 23, &rreqs[23]);
        RMG_MPI_trade( &m0_r[6*corner_node_stride], ACTIVE_THREADS * img3, RMG_MPI_IRECV, 1, 1, -1, pct.grid_comm, 24, &rreqs[24]);
        RMG_MPI_trade( &m0_r[7*corner_node_stride], ACTIVE_THREADS * img3, RMG_MPI_IRECV, 1, 1, 1, pct.grid_comm, 25, &rreqs[25]);

        // The z planes
        RMG_MPI_trade(frdz2n, ACTIVE_THREADS * zlen, RMG_MPI_IRECV, 0, 0, 1, pct.grid_comm, 0, &rreqs[0]);
        RMG_MPI_trade(frdz1n, ACTIVE_THREADS * zlen, RMG_MPI_IRECV, 0, 0, -1, pct.grid_comm, 1, &rreqs[1]);

        // The y planes
        RMG_MPI_trade(frdy2n, ACTIVE_THREADS * ylen, RMG_MPI_IRECV, 0, 1, 0, pct.grid_comm, 2, &rreqs[2]);
        RMG_MPI_trade(frdy1n, ACTIVE_THREADS * ylen, RMG_MPI_IRECV, 0, -1, 0, pct.grid_comm, 3, &rreqs[3]);

        // The x planes
        RMG_MPI_trade(frdx2n, ACTIVE_THREADS * xlen, RMG_MPI_IRECV, 1, 0, 0, pct.grid_comm, 4, &rreqs[4]);
        RMG_MPI_trade(frdx1n, ACTIVE_THREADS * xlen, RMG_MPI_IRECV, -1, 0, 0, pct.grid_comm, 5, &rreqs[5]);

        // And the yz-plane edges
        RMG_MPI_trade(yzpsms_r, ACTIVE_THREADS * yzlen, RMG_MPI_IRECV, 0, 1, -1, pct.grid_comm, 6, &rreqs[6]);
        RMG_MPI_trade(yzpsps_r, ACTIVE_THREADS * yzlen, RMG_MPI_IRECV, 0, 1, 1, pct.grid_comm, 7, &rreqs[7]);
        RMG_MPI_trade(yzmsms_r, ACTIVE_THREADS * yzlen, RMG_MPI_IRECV, 0, -1, -1, pct.grid_comm, 8, &rreqs[8]);
        RMG_MPI_trade(yzmsps_r, ACTIVE_THREADS * yzlen, RMG_MPI_IRECV, 0, -1, 1, pct.grid_comm, 9, &rreqs[9]);

        // And the xy-plane edges
        RMG_MPI_trade(xypsms_r, ACTIVE_THREADS * xylen, RMG_MPI_IRECV, 1, -1, 0, pct.grid_comm, 10, &rreqs[10]);
        RMG_MPI_trade(xypsps_r, ACTIVE_THREADS * xylen, RMG_MPI_IRECV, 1, 1, 0, pct.grid_comm, 11, &rreqs[11]);
        RMG_MPI_trade(xymsms_r, ACTIVE_THREADS * xylen, RMG_MPI_IRECV, -1, -1, 0, pct.grid_comm, 12, &rreqs[12]);
        RMG_MPI_trade(xymsps_r, ACTIVE_THREADS * xylen, RMG_MPI_IRECV, -1, 1, 0, pct.grid_comm, 13, &rreqs[13]);

        // And the xz-plane edges
        RMG_MPI_trade(xzpsms_r, ACTIVE_THREADS * xzlen, RMG_MPI_IRECV, 1, 0, -1, pct.grid_comm, 14, &rreqs[14]);
        RMG_MPI_trade(xzpsps_r, ACTIVE_THREADS * xzlen, RMG_MPI_IRECV, 1, 0, 1, pct.grid_comm, 15, &rreqs[15]);
        RMG_MPI_trade(xzmsms_r, ACTIVE_THREADS * xzlen, RMG_MPI_IRECV, -1, 0, -1, pct.grid_comm, 16, &rreqs[16]);
        RMG_MPI_trade(xzmsps_r, ACTIVE_THREADS * xzlen, RMG_MPI_IRECV, -1, 0, 1, pct.grid_comm, 17, &rreqs[17]);

    }

    // Collect data for the corner nodes
    node_idx = 0;
    for(ix = -1;ix <= 1;ix+=2) 
    {
        c1 = (ix == -1) ? (dimx-images)*incx0 : 0;

        for(iy = -1;iy <= 1;iy+=2) 
        {
            c2 = (iy == -1) ? (dimy-images)*incy0 : 0;

            for(iz = -1;iz <= 1;iz+=2) 
            {
                c3 = (iz == -1) ? (dimz-images) : 0;

                // Pack the send arrays
                idx1 = node_idx*corner_node_stride + tid * img3;
                for(ix1 = 0;ix1 < images;ix1++) 
                {
                    ixs2 = ix1*incx0;

                    for(iy1 = 0;iy1 < images;iy1++) 
                    {
                        iys2 = ixs2 + iy1*incy0;

                        for(iz1 = 0;iz1 < images;iz1++) 
                        {
                            index = iys2 + iz1;
                            m0_s[idx1] = f[index + c1 + c2 + c3];
                            idx1++;
                        }
                    }
                }

                node_idx++;
            }
        }
    }


    // Send them
    thread_barrier_wait();
    if(tid == 0) {

        // Corners
        RMG_MPI_trade( &m0_s[0*corner_node_stride], ACTIVE_THREADS * img3, RMG_MPI_ISEND, 1, 1, 1, pct.grid_comm, 18, &sreqs[18]);
        RMG_MPI_trade( &m0_s[1*corner_node_stride], ACTIVE_THREADS * img3, RMG_MPI_ISEND, 1, 1, -1, pct.grid_comm, 19, &sreqs[19]);
        RMG_MPI_trade( &m0_s[2*corner_node_stride], ACTIVE_THREADS * img3, RMG_MPI_ISEND, 1, -1, 1, pct.grid_comm, 20, &sreqs[20]);
        RMG_MPI_trade( &m0_s[3*corner_node_stride], ACTIVE_THREADS * img3, RMG_MPI_ISEND, 1, -1, -1, pct.grid_comm, 21, &sreqs[21]);
        RMG_MPI_trade( &m0_s[4*corner_node_stride], ACTIVE_THREADS * img3, RMG_MPI_ISEND, -1, 1, 1, pct.grid_comm, 22, &sreqs[22]);
        RMG_MPI_trade( &m0_s[5*corner_node_stride], ACTIVE_THREADS * img3, RMG_MPI_ISEND, -1, 1, -1, pct.grid_comm, 23, &sreqs[23]);
        RMG_MPI_trade( &m0_s[6*corner_node_stride], ACTIVE_THREADS * img3, RMG_MPI_ISEND, -1, -1, 1, pct.grid_comm, 24, &sreqs[24]);
        RMG_MPI_trade( &m0_s[7*corner_node_stride], ACTIVE_THREADS * img3, RMG_MPI_ISEND, -1, -1, -1, pct.grid_comm, 25, &sreqs[25]);

    }


    /* Collect the positive z-plane and negative z-planes */
    c1 = (dimz-images);
    idx = tid * dimx * dimy * images;

    if(images == 2) {

        for (ix = 0; ix < dimx; ix++)
        {
            ixs2 = ix * incx0;
            for (iy = 0; iy < dimy; iy++)
            {
                iys2 = ixs2 + iy * incy0;

                    frdz1[idx] = f[iys2];
                    frdz2[idx] = f[iys2+c1];
                    frdz1[idx+1] = f[iys2+1];
                    frdz2[idx+1] = f[iys2+c1+1];
                    idx+=2;

            }                       /* end for */

        }                           /* end for */

    }
    else{

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

    }


    /* Collect the north and south planes */
    c1 = (dimy - images)*incy0;
    idx = tid * dimx * dimz * images;
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
    idx = tid * dimy * dimz * images;
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


    thread_barrier_wait();
    if(tid == 0) {

        // Send z planes
        RMG_MPI_trade(frdz1, ACTIVE_THREADS * zlen, RMG_MPI_ISEND, 0, 0, -1, pct.grid_comm, 0, &sreqs[0]);
        RMG_MPI_trade(frdz2, ACTIVE_THREADS * zlen, RMG_MPI_ISEND, 0, 0, 1, pct.grid_comm, 1, &sreqs[1]);

        // Send the north and south planes
        RMG_MPI_trade(frdy1, ACTIVE_THREADS * ylen, RMG_MPI_ISEND, 0, -1, 0, pct.grid_comm, 2, &sreqs[2]);
        RMG_MPI_trade(frdy2, ACTIVE_THREADS * ylen, RMG_MPI_ISEND, 0, 1, 0, pct.grid_comm, 3, &sreqs[3]);

        // Send the east and west planes
        RMG_MPI_trade(frdx1, ACTIVE_THREADS * xlen, RMG_MPI_ISEND, -1, 0, 0, pct.grid_comm, 4, &sreqs[4]);
        RMG_MPI_trade(frdx2, ACTIVE_THREADS * xlen, RMG_MPI_ISEND, 1, 0, 0, pct.grid_comm, 5, &sreqs[5]);

    }


    // Pack yz-plane edges
    c1 = (dimy-images)*incy0;
    c2 = (dimz - images);
    idx = tid * dimx * images * images;
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
    c1 = (dimy-images)*incy0;
    c2 = (dimx-images)*incx0;
    idx = tid * dimz * images * images;
    for (ix = 0; ix < images; ix++)
    {
        ixs2 = ix * incx0;
        for(iy = 0;iy < images;iy++) {
            iys2 = ixs2 + iy*incy0;
            for(iz = 0;iz < dimz;iz++) {
                index = iys2 + iz;
                xypsms_s[idx] = f[index + c1]; 
                xypsps_s[idx] = f[index];
                xymsms_s[idx] = f[index + c2 + c1];
                xymsps_s[idx] = f[index + c2];
                idx++;
            }
        }
    }


    // Pack xz plane edges
    c1 = (dimz-images);
    c2 = (dimx-images)*incx0;
    idx = tid * dimy * images * images;
    for (ix = 0; ix < images; ix++)
    {
        ixs2 = ix * incx0;
        for(iy = 0;iy < dimy;iy++) {
            iys2 = ixs2 + iy*incy0;
            for(iz = 0;iz < images;iz++) {
                index = iys2 + iz;
                xzpsms_s[idx] = f[index + c1]; 
                xzpsps_s[idx] = f[index];
                xzmsms_s[idx] = f[index + c2 + c1];
                xzmsps_s[idx] = f[index + c2];
                idx++;
            }
        }
    }


    // Send them
    thread_barrier_wait();
    if(tid == 0) {

        // Send yz-plane edges
        RMG_MPI_trade(yzpsms_s, ACTIVE_THREADS * yzlen, RMG_MPI_ISEND, 0, -1, 1, pct.grid_comm, 6, &sreqs[6]);
        RMG_MPI_trade(yzpsps_s, ACTIVE_THREADS * yzlen, RMG_MPI_ISEND, 0, -1, -1, pct.grid_comm, 7, &sreqs[7]);
        RMG_MPI_trade(yzmsms_s, ACTIVE_THREADS * yzlen, RMG_MPI_ISEND, 0, 1, 1, pct.grid_comm, 8, &sreqs[8]);
        RMG_MPI_trade(yzmsps_s, ACTIVE_THREADS * yzlen, RMG_MPI_ISEND, 0, 1, -1, pct.grid_comm, 9, &sreqs[9]);

        // Send xy-plane edges
        RMG_MPI_trade(xypsms_s, ACTIVE_THREADS * xylen, RMG_MPI_ISEND, -1, 1, 0, pct.grid_comm, 10, &sreqs[10]);
        RMG_MPI_trade(xypsps_s, ACTIVE_THREADS * xylen, RMG_MPI_ISEND, -1, -1, 0, pct.grid_comm, 11, &sreqs[11]);
        RMG_MPI_trade(xymsms_s, ACTIVE_THREADS * xylen, RMG_MPI_ISEND, 1, 1, 0, pct.grid_comm, 12, &sreqs[12]);
        RMG_MPI_trade(xymsps_s, ACTIVE_THREADS * xylen, RMG_MPI_ISEND, 1, -1, 0, pct.grid_comm, 13, &sreqs[13]);

        // Send xz-plane edges
        RMG_MPI_trade(xzpsms_s, ACTIVE_THREADS * xzlen, RMG_MPI_ISEND, -1, 0, 1, pct.grid_comm, 14, &sreqs[14]);
        RMG_MPI_trade(xzpsps_s, ACTIVE_THREADS * xzlen, RMG_MPI_ISEND, -1, 0, -1, pct.grid_comm, 15, &sreqs[15]);
        RMG_MPI_trade(xzmsms_s, ACTIVE_THREADS * xzlen, RMG_MPI_ISEND, 1, 0, 1, pct.grid_comm, 16, &sreqs[16]);
        RMG_MPI_trade(xzmsps_s, ACTIVE_THREADS * xzlen, RMG_MPI_ISEND, 1, 0, -1, pct.grid_comm, 17, &sreqs[17]);

    }

    // Load up w with the interior points
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


    // Wait for all the recvs to finish
    if(tid == 0) {
        retval = MPI_Waitall(26, rreqs, MPI_STATUSES_IGNORE);
        if(retval != MPI_SUCCESS)  error_handler("Error in MPI_Waitall.\n");
    }
    thread_barrier_wait();


    // Unpack yz-plane edges
    c1 = (dimy + images) * incy;
    c2 = (dimz + images);
    idx = tid * dimx * images * images;
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
    c1 = (dimy + images) * incy;
    c2 = (dimx + images) * incx;
    idx = tid * dimz * images * images;
    for (ix = 0; ix < images; ix++)
    {
        ixs2 = ix * incx;
        for(iy = 0;iy < images;iy++) {
            iys2 = ixs2 + iy*incy;
            for(iz = 0;iz < dimz;iz++) {
                index = iys2 + iz + images;
                w[index + c2] = xypsms_r[idx];
                w[index + c2 + c1] = xypsps_r[idx];
                w[index] =  xymsms_r[idx];
                w[index + c1] = xymsps_r[idx];
                idx++;
            }
        }
    }


    // Unpack xz-plane edges
    c1 = (dimz + images);
    c2 = (dimx + images) * incx;
    idx = tid * dimy * images * images;
    for (ix = 0; ix < images; ix++)
    {
        ixs2 = ix * incx;
        for(iy = 0;iy < dimy;iy++) {
            iys2 = ixs2 + (iy + images)*incy;
            for(iz = 0;iz < images;iz++) {
                index = iys2 + iz;
                w[index + c2] = xzpsms_r[idx];
                w[index + c2 + c1] = xzpsps_r[idx];
                w[index] =  xzmsms_r[idx];
                w[index + c1] = xzmsps_r[idx];
                idx++;
            }
        }
    }





    /* Unpack z-planes */
    c1 = dimz + images;
    idx = tid * dimx * dimy * images;
    if(images == 2) {

        for (ix = 0; ix < dimx; ix++)
        {
            ixs2 = (ix + images) * incx;
            for (iy = 0; iy < dimy; iy++)
            {
                iys2 = ixs2 + (iy + images) * incy;
                w[iys2] = frdz1n[idx];
                w[iys2+1] = frdz1n[idx+1];
                w[iys2 + c1] = frdz2n[idx];
                w[iys2 + c1 + 1] = frdz2n[idx+1];
                idx+=2;
            }                       /* end for */

        }                           /* end for */

    }
    else {

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

    }

    /* Unpack the north and south planes */
    c1 = (dimy + images) * incy;
    idx = tid * dimx * dimz * images;
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
    idx = tid * dimy * dimz * images;
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


    // Unpack the corners
    node_idx = 0;
    for(ix = -1;ix <= 1;ix+=2)
    {
        c1 = (ix == -1) ? 0 : (dimx+images)*incx;

        for(iy = -1;iy <= 1;iy+=2)
        {
            c2 = (iy == -1) ? 0 : (dimy+images)*incy;

            for(iz = -1;iz <= 1;iz+=2)
            {
                c3 = (iz == -1) ? 0 : (dimz+images);

                // unpack the recv arrays
                idx1 = node_idx*corner_node_stride + tid * img3;
                for(ix1 = 0;ix1 < images;ix1++)
                {
                    ixs2 = ix1*incx;

                    for(iy1 = 0;iy1 < images;iy1++)
                    {
                        iys2 = ixs2 + iy1*incy;

                        for(iz1 = 0;iz1 < images;iz1++)
                        {
                            index = iys2 + iz1;
                            w[index + c1 + c2 + c3] = m0_r[idx1];
                            idx1++;
                        }
                    }
                }
                node_idx++;

            }
        }
    }


    // Finally wait for all the sends to finish
    if(tid == 0) {
        retval = MPI_Waitall(26, sreqs, MPI_STATUSES_IGNORE);
        if(retval != MPI_SUCCESS)  error_handler("Error in MPI_Waitall.\n");
    }
    thread_barrier_wait();

#if MD_TIMERS
    time2 = my_crtc ();
    rmg_timings (IMAGE_TIME, (time2 - time1));
#endif

} // end trade_imagesx



// Asynchronous image trades for central finite difference operators
void trade_imagesx_central_async (rmg_double_t * f, rmg_double_t * w, int dimx, int dimy, int dimz, int images)
{
    int ix, iy, iz, incx, incy, incx0, incy0, index, tim, ione = 1;
    int ixs, iys, ixs2, iys2, c1, idx;
    int xlen, ylen, zlen;
    int tid=0, retval;
    int ACTIVE_THREADS = 1;

#if MD_TIMERS
    rmg_double_t time1, time2;
    time1 = my_crtc ();
#endif

    if(images > MAX_TRADE_IMAGES) {
        error_handler ("Images count too high in trade_imagesx_async. Modify and recompile may be required.\n");
    }

    tid = get_thread_tid();
    if(tid < 0) tid = 0;
    if(is_loop_over_states()) ACTIVE_THREADS = ct.THREADS_PER_NODE;

    tim = 2 * images;

    incx = (dimy + tim) * (dimz + tim);
    incy = dimz + tim;
    incx0 = dimy * dimz;
    incy0 = dimz;

    zlen = dimx * dimy * images;
    ylen = dimx * dimz * images;
    xlen = dimy * dimz * images;


    // Thread 0 posts all of the receives
    if(tid == 0) {

        // The z planes
        RMG_MPI_trade(frdz2n, ACTIVE_THREADS * zlen, RMG_MPI_IRECV, 0, 0, 1, pct.grid_comm, 0, &rreqs[0]);
        RMG_MPI_trade(frdz1n, ACTIVE_THREADS * zlen, RMG_MPI_IRECV, 0, 0, -1, pct.grid_comm, 1, &rreqs[1]);

        // The y planes
        RMG_MPI_trade(frdy2n, ACTIVE_THREADS * ylen, RMG_MPI_IRECV, 0, 1, 0, pct.grid_comm, 2, &rreqs[2]);
        RMG_MPI_trade(frdy1n, ACTIVE_THREADS * ylen, RMG_MPI_IRECV, 0, -1, 0, pct.grid_comm, 3, &rreqs[3]);

        // The x planes
        RMG_MPI_trade(frdx2n, ACTIVE_THREADS * xlen, RMG_MPI_IRECV, 1, 0, 0, pct.grid_comm, 4, &rreqs[4]);
        RMG_MPI_trade(frdx1n, ACTIVE_THREADS * xlen, RMG_MPI_IRECV, -1, 0, 0, pct.grid_comm, 5, &rreqs[5]);

    }


    /* Collect the positive z-plane and negative z-planes */
    c1 = (dimz-images);
    idx = tid * dimx * dimy * images;
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
    idx = tid * dimx * dimz * images;
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
    idx = tid * dimy * dimz * images;
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


    thread_barrier_wait();
    if(tid == 0) {

        // Send z planes
        RMG_MPI_trade(frdz1, ACTIVE_THREADS * zlen, RMG_MPI_ISEND, 0, 0, -1, pct.grid_comm, 0, &sreqs[0]);
        RMG_MPI_trade(frdz2, ACTIVE_THREADS * zlen, RMG_MPI_ISEND, 0, 0, 1, pct.grid_comm, 1, &sreqs[1]);

        // Send the north and south planes
        RMG_MPI_trade(frdy1, ACTIVE_THREADS * ylen, RMG_MPI_ISEND, 0, -1, 0, pct.grid_comm, 2, &sreqs[2]);
        RMG_MPI_trade(frdy2, ACTIVE_THREADS * ylen, RMG_MPI_ISEND, 0, 1, 0, pct.grid_comm, 3, &sreqs[3]);

        // Send the east and west planes
        RMG_MPI_trade(frdx1, ACTIVE_THREADS * xlen, RMG_MPI_ISEND, -1, 0, 0, pct.grid_comm, 4, &sreqs[4]);
        RMG_MPI_trade(frdx2, ACTIVE_THREADS * xlen, RMG_MPI_ISEND, 1, 0, 0, pct.grid_comm, 5, &sreqs[5]);

    }


    // Load up w with the interior points
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


    // Wait for all the recvs to finish
    if(tid == 0) {
        retval = MPI_Waitall(6, rreqs, MPI_STATUSES_IGNORE);
        if(retval != MPI_SUCCESS)  error_handler("Error in MPI_Waitall.\n");
    }
    thread_barrier_wait();


    /* Unpack z-planes */
    c1 = dimz + images;
    idx = tid * dimx * dimy * images;
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
    idx = tid * dimx * dimz * images;
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
    idx = tid * dimy * dimz * images;
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
    if(tid == 0) {
        retval = MPI_Waitall(6, sreqs, MPI_STATUSES_IGNORE);
        if(retval != MPI_SUCCESS)  error_handler("Error in MPI_Waitall.\n");
    }
    thread_barrier_wait();

#if MD_TIMERS
    time2 = my_crtc ();
    rmg_timings (IMAGE_TIME, (time2 - time1));
#endif

} // end trade_imagesx_central_async



// Used by coarse level multigrid routines which assume that the central data is already present
// and that the f array is sized [dimx+2][dimy+2][dimz+2]
//
void trade_images1_async (rmg_double_t * f, int dimx, int dimy, int dimz)
{
    int ix, iy, iz, incx, incy, index;
    int ixs2, iys2, c1, c2, c3, idx, idx1;
    int xlen, ylen, zlen, yzlen, xylen, xzlen;
    int tid=0, corner_node_stride, node_idx, retval;
    int ACTIVE_THREADS = 1;

#if MD_TIMERS
    rmg_double_t time1, time2;
    time1 = my_crtc ();
#endif

    tid = get_thread_tid();
    if(tid < 0) tid = 0;
    if(is_loop_over_states()) ACTIVE_THREADS = ct.THREADS_PER_NODE;

    corner_node_stride = ct.THREADS_PER_NODE * MAX_IMG3;

    incx = (dimy + 2) * (dimz + 2);
    incy = dimz + 2;

    zlen = dimx * dimy;
    ylen = dimx * dimz;
    xlen = dimy * dimz;

    yzlen = dimx;
    xylen = dimz;
    xzlen = dimy;

    // Thread 0 posts all of the receives
    if(tid == 0) {

        // Corner nodes
        RMG_MPI_trade( &m0_r[0*corner_node_stride], ACTIVE_THREADS, RMG_MPI_IRECV, -1, -1, -1, pct.grid_comm, 18, &rreqs[18]);
        RMG_MPI_trade( &m0_r[1*corner_node_stride], ACTIVE_THREADS, RMG_MPI_IRECV, -1, -1, 1, pct.grid_comm, 19, &rreqs[19]);
        RMG_MPI_trade( &m0_r[2*corner_node_stride], ACTIVE_THREADS, RMG_MPI_IRECV, -1, 1, -1, pct.grid_comm, 20, &rreqs[20]);
        RMG_MPI_trade( &m0_r[3*corner_node_stride], ACTIVE_THREADS, RMG_MPI_IRECV, -1, 1, 1, pct.grid_comm, 21, &rreqs[21]);
        RMG_MPI_trade( &m0_r[4*corner_node_stride], ACTIVE_THREADS, RMG_MPI_IRECV, 1, -1, -1, pct.grid_comm, 22, &rreqs[22]);
        RMG_MPI_trade( &m0_r[5*corner_node_stride], ACTIVE_THREADS, RMG_MPI_IRECV, 1, -1, 1, pct.grid_comm, 23, &rreqs[23]);
        RMG_MPI_trade( &m0_r[6*corner_node_stride], ACTIVE_THREADS, RMG_MPI_IRECV, 1, 1, -1, pct.grid_comm, 24, &rreqs[24]);
        RMG_MPI_trade( &m0_r[7*corner_node_stride], ACTIVE_THREADS, RMG_MPI_IRECV, 1, 1, 1, pct.grid_comm, 25, &rreqs[25]);

        // The z planes
        RMG_MPI_trade(frdz2n, ACTIVE_THREADS * zlen, RMG_MPI_IRECV, 0, 0, 1, pct.grid_comm, 0, &rreqs[0]);
        RMG_MPI_trade(frdz1n, ACTIVE_THREADS * zlen, RMG_MPI_IRECV, 0, 0, -1, pct.grid_comm, 1, &rreqs[1]);

        // The y planes
        RMG_MPI_trade(frdy2n, ACTIVE_THREADS * ylen, RMG_MPI_IRECV, 0, 1, 0, pct.grid_comm, 2, &rreqs[2]);
        RMG_MPI_trade(frdy1n, ACTIVE_THREADS * ylen, RMG_MPI_IRECV, 0, -1, 0, pct.grid_comm, 3, &rreqs[3]);

        // The x planes
        RMG_MPI_trade(frdx2n, ACTIVE_THREADS * xlen, RMG_MPI_IRECV, 1, 0, 0, pct.grid_comm, 4, &rreqs[4]);
        RMG_MPI_trade(frdx1n, ACTIVE_THREADS * xlen, RMG_MPI_IRECV, -1, 0, 0, pct.grid_comm, 5, &rreqs[5]);

        // And the yz-plane edges
        RMG_MPI_trade(yzpsms_r, ACTIVE_THREADS * yzlen, RMG_MPI_IRECV, 0, 1, -1, pct.grid_comm, 6, &rreqs[6]);
        RMG_MPI_trade(yzpsps_r, ACTIVE_THREADS * yzlen, RMG_MPI_IRECV, 0, 1, 1, pct.grid_comm, 7, &rreqs[7]);
        RMG_MPI_trade(yzmsms_r, ACTIVE_THREADS * yzlen, RMG_MPI_IRECV, 0, -1, -1, pct.grid_comm, 8, &rreqs[8]);
        RMG_MPI_trade(yzmsps_r, ACTIVE_THREADS * yzlen, RMG_MPI_IRECV, 0, -1, 1, pct.grid_comm, 9, &rreqs[9]);

        // And the xy-plane edges
        RMG_MPI_trade(xypsms_r, ACTIVE_THREADS * xylen, RMG_MPI_IRECV, 1, -1, 0, pct.grid_comm, 10, &rreqs[10]);
        RMG_MPI_trade(xypsps_r, ACTIVE_THREADS * xylen, RMG_MPI_IRECV, 1, 1, 0, pct.grid_comm, 11, &rreqs[11]);
        RMG_MPI_trade(xymsms_r, ACTIVE_THREADS * xylen, RMG_MPI_IRECV, -1, -1, 0, pct.grid_comm, 12, &rreqs[12]);
        RMG_MPI_trade(xymsps_r, ACTIVE_THREADS * xylen, RMG_MPI_IRECV, -1, 1, 0, pct.grid_comm, 13, &rreqs[13]);

        // And the xz-plane edges
        RMG_MPI_trade(xzpsms_r, ACTIVE_THREADS * xzlen, RMG_MPI_IRECV, 1, 0, -1, pct.grid_comm, 14, &rreqs[14]);
        RMG_MPI_trade(xzpsps_r, ACTIVE_THREADS * xzlen, RMG_MPI_IRECV, 1, 0, 1, pct.grid_comm, 15, &rreqs[15]);
        RMG_MPI_trade(xzmsms_r, ACTIVE_THREADS * xzlen, RMG_MPI_IRECV, -1, 0, -1, pct.grid_comm, 16, &rreqs[16]);
        RMG_MPI_trade(xzmsps_r, ACTIVE_THREADS * xzlen, RMG_MPI_IRECV, -1, 0, 1, pct.grid_comm, 17, &rreqs[17]);

    }

    // Collect data for the corner nodes
    node_idx = 0;
    for(ix = -1;ix <= 1;ix+=2) 
    {
        c1 = (ix == -1) ? (dimx)*incx : 1;

        for(iy = -1;iy <= 1;iy+=2) 
        {
            c2 = (iy == -1) ? (dimy)*incy : 1;

            for(iz = -1;iz <= 1;iz+=2) 
            {
                c3 = (iz == -1) ? (dimz) : 1;

                // Pack the send arrays
                idx1 = node_idx*corner_node_stride + tid;
                m0_s[idx1] = f[c1 + c2 + c3];
                node_idx++;

            }
        }
    }


    // Send them
    thread_barrier_wait();
    if(tid == 0) {

        // Corners
        RMG_MPI_trade( &m0_s[0*corner_node_stride], ACTIVE_THREADS, RMG_MPI_ISEND, 1, 1, 1, pct.grid_comm, 18, &sreqs[18]);
        RMG_MPI_trade( &m0_s[1*corner_node_stride], ACTIVE_THREADS, RMG_MPI_ISEND, 1, 1, -1, pct.grid_comm, 19, &sreqs[19]);
        RMG_MPI_trade( &m0_s[2*corner_node_stride], ACTIVE_THREADS, RMG_MPI_ISEND, 1, -1, 1, pct.grid_comm, 20, &sreqs[20]);
        RMG_MPI_trade( &m0_s[3*corner_node_stride], ACTIVE_THREADS, RMG_MPI_ISEND, 1, -1, -1, pct.grid_comm, 21, &sreqs[21]);
        RMG_MPI_trade( &m0_s[4*corner_node_stride], ACTIVE_THREADS, RMG_MPI_ISEND, -1, 1, 1, pct.grid_comm, 22, &sreqs[22]);
        RMG_MPI_trade( &m0_s[5*corner_node_stride], ACTIVE_THREADS, RMG_MPI_ISEND, -1, 1, -1, pct.grid_comm, 23, &sreqs[23]);
        RMG_MPI_trade( &m0_s[6*corner_node_stride], ACTIVE_THREADS, RMG_MPI_ISEND, -1, -1, 1, pct.grid_comm, 24, &sreqs[24]);
        RMG_MPI_trade( &m0_s[7*corner_node_stride], ACTIVE_THREADS, RMG_MPI_ISEND, -1, -1, -1, pct.grid_comm, 25, &sreqs[25]);

    }


    /* Collect the positive z-plane and negative z-planes */
    c1 = dimz-1;
    idx = tid * dimx * dimy;
    for (ix = 1; ix < dimx + 1; ix++)
    {
        ixs2 = ix * incx;
        for (iy = 1; iy < dimy + 1; iy++)
        {
            iys2 = ixs2 + iy * incy;
            for (iz = 1; iz < 2; iz++)
            {

                index = iys2 + iz;

                frdz1[idx] = f[index];
                frdz2[idx] = f[index+c1];
                idx++;

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */


    /* Collect the north and south planes */
    c1 = (dimy-1)*incy;
    idx = tid * dimx * dimz;
    for (ix = 1; ix < dimx + 1; ix++)
    {
        ixs2 = ix * incx;
        for (iy = 1; iy < 2; iy++)
        {

            iys2 = ixs2 + iy * incy;
            for (iz = 1; iz < dimz + 1; iz++)
            {

                index = iys2 + iz;

                frdy1[idx] = f[index];
                frdy2[idx] = f[index + c1];
                idx++;

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */


    /* Collect the east and west planes */
    c1 = (dimx-1)*incx;
    idx = tid * dimy * dimz;
    for (ix = 1; ix < 2; ix++)
    {
        ixs2 = ix * incx;
        for (iy = 1; iy < dimy + 1; iy++)
        {
            iys2 = ixs2 + iy * incy;
            for (iz = 1; iz < dimz + 1; iz++)
            {

                index = iys2 + iz;

                frdx1[idx] = f[index];
                frdx2[idx] = f[index + c1];
                idx++;

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */


    thread_barrier_wait();
    if(tid == 0) {

        // Send z planes
        RMG_MPI_trade(frdz1, ACTIVE_THREADS * zlen, RMG_MPI_ISEND, 0, 0, -1, pct.grid_comm, 0, &sreqs[0]);
        RMG_MPI_trade(frdz2, ACTIVE_THREADS * zlen, RMG_MPI_ISEND, 0, 0, 1, pct.grid_comm, 1, &sreqs[1]);

        // Send the north and south planes
        RMG_MPI_trade(frdy1, ACTIVE_THREADS * ylen, RMG_MPI_ISEND, 0, -1, 0, pct.grid_comm, 2, &sreqs[2]);
        RMG_MPI_trade(frdy2, ACTIVE_THREADS * ylen, RMG_MPI_ISEND, 0, 1, 0, pct.grid_comm, 3, &sreqs[3]);

        // Send the east and west planes
        RMG_MPI_trade(frdx1, ACTIVE_THREADS * xlen, RMG_MPI_ISEND, -1, 0, 0, pct.grid_comm, 4, &sreqs[4]);
        RMG_MPI_trade(frdx2, ACTIVE_THREADS * xlen, RMG_MPI_ISEND, 1, 0, 0, pct.grid_comm, 5, &sreqs[5]);

    }


    // Pack yz-plane edges
    c1 = (dimy-1)*incy;
    c2 = dimz-1;
    idx = tid * dimx;
    for (ix = 1; ix < dimx + 1; ix++)
    {
        ixs2 = ix * incx;
        for(iy = 1;iy < 2;iy++) {
            iys2 = ixs2 + iy*incy;
            for(iz = 1;iz < 2;iz++) {
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
    c1 = (dimy-1)*incy;
    c2 = (dimx-1)*incx;
    idx = tid * dimz;
    for (ix = 1; ix < 2; ix++)
    {
        ixs2 = ix * incx;
        for(iy = 1;iy < 2;iy++) {
            iys2 = ixs2 + iy*incy;
            for(iz = 1;iz < dimz + 1;iz++) {
                index = iys2 + iz;
                xypsms_s[idx] = f[index + c1]; 
                xypsps_s[idx] = f[index];
                xymsms_s[idx] = f[index + c2 + c1];
                xymsps_s[idx] = f[index + c2];
                idx++;
            }
        }
    }


    // Pack xz plane edges
    c1 = (dimz-1);
    c2 = (dimx-1)*incx;
    idx = tid * dimy;
    for (ix = 1; ix < 2; ix++)
    {
        ixs2 = ix * incx;
        for(iy = 1;iy < dimy + 1;iy++) {
            iys2 = ixs2 + iy*incy;
            for(iz = 1;iz < 2;iz++) {
                index = iys2 + iz;
                xzpsms_s[idx] = f[index + c1]; 
                xzpsps_s[idx] = f[index];
                xzmsms_s[idx] = f[index + c2 + c1];
                xzmsps_s[idx] = f[index + c2];
                idx++;
            }
        }
    }


    // Send them
    thread_barrier_wait();
    if(tid == 0) {

        // Send yz-plane edges
        RMG_MPI_trade(yzpsms_s, ACTIVE_THREADS * yzlen, RMG_MPI_ISEND, 0, -1, 1, pct.grid_comm, 6, &sreqs[6]);
        RMG_MPI_trade(yzpsps_s, ACTIVE_THREADS * yzlen, RMG_MPI_ISEND, 0, -1, -1, pct.grid_comm, 7, &sreqs[7]);
        RMG_MPI_trade(yzmsms_s, ACTIVE_THREADS * yzlen, RMG_MPI_ISEND, 0, 1, 1, pct.grid_comm, 8, &sreqs[8]);
        RMG_MPI_trade(yzmsps_s, ACTIVE_THREADS * yzlen, RMG_MPI_ISEND, 0, 1, -1, pct.grid_comm, 9, &sreqs[9]);

        // Send xy-plane edges
        RMG_MPI_trade(xypsms_s, ACTIVE_THREADS * xylen, RMG_MPI_ISEND, -1, 1, 0, pct.grid_comm, 10, &sreqs[10]);
        RMG_MPI_trade(xypsps_s, ACTIVE_THREADS * xylen, RMG_MPI_ISEND, -1, -1, 0, pct.grid_comm, 11, &sreqs[11]);
        RMG_MPI_trade(xymsms_s, ACTIVE_THREADS * xylen, RMG_MPI_ISEND, 1, 1, 0, pct.grid_comm, 12, &sreqs[12]);
        RMG_MPI_trade(xymsps_s, ACTIVE_THREADS * xylen, RMG_MPI_ISEND, 1, -1, 0, pct.grid_comm, 13, &sreqs[13]);

        // Send xz-plane edges
        RMG_MPI_trade(xzpsms_s, ACTIVE_THREADS * xzlen, RMG_MPI_ISEND, -1, 0, 1, pct.grid_comm, 14, &sreqs[14]);
        RMG_MPI_trade(xzpsps_s, ACTIVE_THREADS * xzlen, RMG_MPI_ISEND, -1, 0, -1, pct.grid_comm, 15, &sreqs[15]);
        RMG_MPI_trade(xzmsms_s, ACTIVE_THREADS * xzlen, RMG_MPI_ISEND, 1, 0, 1, pct.grid_comm, 16, &sreqs[16]);
        RMG_MPI_trade(xzmsps_s, ACTIVE_THREADS * xzlen, RMG_MPI_ISEND, 1, 0, -1, pct.grid_comm, 17, &sreqs[17]);

    }

    // Wait for all the recvs to finish
    if(tid == 0) {
        retval = MPI_Waitall(26, rreqs, MPI_STATUSES_IGNORE);
        if(retval != MPI_SUCCESS)  error_handler("Error in MPI_Waitall.\n");
    }
    thread_barrier_wait();


    // Unpack yz-plane edges
    c1 = (dimy + 1) * incy;
    c2 = (dimz + 1);
    idx = tid * dimx;
    for (ix = 0; ix < dimx; ix++)
    {
        ixs2 = (ix + 1) * incx;
        for(iy = 0;iy < 1;iy++) {
            iys2 = ixs2 + iy*incy;
            for(iz = 0;iz < 1;iz++) {
                index = iys2 + iz;
                f[index + c1] = yzpsms_r[idx];
                f[index + c2 + c1] = yzpsps_r[idx];
                f[index] =  yzmsms_r[idx];
                f[index + c2] = yzmsps_r[idx];
                idx++;
            }
        }
    }



    // Unpack xy-plane edges
    c1 = (dimy + 1) * incy;
    c2 = (dimx + 1) * incx;
    idx = tid * dimz;
    for (ix = 0; ix < 1; ix++)
    {
        ixs2 = ix * incx;
        for(iy = 0;iy < 1;iy++) {
            iys2 = ixs2 + iy*incy;
            for(iz = 0;iz < dimz;iz++) {
                index = iys2 + iz + 1;
                f[index + c2] = xypsms_r[idx];
                f[index + c2 + c1] = xypsps_r[idx];
                f[index] =  xymsms_r[idx];
                f[index + c1] = xymsps_r[idx];
                idx++;
            }
        }
    }


    // Unpack xz-plane edges
    c1 = (dimz + 1);
    c2 = (dimx + 1) * incx;
    idx = tid * dimy;
    for (ix = 0; ix < 1; ix++)
    {
        ixs2 = ix * incx;
        for(iy = 0;iy < dimy;iy++) {
            iys2 = ixs2 + (iy + 1)*incy;
            for(iz = 0;iz < 1;iz++) {
                index = iys2 + iz;
                f[index + c2] = xzpsms_r[idx];
                f[index + c2 + c1] = xzpsps_r[idx];
                f[index] =  xzmsms_r[idx];
                f[index + c1] = xzmsps_r[idx];
                idx++;
            }
        }
    }


    /* Unpack z-planes */
    c1 = dimz + 1;
    idx = tid * dimx * dimy;
    for (ix = 0; ix < dimx; ix++)
    {
        ixs2 = (ix + 1) * incx;
        for (iy = 0; iy < dimy; iy++)
        {
            iys2 = ixs2 + (iy + 1) * incy;
            for (iz = 0; iz < 1; iz++)
            {

                index = iys2 + iz;
                f[index] = frdz1n[idx];
                f[index + c1] = frdz2n[idx];
                idx++;

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */


    /* Unpack the north and south planes */
    c1 = (dimy + 1) * incy;
    idx = tid * dimx * dimz;
    for (ix = 0; ix < dimx; ix++)
    {
        ixs2 = (ix + 1) * incx;
        for (iy = 0; iy < 1; iy++)
        {
            iys2 = ixs2 + iy * incy;
            for (iz = 0; iz < dimz; iz++)
            {
                index = iys2 + iz + 1;
                f[index] = frdy1n[idx];
                f[index + c1] = frdy2n[idx];
                idx++;

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */


    /* Unpack the east and west planes */
    c1 = (dimx + 1) * incx;
    idx = tid * dimy * dimz;
    for (ix = 0; ix < 1; ix++)
    {
        ixs2 = ix * incx;
        for (iy = 0; iy < dimy; iy++)
        {
            iys2 = ixs2 + (iy + 1) * incy;
            for (iz = 0; iz < dimz; iz++)
            {

                index = iys2 + iz + 1;
                f[index] = frdx1n[idx];
                f[index + c1] = frdx2n[idx];
                idx++;

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */


    // Unpack the corners
    node_idx = 0;
    for(ix = -1;ix <= 1;ix+=2)
    {
        c1 = (ix == -1) ? 0 : (dimx+1)*incx;

        for(iy = -1;iy <= 1;iy+=2)
        {
            c2 = (iy == -1) ? 0 : (dimy+1)*incy;

            for(iz = -1;iz <= 1;iz+=2)
            {
                c3 = (iz == -1) ? 0 : (dimz+1);

                // unpack the recv arrays
                idx1 = node_idx*corner_node_stride + tid;
                f[c1 + c2 + c3] = m0_r[idx1];
                node_idx++;

            }
        }
    }


    // Finally wait for all the sends to finish
    if(tid == 0) {
        retval = MPI_Waitall(26, sreqs, MPI_STATUSES_IGNORE);
        if(retval != MPI_SUCCESS)  error_handler("Error in MPI_Waitall.\n");
    }
    thread_barrier_wait();

#if MD_TIMERS
    time2 = my_crtc ();
    rmg_timings (IMAGE_TIME, (time2 - time1));
#endif

} // end trade_images1_async



// This function is used to setup the MPI request
void RMG_MPI_trade(rmg_double_t *buf, int count, int type, int pe_x_offset, int pe_y_offset, int pe_z_offset, MPI_Comm comm, int tag, MPI_Request *req)
{
    int tid=0, ntag, target;

    // Incement offsets so they can act as array indices into send and recv lists
    pe_x_offset++;
    pe_y_offset++;
    pe_z_offset++;

    // Tag is based on tid in the lower 8 bits which gives us up to 256 threads
    tid = get_thread_tid();
    if(tid == -1) tid = 0;


    ntag = (tag<<8) + tid;
    target = target_node[pe_x_offset][pe_y_offset][pe_z_offset];

    if(type == RMG_MPI_IRECV) {
        MPI_Irecv(buf, count, MPI_DOUBLE, target,
                   ntag, comm, req);
    }
    else {
        MPI_Isend(buf, count, MPI_DOUBLE, target,
                       ntag, comm, req);
    }

}

// This function is used to setup the MPI request
void RMG_MPI_trade_f(rmg_float_t *buf, int count, int type, int pe_x_offset, int pe_y_offset, int pe_z_offset, MPI_Comm comm, int tag, MPI_Request *req)
{
    int tid=0, ntag, target;

    // Incement offsets so they can act as array indices into send and recv lists
    pe_x_offset++;
    pe_y_offset++;
    pe_z_offset++;

    // Tag is based on tid in the lower 8 bits which gives us up to 256 threads
    tid = get_thread_tid();
    if(tid == -1) tid = 0;


    ntag = (tag<<8) + tid;
    target = target_node[pe_x_offset][pe_y_offset][pe_z_offset];

    if(type == RMG_MPI_IRECV) {
        MPI_Irecv(buf, count, MPI_FLOAT, target,
                   ntag, comm, req);
    }
    else {
        MPI_Isend(buf, count, MPI_FLOAT, target,
                       ntag, comm, req);
    }

}


// Allocates memory via MPI_Alloc_mem for use on systems with RDMA capability
void init_trade_imagesx_async(void) 
{
    int retval;
    int ix, iy, iz;
    int pe_x, pe_y, pe_z;
    int t_pe_x, t_pe_y, t_pe_z;
    int grid_xp, grid_yp, grid_zp;


    printf("Using Async trade_images with max images = %d.\n", MAX_TRADE_IMAGES);

    grid_xp = get_FPX0_GRID() + 2*MAX_TRADE_IMAGES;
    grid_yp = get_FPY0_GRID() + 2*MAX_TRADE_IMAGES;
    grid_zp = get_FPZ0_GRID() + 2*MAX_TRADE_IMAGES;
    if(grid_xp > grid_yp) {
        GRID_MAX1 = grid_xp;
        if(grid_yp > grid_zp) {
            GRID_MAX2 = grid_yp;
        }
        else {
            GRID_MAX2 = grid_zp;
        }
     }
     else {
        GRID_MAX1 = grid_yp;
          if(grid_xp > grid_zp) {
              GRID_MAX2 = grid_xp;
          }
          else {
              GRID_MAX2 = grid_zp;
          }
     }

     if(GRID_MAX1 > GRID_MAX2) {
         TRADE_GRID_EDGES = GRID_MAX1;
     }
     else {
         TRADE_GRID_EDGES = GRID_MAX2;
     }

    // Set up the target node array
    pe2xyz(pct.gridpe, &pe_x, &pe_y, &pe_z);
    for(ix = -1;ix <= 1;ix++) {

        t_pe_x = (pe_x + ix + get_PE_X()) % get_PE_X();

        for(iy = -1;iy <= 1;iy++) {

            t_pe_y = (pe_y + iy + get_PE_Y()) % get_PE_Y();

            for(iz = -1;iz <= 1;iz++) {

                t_pe_z = (pe_z + iz + get_PE_Z()) % get_PE_Z();
                XYZ2PE (t_pe_x, t_pe_y, t_pe_z, target_node[ix+1][iy+1][iz+1]);

            }
        }
    } // end for


    // Allocate memory buffers using MPI_Alloc_mem
    retval = MPI_Alloc_mem(sizeof(rmg_double_t) * MAX_TRADE_IMAGES * ct.THREADS_PER_NODE * GRID_MAX1 * GRID_MAX2 , MPI_INFO_NULL, &frdx1);
    if(retval != MPI_SUCCESS) {
         error_handler("Error in MPI_Alloc_mem.\n");
    }
    frdx1_f = (rmg_float_t *)frdx1;
    retval = MPI_Alloc_mem(sizeof(rmg_double_t) * MAX_TRADE_IMAGES * ct.THREADS_PER_NODE * GRID_MAX1 * GRID_MAX2 , MPI_INFO_NULL, &frdx2);
    if(retval != MPI_SUCCESS) {
         error_handler("Error in MPI_Alloc_mem.\n");
    }
    frdx2_f = (rmg_float_t *)frdx2;
    retval = MPI_Alloc_mem(sizeof(rmg_double_t) * MAX_TRADE_IMAGES * ct.THREADS_PER_NODE * GRID_MAX1 * GRID_MAX2 , MPI_INFO_NULL, &frdy1);
    if(retval != MPI_SUCCESS) {
         error_handler("Error in MPI_Alloc_mem.\n");
    }
    frdy1_f = (rmg_float_t *)frdy1;
    retval = MPI_Alloc_mem(sizeof(rmg_double_t) * MAX_TRADE_IMAGES * ct.THREADS_PER_NODE * GRID_MAX1 * GRID_MAX2 , MPI_INFO_NULL, &frdy2);
    if(retval != MPI_SUCCESS) {
         error_handler("Error in MPI_Alloc_mem.\n");
    }
    frdy2_f = (rmg_float_t *)frdy2;
    retval = MPI_Alloc_mem(sizeof(rmg_double_t) * MAX_TRADE_IMAGES * ct.THREADS_PER_NODE * GRID_MAX1 * GRID_MAX2 , MPI_INFO_NULL, &frdz1);
    if(retval != MPI_SUCCESS) {
         error_handler("Error in MPI_Alloc_mem.\n");
    }
    frdz1_f = (rmg_float_t *)frdz1;
    retval = MPI_Alloc_mem(sizeof(rmg_double_t) * MAX_TRADE_IMAGES * ct.THREADS_PER_NODE * GRID_MAX1 * GRID_MAX2 , MPI_INFO_NULL, &frdz2);
    if(retval != MPI_SUCCESS) {
         error_handler("Error in MPI_Alloc_mem.\n");
    }
    frdz2_f = (rmg_float_t *)frdz2;
    retval = MPI_Alloc_mem(sizeof(rmg_double_t) * MAX_TRADE_IMAGES * ct.THREADS_PER_NODE * GRID_MAX1 * GRID_MAX2 , MPI_INFO_NULL, &frdx1n);
    if(retval != MPI_SUCCESS) {
         error_handler("Error in MPI_Alloc_mem.\n");
    }
    frdx1n_f = (rmg_float_t *)frdx1n;
    retval = MPI_Alloc_mem(sizeof(rmg_double_t) * MAX_TRADE_IMAGES * ct.THREADS_PER_NODE * GRID_MAX1 * GRID_MAX2 , MPI_INFO_NULL, &frdx2n);
    if(retval != MPI_SUCCESS) {
         error_handler("Error in MPI_Alloc_mem.\n");
    }
    frdx2n_f = (rmg_float_t *)frdx2n;
    retval = MPI_Alloc_mem(sizeof(rmg_double_t) * MAX_TRADE_IMAGES * ct.THREADS_PER_NODE * GRID_MAX1 * GRID_MAX2 , MPI_INFO_NULL, &frdy1n);
    if(retval != MPI_SUCCESS) {
         error_handler("Error in MPI_Alloc_mem.\n");
    }
    frdy1n_f = (rmg_float_t *)frdy1n;
    retval = MPI_Alloc_mem(sizeof(rmg_double_t) * MAX_TRADE_IMAGES * ct.THREADS_PER_NODE * GRID_MAX1 * GRID_MAX2 , MPI_INFO_NULL, &frdy2n);
    if(retval != MPI_SUCCESS) {
         error_handler("Error in MPI_Alloc_mem.\n");
    }
    frdy2n_f = (rmg_float_t *)frdy2n;
    retval = MPI_Alloc_mem(sizeof(rmg_double_t) * MAX_TRADE_IMAGES * ct.THREADS_PER_NODE * GRID_MAX1 * GRID_MAX2 , MPI_INFO_NULL, &frdz1n);
    if(retval != MPI_SUCCESS) {
         error_handler("Error in MPI_Alloc_mem.\n");
    }
    frdz1n_f = (rmg_float_t *)frdz1n;
    retval = MPI_Alloc_mem(sizeof(rmg_double_t) * MAX_TRADE_IMAGES * ct.THREADS_PER_NODE * GRID_MAX1 * GRID_MAX2 , MPI_INFO_NULL, &frdz2n);
    if(retval != MPI_SUCCESS) {
         error_handler("Error in MPI_Alloc_mem.\n");
    }
    frdz2n_f = (rmg_float_t *)frdz2n;
    retval = MPI_Alloc_mem(sizeof(rmg_double_t) * 8 * MAX_IMG2 * ct.THREADS_PER_NODE * TRADE_GRID_EDGES , MPI_INFO_NULL, &yzpsms_r);
    if(retval != MPI_SUCCESS) {
         error_handler("Error in MPI_Alloc_mem.\n");
    }
    yzpsps_r = yzpsms_r + MAX_IMG2 * ct.THREADS_PER_NODE * TRADE_GRID_EDGES;
    yzmsms_r = yzpsps_r + MAX_IMG2 * ct.THREADS_PER_NODE * TRADE_GRID_EDGES;
    yzmsps_r = yzmsms_r + MAX_IMG2 * ct.THREADS_PER_NODE * TRADE_GRID_EDGES;
    yzpsms_s = yzmsps_r + MAX_IMG2 * ct.THREADS_PER_NODE * TRADE_GRID_EDGES;
    yzpsps_s = yzpsms_s + MAX_IMG2 * ct.THREADS_PER_NODE * TRADE_GRID_EDGES;
    yzmsms_s = yzpsps_s + MAX_IMG2 * ct.THREADS_PER_NODE * TRADE_GRID_EDGES;
    yzmsps_s = yzmsms_s + MAX_IMG2 * ct.THREADS_PER_NODE * TRADE_GRID_EDGES;

    yzpsms_r_f = (rmg_float_t *)yzpsms_r;
    yzpsps_r_f = yzpsms_r_f + MAX_IMG2 * ct.THREADS_PER_NODE * TRADE_GRID_EDGES;
    yzmsms_r_f = yzpsps_r_f + MAX_IMG2 * ct.THREADS_PER_NODE * TRADE_GRID_EDGES;
    yzmsps_r_f = yzmsms_r_f + MAX_IMG2 * ct.THREADS_PER_NODE * TRADE_GRID_EDGES;
    yzpsms_s_f = yzmsps_r_f + MAX_IMG2 * ct.THREADS_PER_NODE * TRADE_GRID_EDGES;
    yzpsps_s_f = yzpsms_s_f + MAX_IMG2 * ct.THREADS_PER_NODE * TRADE_GRID_EDGES;
    yzmsms_s_f = yzpsps_s_f + MAX_IMG2 * ct.THREADS_PER_NODE * TRADE_GRID_EDGES;
    yzmsps_s_f = yzmsms_s_f + MAX_IMG2 * ct.THREADS_PER_NODE * TRADE_GRID_EDGES;

    retval = MPI_Alloc_mem(sizeof(rmg_double_t) * 8 * MAX_IMG2 * ct.THREADS_PER_NODE * TRADE_GRID_EDGES , MPI_INFO_NULL, &xypsms_r);
    if(retval != MPI_SUCCESS) {
         error_handler("Error in MPI_Alloc_mem.\n");
    }
    xypsps_r = xypsms_r + MAX_IMG2 * ct.THREADS_PER_NODE * TRADE_GRID_EDGES;
    xymsms_r = xypsps_r + MAX_IMG2 * ct.THREADS_PER_NODE * TRADE_GRID_EDGES;
    xymsps_r = xymsms_r + MAX_IMG2 * ct.THREADS_PER_NODE * TRADE_GRID_EDGES;
    xypsms_s = xymsps_r + MAX_IMG2 * ct.THREADS_PER_NODE * TRADE_GRID_EDGES;
    xypsps_s = xypsms_s + MAX_IMG2 * ct.THREADS_PER_NODE * TRADE_GRID_EDGES;
    xymsms_s = xypsps_s + MAX_IMG2 * ct.THREADS_PER_NODE * TRADE_GRID_EDGES;
    xymsps_s = xymsms_s + MAX_IMG2 * ct.THREADS_PER_NODE * TRADE_GRID_EDGES;

    xypsms_r_f = (rmg_float_t *)xypsms_r;
    xypsps_r_f = xypsms_r_f + MAX_IMG2 * ct.THREADS_PER_NODE * TRADE_GRID_EDGES;
    xymsms_r_f = xypsps_r_f + MAX_IMG2 * ct.THREADS_PER_NODE * TRADE_GRID_EDGES;
    xymsps_r_f = xymsms_r_f + MAX_IMG2 * ct.THREADS_PER_NODE * TRADE_GRID_EDGES;
    xypsms_s_f = xymsps_r_f + MAX_IMG2 * ct.THREADS_PER_NODE * TRADE_GRID_EDGES;
    xypsps_s_f = xypsms_s_f + MAX_IMG2 * ct.THREADS_PER_NODE * TRADE_GRID_EDGES;
    xymsms_s_f = xypsps_s_f + MAX_IMG2 * ct.THREADS_PER_NODE * TRADE_GRID_EDGES;
    xymsps_s_f = xymsms_s_f + MAX_IMG2 * ct.THREADS_PER_NODE * TRADE_GRID_EDGES;


    retval = MPI_Alloc_mem(sizeof(rmg_double_t) * 8 * MAX_IMG2 * ct.THREADS_PER_NODE * TRADE_GRID_EDGES , MPI_INFO_NULL, &xzpsms_r);
    if(retval != MPI_SUCCESS) {
         error_handler("Error in MPI_Alloc_mem.\n");
    }
    xzpsps_r = xzpsms_r + MAX_IMG2 * ct.THREADS_PER_NODE * TRADE_GRID_EDGES;
    xzmsms_r = xzpsps_r + MAX_IMG2 * ct.THREADS_PER_NODE * TRADE_GRID_EDGES;
    xzmsps_r = xzmsms_r + MAX_IMG2 * ct.THREADS_PER_NODE * TRADE_GRID_EDGES;
    xzpsms_s = xzmsps_r + MAX_IMG2 * ct.THREADS_PER_NODE * TRADE_GRID_EDGES;
    xzpsps_s = xzpsms_s + MAX_IMG2 * ct.THREADS_PER_NODE * TRADE_GRID_EDGES;
    xzmsms_s = xzpsps_s + MAX_IMG2 * ct.THREADS_PER_NODE * TRADE_GRID_EDGES;
    xzmsps_s = xzmsms_s + MAX_IMG2 * ct.THREADS_PER_NODE * TRADE_GRID_EDGES;

    xzpsms_r_f = (rmg_float_t *)xzpsms_r;
    xzpsps_r_f = xzpsms_r_f + MAX_IMG2 * ct.THREADS_PER_NODE * TRADE_GRID_EDGES;
    xzmsms_r_f = xzpsps_r_f + MAX_IMG2 * ct.THREADS_PER_NODE * TRADE_GRID_EDGES;
    xzmsps_r_f = xzmsms_r_f + MAX_IMG2 * ct.THREADS_PER_NODE * TRADE_GRID_EDGES;
    xzpsms_s_f = xzmsps_r_f + MAX_IMG2 * ct.THREADS_PER_NODE * TRADE_GRID_EDGES;
    xzpsps_s_f = xzpsms_s_f + MAX_IMG2 * ct.THREADS_PER_NODE * TRADE_GRID_EDGES;
    xzmsms_s_f = xzpsps_s_f + MAX_IMG2 * ct.THREADS_PER_NODE * TRADE_GRID_EDGES;
    xzmsps_s_f = xzmsms_s_f + MAX_IMG2 * ct.THREADS_PER_NODE * TRADE_GRID_EDGES;

    retval = MPI_Alloc_mem(sizeof(rmg_double_t) * 8 * MAX_IMG3 * ct.THREADS_PER_NODE , MPI_INFO_NULL, &m0_r);
    if(retval != MPI_SUCCESS) {
         error_handler("Error in MPI_Alloc_mem.\n");
    }
    m0_r_f = (rmg_float_t *)m0_r;

    retval = MPI_Alloc_mem(sizeof(rmg_double_t) * 8 * MAX_IMG3 * ct.THREADS_PER_NODE , MPI_INFO_NULL, &m0_s);
    if(retval != MPI_SUCCESS) {
         error_handler("Error in MPI_Alloc_mem.\n");
    }
    m0_s_f = (rmg_float_t *)m0_s;
}


void trade_imagesx_async_f (rmg_float_t * f, rmg_float_t * w, int dimx, int dimy, int dimz, int images)
{
    int ix, iy, iz, ix1, iy1, iz1, incx, incy, incx0, incy0, index, tim, ione = 1;
    int ixs, iys, ixs2, iys2, c1, c2, c3, idx, idx1, img3;
    int xlen, ylen, zlen, stop, yzlen, xylen, xzlen;
    int tid=0, corner_node_stride, node_idx, retval;
    int ACTIVE_THREADS = 1;

#if MD_TIMERS
    rmg_double_t time1, time2;
    time1 = my_crtc ();
#endif

    if(images > MAX_TRADE_IMAGES) {
        error_handler ("Images count too high in trade_imagesx_async. Modify and recompile may be required.\n");
    }

    tid = get_thread_tid();
    if(tid < 0) tid = 0;
    if(is_loop_over_states()) ACTIVE_THREADS = ct.THREADS_PER_NODE;

    corner_node_stride = ct.THREADS_PER_NODE * MAX_IMG3; 
    tim = 2 * images;
    img3 = images*images*images;

    incx = (dimy + tim) * (dimz + tim);
    incy = dimz + tim;
    incx0 = dimy * dimz;
    incy0 = dimz;

    zlen = dimx * dimy * images;
    ylen = dimx * dimz * images;
    xlen = dimy * dimz * images;

    yzlen = images * images * dimx;
    xylen = images * images * dimz;
    xzlen = images * images * dimy;

    // Thread 0 posts all of the receives
    if(tid == 0) {

        // Corner nodes
        RMG_MPI_trade_f( &m0_r_f[0*corner_node_stride], ACTIVE_THREADS * img3, RMG_MPI_IRECV, -1, -1, -1, pct.grid_comm, 18, &rreqs[18]);
        RMG_MPI_trade_f( &m0_r_f[1*corner_node_stride], ACTIVE_THREADS * img3, RMG_MPI_IRECV, -1, -1, 1, pct.grid_comm, 19, &rreqs[19]);
        RMG_MPI_trade_f( &m0_r_f[2*corner_node_stride], ACTIVE_THREADS * img3, RMG_MPI_IRECV, -1, 1, -1, pct.grid_comm, 20, &rreqs[20]);
        RMG_MPI_trade_f( &m0_r_f[3*corner_node_stride], ACTIVE_THREADS * img3, RMG_MPI_IRECV, -1, 1, 1, pct.grid_comm, 21, &rreqs[21]);
        RMG_MPI_trade_f( &m0_r_f[4*corner_node_stride], ACTIVE_THREADS * img3, RMG_MPI_IRECV, 1, -1, -1, pct.grid_comm, 22, &rreqs[22]);
        RMG_MPI_trade_f( &m0_r_f[5*corner_node_stride], ACTIVE_THREADS * img3, RMG_MPI_IRECV, 1, -1, 1, pct.grid_comm, 23, &rreqs[23]);
        RMG_MPI_trade_f( &m0_r_f[6*corner_node_stride], ACTIVE_THREADS * img3, RMG_MPI_IRECV, 1, 1, -1, pct.grid_comm, 24, &rreqs[24]);
        RMG_MPI_trade_f( &m0_r_f[7*corner_node_stride], ACTIVE_THREADS * img3, RMG_MPI_IRECV, 1, 1, 1, pct.grid_comm, 25, &rreqs[25]);

        // The z planes
        RMG_MPI_trade_f(frdz2n_f, ACTIVE_THREADS * zlen, RMG_MPI_IRECV, 0, 0, 1, pct.grid_comm, 0, &rreqs[0]);
        RMG_MPI_trade_f(frdz1n_f, ACTIVE_THREADS * zlen, RMG_MPI_IRECV, 0, 0, -1, pct.grid_comm, 1, &rreqs[1]);

        // The y planes
        RMG_MPI_trade_f(frdy2n_f, ACTIVE_THREADS * ylen, RMG_MPI_IRECV, 0, 1, 0, pct.grid_comm, 2, &rreqs[2]);
        RMG_MPI_trade_f(frdy1n_f, ACTIVE_THREADS * ylen, RMG_MPI_IRECV, 0, -1, 0, pct.grid_comm, 3, &rreqs[3]);

        // The x planes
        RMG_MPI_trade_f(frdx2n_f, ACTIVE_THREADS * xlen, RMG_MPI_IRECV, 1, 0, 0, pct.grid_comm, 4, &rreqs[4]);
        RMG_MPI_trade_f(frdx1n_f, ACTIVE_THREADS * xlen, RMG_MPI_IRECV, -1, 0, 0, pct.grid_comm, 5, &rreqs[5]);

        // And the yz-plane edges
        RMG_MPI_trade_f(yzpsms_r_f, ACTIVE_THREADS * yzlen, RMG_MPI_IRECV, 0, 1, -1, pct.grid_comm, 6, &rreqs[6]);
        RMG_MPI_trade_f(yzpsps_r_f, ACTIVE_THREADS * yzlen, RMG_MPI_IRECV, 0, 1, 1, pct.grid_comm, 7, &rreqs[7]);
        RMG_MPI_trade_f(yzmsms_r_f, ACTIVE_THREADS * yzlen, RMG_MPI_IRECV, 0, -1, -1, pct.grid_comm, 8, &rreqs[8]);
        RMG_MPI_trade_f(yzmsps_r_f, ACTIVE_THREADS * yzlen, RMG_MPI_IRECV, 0, -1, 1, pct.grid_comm, 9, &rreqs[9]);

        // And the xy-plane edges
        RMG_MPI_trade_f(xypsms_r_f, ACTIVE_THREADS * xylen, RMG_MPI_IRECV, 1, -1, 0, pct.grid_comm, 10, &rreqs[10]);
        RMG_MPI_trade_f(xypsps_r_f, ACTIVE_THREADS * xylen, RMG_MPI_IRECV, 1, 1, 0, pct.grid_comm, 11, &rreqs[11]);
        RMG_MPI_trade_f(xymsms_r_f, ACTIVE_THREADS * xylen, RMG_MPI_IRECV, -1, -1, 0, pct.grid_comm, 12, &rreqs[12]);
        RMG_MPI_trade_f(xymsps_r_f, ACTIVE_THREADS * xylen, RMG_MPI_IRECV, -1, 1, 0, pct.grid_comm, 13, &rreqs[13]);

        // And the xz-plane edges
        RMG_MPI_trade_f(xzpsms_r_f, ACTIVE_THREADS * xzlen, RMG_MPI_IRECV, 1, 0, -1, pct.grid_comm, 14, &rreqs[14]);
        RMG_MPI_trade_f(xzpsps_r_f, ACTIVE_THREADS * xzlen, RMG_MPI_IRECV, 1, 0, 1, pct.grid_comm, 15, &rreqs[15]);
        RMG_MPI_trade_f(xzmsms_r_f, ACTIVE_THREADS * xzlen, RMG_MPI_IRECV, -1, 0, -1, pct.grid_comm, 16, &rreqs[16]);
        RMG_MPI_trade_f(xzmsps_r_f, ACTIVE_THREADS * xzlen, RMG_MPI_IRECV, -1, 0, 1, pct.grid_comm, 17, &rreqs[17]);

    }

    // Collect data for the corner nodes
    node_idx = 0;
    for(ix = -1;ix <= 1;ix+=2) 
    {
        c1 = (ix == -1) ? (dimx-images)*incx0 : 0;

        for(iy = -1;iy <= 1;iy+=2) 
        {
            c2 = (iy == -1) ? (dimy-images)*incy0 : 0;

            for(iz = -1;iz <= 1;iz+=2) 
            {
                c3 = (iz == -1) ? (dimz-images) : 0;

                // Pack the send arrays
                idx1 = node_idx*corner_node_stride + tid * img3;
                for(ix1 = 0;ix1 < images;ix1++) 
                {
                    ixs2 = ix1*incx0;

                    for(iy1 = 0;iy1 < images;iy1++) 
                    {
                        iys2 = ixs2 + iy1*incy0;

                        for(iz1 = 0;iz1 < images;iz1++) 
                        {
                            index = iys2 + iz1;
                            m0_s_f[idx1] = f[index + c1 + c2 + c3];
                            idx1++;
                        }
                    }
                }

                node_idx++;
            }
        }
    }


    // Send them
    thread_barrier_wait();
    if(tid == 0) {

        // Corners
        RMG_MPI_trade_f( &m0_s_f[0*corner_node_stride], ACTIVE_THREADS * img3, RMG_MPI_ISEND, 1, 1, 1, pct.grid_comm, 18, &sreqs[18]);
        RMG_MPI_trade_f( &m0_s_f[1*corner_node_stride], ACTIVE_THREADS * img3, RMG_MPI_ISEND, 1, 1, -1, pct.grid_comm, 19, &sreqs[19]);
        RMG_MPI_trade_f( &m0_s_f[2*corner_node_stride], ACTIVE_THREADS * img3, RMG_MPI_ISEND, 1, -1, 1, pct.grid_comm, 20, &sreqs[20]);
        RMG_MPI_trade_f( &m0_s_f[3*corner_node_stride], ACTIVE_THREADS * img3, RMG_MPI_ISEND, 1, -1, -1, pct.grid_comm, 21, &sreqs[21]);
        RMG_MPI_trade_f( &m0_s_f[4*corner_node_stride], ACTIVE_THREADS * img3, RMG_MPI_ISEND, -1, 1, 1, pct.grid_comm, 22, &sreqs[22]);
        RMG_MPI_trade_f( &m0_s_f[5*corner_node_stride], ACTIVE_THREADS * img3, RMG_MPI_ISEND, -1, 1, -1, pct.grid_comm, 23, &sreqs[23]);
        RMG_MPI_trade_f( &m0_s_f[6*corner_node_stride], ACTIVE_THREADS * img3, RMG_MPI_ISEND, -1, -1, 1, pct.grid_comm, 24, &sreqs[24]);
        RMG_MPI_trade_f( &m0_s_f[7*corner_node_stride], ACTIVE_THREADS * img3, RMG_MPI_ISEND, -1, -1, -1, pct.grid_comm, 25, &sreqs[25]);

    }


    /* Collect the positive z-plane and negative z-planes. Unroll the most common case. */
    c1 = (dimz-images);
    idx = tid * dimx * dimy * images;
    if(images == 2) {

        for (ix = 0; ix < dimx; ix++)
        {
            ixs2 = ix * incx0;
            for (iy = 0; iy < dimy; iy++)
            {
                iys2 = ixs2 + iy * incy0;

                    frdz1_f[idx] = f[iys2];
                    frdz1_f[idx+1] = f[iys2+1];
                    frdz2_f[idx] = f[iys2+c1];
                    frdz2_f[idx+1] = f[iys2+c1+1];
                    idx+=2;


            }                       /* end for */

        }                           /* end for */

    }
    else {

        for (ix = 0; ix < dimx; ix++)
        {
            ixs2 = ix * incx0;
            for (iy = 0; iy < dimy; iy++)
            {
                iys2 = ixs2 + iy * incy0;
                for (iz = 0; iz < images; iz++)
                {

                    index = iys2 + iz;

                    frdz1_f[idx] = f[index];
                    frdz2_f[idx] = f[index+c1];
                    idx++;

                }                   /* end for */

            }                       /* end for */

        }                           /* end for */

    }

    /* Collect the north and south planes */
    c1 = (dimy - images)*incy0;
    idx = tid * dimx * dimz * images;
    for (ix = 0; ix < dimx; ix++)
    {
        ixs2 = ix * incx0;
        for (iy = 0; iy < images; iy++)
        {

            iys2 = ixs2 + iy * incy0;
            for (iz = 0; iz < dimz; iz++)
            {

                index = iys2 + iz;

                frdy1_f[idx] = f[index];
                frdy2_f[idx] = f[index + c1];
                idx++;

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */


    /* Collect the east and west planes */
    c1 = (dimx - images)*incx0;
    idx = tid * dimy * dimz * images;
    for (ix = 0; ix < images; ix++)
    {
        ixs2 = ix * incx0;
        for (iy = 0; iy < dimy; iy++)
        {
            iys2 = ixs2 + iy * incy0;
            for (iz = 0; iz < dimz; iz++)
            {

                index = iys2 + iz;

                frdx1_f[idx] = f[index];
                frdx2_f[idx] = f[index + c1];
                idx++;

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */


    thread_barrier_wait();
    if(tid == 0) {

        // Send z planes
        RMG_MPI_trade_f(frdz1_f, ACTIVE_THREADS * zlen, RMG_MPI_ISEND, 0, 0, -1, pct.grid_comm, 0, &sreqs[0]);
        RMG_MPI_trade_f(frdz2_f, ACTIVE_THREADS * zlen, RMG_MPI_ISEND, 0, 0, 1, pct.grid_comm, 1, &sreqs[1]);

        // Send the north and south planes
        RMG_MPI_trade_f(frdy1_f, ACTIVE_THREADS * ylen, RMG_MPI_ISEND, 0, -1, 0, pct.grid_comm, 2, &sreqs[2]);
        RMG_MPI_trade_f(frdy2_f, ACTIVE_THREADS * ylen, RMG_MPI_ISEND, 0, 1, 0, pct.grid_comm, 3, &sreqs[3]);

        // Send the east and west planes
        RMG_MPI_trade_f(frdx1_f, ACTIVE_THREADS * xlen, RMG_MPI_ISEND, -1, 0, 0, pct.grid_comm, 4, &sreqs[4]);
        RMG_MPI_trade_f(frdx2_f, ACTIVE_THREADS * xlen, RMG_MPI_ISEND, 1, 0, 0, pct.grid_comm, 5, &sreqs[5]);

    }


    // Pack yz-plane edges
    c1 = (dimy-images)*incy0;
    c2 = (dimz - images);
    idx = tid * dimx * images * images;
    for (ix = 0; ix < dimx; ix++)
    {
        ixs2 = ix * incx0;
        for(iy = 0;iy < images;iy++) {
            iys2 = ixs2 + iy*incy0;
            for(iz = 0;iz < images;iz++) {
                index = iys2 + iz;
                yzpsms_s_f[idx] = f[index + c2]; 
                yzpsps_s_f[idx] = f[index];
                yzmsms_s_f[idx] = f[index + c2 + c1];
                yzmsps_s_f[idx] = f[index + c1];
                idx++;
            }
        }
    }


    // Pack xy plane edges
    c1 = (dimy-images)*incy0;
    c2 = (dimx-images)*incx0;
    idx = tid * dimz * images * images;
    for (ix = 0; ix < images; ix++)
    {
        ixs2 = ix * incx0;
        for(iy = 0;iy < images;iy++) {
            iys2 = ixs2 + iy*incy0;
            for(iz = 0;iz < dimz;iz++) {
                index = iys2 + iz;
                xypsms_s_f[idx] = f[index + c1]; 
                xypsps_s_f[idx] = f[index];
                xymsms_s_f[idx] = f[index + c2 + c1];
                xymsps_s_f[idx] = f[index + c2];
                idx++;
            }
        }
    }


    // Pack xz plane edges
    c1 = (dimz-images);
    c2 = (dimx-images)*incx0;
    idx = tid * dimy * images * images;
    for (ix = 0; ix < images; ix++)
    {
        ixs2 = ix * incx0;
        for(iy = 0;iy < dimy;iy++) {
            iys2 = ixs2 + iy*incy0;
            for(iz = 0;iz < images;iz++) {
                index = iys2 + iz;
                xzpsms_s_f[idx] = f[index + c1]; 
                xzpsps_s_f[idx] = f[index];
                xzmsms_s_f[idx] = f[index + c2 + c1];
                xzmsps_s_f[idx] = f[index + c2];
                idx++;
            }
        }
    }


    // Send them
    thread_barrier_wait();
    if(tid == 0) {

        // Send yz-plane edges
        RMG_MPI_trade_f(yzpsms_s_f, ACTIVE_THREADS * yzlen, RMG_MPI_ISEND, 0, -1, 1, pct.grid_comm, 6, &sreqs[6]);
        RMG_MPI_trade_f(yzpsps_s_f, ACTIVE_THREADS * yzlen, RMG_MPI_ISEND, 0, -1, -1, pct.grid_comm, 7, &sreqs[7]);
        RMG_MPI_trade_f(yzmsms_s_f, ACTIVE_THREADS * yzlen, RMG_MPI_ISEND, 0, 1, 1, pct.grid_comm, 8, &sreqs[8]);
        RMG_MPI_trade_f(yzmsps_s_f, ACTIVE_THREADS * yzlen, RMG_MPI_ISEND, 0, 1, -1, pct.grid_comm, 9, &sreqs[9]);

        // Send xy-plane edges
        RMG_MPI_trade_f(xypsms_s_f, ACTIVE_THREADS * xylen, RMG_MPI_ISEND, -1, 1, 0, pct.grid_comm, 10, &sreqs[10]);
        RMG_MPI_trade_f(xypsps_s_f, ACTIVE_THREADS * xylen, RMG_MPI_ISEND, -1, -1, 0, pct.grid_comm, 11, &sreqs[11]);
        RMG_MPI_trade_f(xymsms_s_f, ACTIVE_THREADS * xylen, RMG_MPI_ISEND, 1, 1, 0, pct.grid_comm, 12, &sreqs[12]);
        RMG_MPI_trade_f(xymsps_s_f, ACTIVE_THREADS * xylen, RMG_MPI_ISEND, 1, -1, 0, pct.grid_comm, 13, &sreqs[13]);

        // Send xz-plane edges
        RMG_MPI_trade_f(xzpsms_s_f, ACTIVE_THREADS * xzlen, RMG_MPI_ISEND, -1, 0, 1, pct.grid_comm, 14, &sreqs[14]);
        RMG_MPI_trade_f(xzpsps_s_f, ACTIVE_THREADS * xzlen, RMG_MPI_ISEND, -1, 0, -1, pct.grid_comm, 15, &sreqs[15]);
        RMG_MPI_trade_f(xzmsms_s_f, ACTIVE_THREADS * xzlen, RMG_MPI_ISEND, 1, 0, 1, pct.grid_comm, 16, &sreqs[16]);
        RMG_MPI_trade_f(xzmsps_s_f, ACTIVE_THREADS * xzlen, RMG_MPI_ISEND, 1, 0, -1, pct.grid_comm, 17, &sreqs[17]);

    }

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
    if(tid == 0) {
        retval = MPI_Waitall(26, rreqs, MPI_STATUSES_IGNORE);
        if(retval != MPI_SUCCESS)  error_handler("Error in MPI_Waitall.\n");
    }
    thread_barrier_wait();


    // Unpack yz-plane edges
    c1 = (dimy + images) * incy;
    c2 = (dimz + images);
    idx = tid * dimx * images * images;
    for (ix = 0; ix < dimx; ix++)
    {
        ixs2 = (ix + images) * incx;
        for(iy = 0;iy < images;iy++) {
            iys2 = ixs2 + iy*incy;
            for(iz = 0;iz < images;iz++) {
                index = iys2 + iz;
                w[index + c1] = yzpsms_r_f[idx];
                w[index + c2 + c1] = yzpsps_r_f[idx];
                w[index] =  yzmsms_r_f[idx];
                w[index + c2] = yzmsps_r_f[idx];
                idx++;
            }
        }
    }



    // Unpack xy-plane edges
    c1 = (dimy + images) * incy;
    c2 = (dimx + images) * incx;
    idx = tid * dimz * images * images;
    for (ix = 0; ix < images; ix++)
    {
        ixs2 = ix * incx;
        for(iy = 0;iy < images;iy++) {
            iys2 = ixs2 + iy*incy;
            for(iz = 0;iz < dimz;iz++) {
                index = iys2 + iz + images;
                w[index + c2] = xypsms_r_f[idx];
                w[index + c2 + c1] = xypsps_r_f[idx];
                w[index] =  xymsms_r_f[idx];
                w[index + c1] = xymsps_r_f[idx];
                idx++;
            }
        }
    }


    // Unpack xz-plane edges
    c1 = (dimz + images);
    c2 = (dimx + images) * incx;
    idx = tid * dimy * images * images;
    for (ix = 0; ix < images; ix++)
    {
        ixs2 = ix * incx;
        for(iy = 0;iy < dimy;iy++) {
            iys2 = ixs2 + (iy + images)*incy;
            for(iz = 0;iz < images;iz++) {
                index = iys2 + iz;
                w[index + c2] = xzpsms_r_f[idx];
                w[index + c2 + c1] = xzpsps_r_f[idx];
                w[index] =  xzmsms_r_f[idx];
                w[index + c1] = xzmsps_r_f[idx];
                idx++;
            }
        }
    }


    /* Unpack z-planes. Unroll the most common case. */
    c1 = dimz + images;
    idx = tid * dimx * dimy * images;
    if(images == 2) {

        for (ix = 0; ix < dimx; ix++)
        {
            ixs2 = (ix + images) * incx;
            for (iy = 0; iy < dimy; iy++)
            {
                iys2 = ixs2 + (iy + images) * incy;

                    w[iys2] = frdz1n_f[idx];
                    w[iys2+1] = frdz1n_f[idx+1];
                    w[iys2 + c1] = frdz2n_f[idx];
                    w[iys2 + c1 + 1] = frdz2n_f[idx+1];
                    idx+=2;

            }                       /* end for */

        }                           /* end for */

    }
    else {

        for (ix = 0; ix < dimx; ix++)
        {
            ixs2 = (ix + images) * incx;
            for (iy = 0; iy < dimy; iy++)
            {
                iys2 = ixs2 + (iy + images) * incy;
                for (iz = 0; iz < images; iz++)
                {

                    index = iys2 + iz;
                    w[index] = frdz1n_f[idx];
                    w[index + c1] = frdz2n_f[idx];
                    idx++;

                }                   /* end for */

            }                       /* end for */

        }                           /* end for */

    }

    /* Unpack the north and south planes */
    c1 = (dimy + images) * incy;
    idx = tid * dimx * dimz * images;
    for (ix = 0; ix < dimx; ix++)
    {
        ixs2 = (ix + images) * incx;
        for (iy = 0; iy < images; iy++)
        {
            iys2 = ixs2 + iy * incy;
            for (iz = 0; iz < dimz; iz++)
            {
                index = iys2 + iz + images;
                w[index] = frdy1n_f[idx];
                w[index + c1] = frdy2n_f[idx];
                idx++;

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */


    /* Unpack the east and west planes */
    c1 = (dimx + images) * incx;
    idx = tid * dimy * dimz * images;
    for (ix = 0; ix < images; ix++)
    {
        ixs2 = ix * incx;
        for (iy = 0; iy < dimy; iy++)
        {
            iys2 = ixs2 + (iy + images) * incy;
            for (iz = 0; iz < dimz; iz++)
            {

                index = iys2 + iz + images;
                w[index] = frdx1n_f[idx];
                w[index + c1] = frdx2n_f[idx];
                idx++;

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */


    // Unpack the corners
    node_idx = 0;
    for(ix = -1;ix <= 1;ix+=2)
    {
        c1 = (ix == -1) ? 0 : (dimx+images)*incx;

        for(iy = -1;iy <= 1;iy+=2)
        {
            c2 = (iy == -1) ? 0 : (dimy+images)*incy;

            for(iz = -1;iz <= 1;iz+=2)
            {
                c3 = (iz == -1) ? 0 : (dimz+images);

                // unpack the recv arrays
                idx1 = node_idx*corner_node_stride + tid * img3;
                for(ix1 = 0;ix1 < images;ix1++)
                {
                    ixs2 = ix1*incx;

                    for(iy1 = 0;iy1 < images;iy1++)
                    {
                        iys2 = ixs2 + iy1*incy;

                        for(iz1 = 0;iz1 < images;iz1++)
                        {
                            index = iys2 + iz1;
                            w[index + c1 + c2 + c3] = m0_r_f[idx1];
                            idx1++;
                        }
                    }
                }
                node_idx++;

            }
        }
    }


    // Finally wait for all the sends to finish
    if(tid == 0) {
        retval = MPI_Waitall(26, sreqs, MPI_STATUSES_IGNORE);
        if(retval != MPI_SUCCESS)  error_handler("Error in MPI_Waitall.\n");
    }
    thread_barrier_wait();

#if MD_TIMERS
    time2 = my_crtc ();
    rmg_timings (IMAGE_TIME, (time2 - time1));
#endif

} // end trade_imagesx_async_f



// Asynchronous image trades for central finite difference operators
void trade_imagesx_central_async_f (rmg_float_t * f, rmg_float_t * w, int dimx, int dimy, int dimz, int images)
{
    int ix, iy, iz, incx, incy, incx0, incy0, index, tim, ione = 1;
    int ixs, iys, ixs2, iys2, c1, idx;
    int xlen, ylen, zlen;
    int tid=0, retval;
    int ACTIVE_THREADS = 1;

#if MD_TIMERS
    rmg_double_t time1, time2;
    time1 = my_crtc ();
#endif

    if(images > MAX_TRADE_IMAGES) {
        error_handler ("Images count too high in trade_imagesx_async. Modify and recompile may be required.\n");
    }

    tid = get_thread_tid();
    if(tid < 0) tid = 0;
    if(is_loop_over_states()) ACTIVE_THREADS = ct.THREADS_PER_NODE;

    tim = 2 * images;

    incx = (dimy + tim) * (dimz + tim);
    incy = dimz + tim;
    incx0 = dimy * dimz;
    incy0 = dimz;

    zlen = dimx * dimy * images;
    ylen = dimx * dimz * images;
    xlen = dimy * dimz * images;


    // Thread 0 posts all of the receives
    if(tid == 0) {

        // The z planes
        RMG_MPI_trade_f(frdz2n_f, ACTIVE_THREADS * zlen, RMG_MPI_IRECV, 0, 0, 1, pct.grid_comm, 0, &rreqs[0]);
        RMG_MPI_trade_f(frdz1n_f, ACTIVE_THREADS * zlen, RMG_MPI_IRECV, 0, 0, -1, pct.grid_comm, 1, &rreqs[1]);

        // The y planes
        RMG_MPI_trade_f(frdy2n_f, ACTIVE_THREADS * ylen, RMG_MPI_IRECV, 0, 1, 0, pct.grid_comm, 2, &rreqs[2]);
        RMG_MPI_trade_f(frdy1n_f, ACTIVE_THREADS * ylen, RMG_MPI_IRECV, 0, -1, 0, pct.grid_comm, 3, &rreqs[3]);

        // The x planes
        RMG_MPI_trade_f(frdx2n_f, ACTIVE_THREADS * xlen, RMG_MPI_IRECV, 1, 0, 0, pct.grid_comm, 4, &rreqs[4]);
        RMG_MPI_trade_f(frdx1n_f, ACTIVE_THREADS * xlen, RMG_MPI_IRECV, -1, 0, 0, pct.grid_comm, 5, &rreqs[5]);

    }


    /* Collect the positive z-plane and negative z-planes */
    c1 = (dimz-images);
    idx = tid * dimx * dimy * images;
    for (ix = 0; ix < dimx; ix++)
    {
        ixs2 = ix * incx0;
        for (iy = 0; iy < dimy; iy++)
        {
            iys2 = ixs2 + iy * incy0;
            for (iz = 0; iz < images; iz++)
            {

                index = iys2 + iz;

                frdz1_f[idx] = f[index];
                frdz2_f[idx] = f[index+c1];
                idx++;

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */


    /* Collect the north and south planes */
    c1 = (dimy - images)*incy0;
    idx = tid * dimx * dimz * images;
    for (ix = 0; ix < dimx; ix++)
    {
        ixs2 = ix * incx0;
        for (iy = 0; iy < images; iy++)
        {

            iys2 = ixs2 + iy * incy0;
            for (iz = 0; iz < dimz; iz++)
            {

                index = iys2 + iz;

                frdy1_f[idx] = f[index];
                frdy2_f[idx] = f[index + c1];
                idx++;

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */


    /* Collect the east and west planes */
    c1 = (dimx - images)*incx0;
    idx = tid * dimy * dimz * images;
    for (ix = 0; ix < images; ix++)
    {
        ixs2 = ix * incx0;
        for (iy = 0; iy < dimy; iy++)
        {
            iys2 = ixs2 + iy * incy0;
            for (iz = 0; iz < dimz; iz++)
            {

                index = iys2 + iz;

                frdx1_f[idx] = f[index];
                frdx2_f[idx] = f[index + c1];
                idx++;

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */


    thread_barrier_wait();
    if(tid == 0) {

        // Send z planes
        RMG_MPI_trade_f(frdz1_f, ACTIVE_THREADS * zlen, RMG_MPI_ISEND, 0, 0, -1, pct.grid_comm, 0, &sreqs[0]);
        RMG_MPI_trade_f(frdz2_f, ACTIVE_THREADS * zlen, RMG_MPI_ISEND, 0, 0, 1, pct.grid_comm, 1, &sreqs[1]);

        // Send the north and south planes
        RMG_MPI_trade_f(frdy1_f, ACTIVE_THREADS * ylen, RMG_MPI_ISEND, 0, -1, 0, pct.grid_comm, 2, &sreqs[2]);
        RMG_MPI_trade_f(frdy2_f, ACTIVE_THREADS * ylen, RMG_MPI_ISEND, 0, 1, 0, pct.grid_comm, 3, &sreqs[3]);

        // Send the east and west planes
        RMG_MPI_trade_f(frdx1_f, ACTIVE_THREADS * xlen, RMG_MPI_ISEND, -1, 0, 0, pct.grid_comm, 4, &sreqs[4]);
        RMG_MPI_trade_f(frdx2_f, ACTIVE_THREADS * xlen, RMG_MPI_ISEND, 1, 0, 0, pct.grid_comm, 5, &sreqs[5]);

    }


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
    if(tid == 0) {
        retval = MPI_Waitall(6, rreqs, MPI_STATUSES_IGNORE);
        if(retval != MPI_SUCCESS)  error_handler("Error in MPI_Waitall.\n");
    }
    thread_barrier_wait();


    /* Unpack z-planes */
    c1 = dimz + images;
    idx = tid * dimx * dimy * images;
    for (ix = 0; ix < dimx; ix++)
    {
        ixs2 = (ix + images) * incx;
        for (iy = 0; iy < dimy; iy++)
        {
            iys2 = ixs2 + (iy + images) * incy;
            for (iz = 0; iz < images; iz++)
            {

                index = iys2 + iz;
                w[index] = frdz1n_f[idx];
                w[index + c1] = frdz2n_f[idx];
                idx++;

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */


    /* Unpack the north and south planes */
    c1 = (dimy + images) * incy;
    idx = tid * dimx * dimz * images;
    for (ix = 0; ix < dimx; ix++)
    {
        ixs2 = (ix + images) * incx;
        for (iy = 0; iy < images; iy++)
        {
            iys2 = ixs2 + iy * incy;
            for (iz = 0; iz < dimz; iz++)
            {
                index = iys2 + iz + images;
                w[index] = frdy1n_f[idx];
                w[index + c1] = frdy2n_f[idx];
                idx++;

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */


    /* Unpack the east and west planes */
    c1 = (dimx + images) * incx;
    idx = tid * dimy * dimz * images;
    for (ix = 0; ix < images; ix++)
    {
        ixs2 = ix * incx;
        for (iy = 0; iy < dimy; iy++)
        {
            iys2 = ixs2 + (iy + images) * incy;
            for (iz = 0; iz < dimz; iz++)
            {

                index = iys2 + iz + images;
                w[index] = frdx1n_f[idx];
                w[index + c1] = frdx2n_f[idx];
                idx++;

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */


    // Finally wait for all the sends to finish
    if(tid == 0) {
        retval = MPI_Waitall(6, sreqs, MPI_STATUSES_IGNORE);
        if(retval != MPI_SUCCESS)  error_handler("Error in MPI_Waitall.\n");
    }
    thread_barrier_wait();

#if MD_TIMERS
    time2 = my_crtc ();
    rmg_timings (IMAGE_TIME, (time2 - time1));
#endif

} // end trade_imagesx_central_async



// Used by coarse level multigrid routines which assume that the central data is already present
// and that the f array is sized [dimx+2][dimy+2][dimz+2]
//
void trade_images1_async_f (rmg_float_t * f, int dimx, int dimy, int dimz)
{
    int ix, iy, iz, incx, incy, index;
    int ixs2, iys2, c1, c2, c3, idx, idx1;
    int xlen, ylen, zlen, yzlen, xylen, xzlen;
    int tid=0, corner_node_stride, node_idx, retval;
    int ACTIVE_THREADS = 1;

#if MD_TIMERS
    rmg_double_t time1, time2;
    time1 = my_crtc ();
#endif

    tid = get_thread_tid();
    if(tid < 0) tid = 0;
    if(is_loop_over_states()) ACTIVE_THREADS = ct.THREADS_PER_NODE;

    corner_node_stride = ct.THREADS_PER_NODE * MAX_IMG3;

    incx = (dimy + 2) * (dimz + 2);
    incy = dimz + 2;

    zlen = dimx * dimy;
    ylen = dimx * dimz;
    xlen = dimy * dimz;

    yzlen = dimx;
    xylen = dimz;
    xzlen = dimy;

    // Thread 0 posts all of the receives
    if(tid == 0) {

        // Corner nodes
        RMG_MPI_trade_f( &m0_r_f[0*corner_node_stride], ACTIVE_THREADS, RMG_MPI_IRECV, -1, -1, -1, pct.grid_comm, 18, &rreqs[18]);
        RMG_MPI_trade_f( &m0_r_f[1*corner_node_stride], ACTIVE_THREADS, RMG_MPI_IRECV, -1, -1, 1, pct.grid_comm, 19, &rreqs[19]);
        RMG_MPI_trade_f( &m0_r_f[2*corner_node_stride], ACTIVE_THREADS, RMG_MPI_IRECV, -1, 1, -1, pct.grid_comm, 20, &rreqs[20]);
        RMG_MPI_trade_f( &m0_r_f[3*corner_node_stride], ACTIVE_THREADS, RMG_MPI_IRECV, -1, 1, 1, pct.grid_comm, 21, &rreqs[21]);
        RMG_MPI_trade_f( &m0_r_f[4*corner_node_stride], ACTIVE_THREADS, RMG_MPI_IRECV, 1, -1, -1, pct.grid_comm, 22, &rreqs[22]);
        RMG_MPI_trade_f( &m0_r_f[5*corner_node_stride], ACTIVE_THREADS, RMG_MPI_IRECV, 1, -1, 1, pct.grid_comm, 23, &rreqs[23]);
        RMG_MPI_trade_f( &m0_r_f[6*corner_node_stride], ACTIVE_THREADS, RMG_MPI_IRECV, 1, 1, -1, pct.grid_comm, 24, &rreqs[24]);
        RMG_MPI_trade_f( &m0_r_f[7*corner_node_stride], ACTIVE_THREADS, RMG_MPI_IRECV, 1, 1, 1, pct.grid_comm, 25, &rreqs[25]);

        // The z planes
        RMG_MPI_trade_f(frdz2n_f, ACTIVE_THREADS * zlen, RMG_MPI_IRECV, 0, 0, 1, pct.grid_comm, 0, &rreqs[0]);
        RMG_MPI_trade_f(frdz1n_f, ACTIVE_THREADS * zlen, RMG_MPI_IRECV, 0, 0, -1, pct.grid_comm, 1, &rreqs[1]);

        // The y planes
        RMG_MPI_trade_f(frdy2n_f, ACTIVE_THREADS * ylen, RMG_MPI_IRECV, 0, 1, 0, pct.grid_comm, 2, &rreqs[2]);
        RMG_MPI_trade_f(frdy1n_f, ACTIVE_THREADS * ylen, RMG_MPI_IRECV, 0, -1, 0, pct.grid_comm, 3, &rreqs[3]);

        // The x planes
        RMG_MPI_trade_f(frdx2n_f, ACTIVE_THREADS * xlen, RMG_MPI_IRECV, 1, 0, 0, pct.grid_comm, 4, &rreqs[4]);
        RMG_MPI_trade_f(frdx1n_f, ACTIVE_THREADS * xlen, RMG_MPI_IRECV, -1, 0, 0, pct.grid_comm, 5, &rreqs[5]);

        // And the yz-plane edges
        RMG_MPI_trade_f(yzpsms_r_f, ACTIVE_THREADS * yzlen, RMG_MPI_IRECV, 0, 1, -1, pct.grid_comm, 6, &rreqs[6]);
        RMG_MPI_trade_f(yzpsps_r_f, ACTIVE_THREADS * yzlen, RMG_MPI_IRECV, 0, 1, 1, pct.grid_comm, 7, &rreqs[7]);
        RMG_MPI_trade_f(yzmsms_r_f, ACTIVE_THREADS * yzlen, RMG_MPI_IRECV, 0, -1, -1, pct.grid_comm, 8, &rreqs[8]);
        RMG_MPI_trade_f(yzmsps_r_f, ACTIVE_THREADS * yzlen, RMG_MPI_IRECV, 0, -1, 1, pct.grid_comm, 9, &rreqs[9]);

        // And the xy-plane edges
        RMG_MPI_trade_f(xypsms_r_f, ACTIVE_THREADS * xylen, RMG_MPI_IRECV, 1, -1, 0, pct.grid_comm, 10, &rreqs[10]);
        RMG_MPI_trade_f(xypsps_r_f, ACTIVE_THREADS * xylen, RMG_MPI_IRECV, 1, 1, 0, pct.grid_comm, 11, &rreqs[11]);
        RMG_MPI_trade_f(xymsms_r_f, ACTIVE_THREADS * xylen, RMG_MPI_IRECV, -1, -1, 0, pct.grid_comm, 12, &rreqs[12]);
        RMG_MPI_trade_f(xymsps_r_f, ACTIVE_THREADS * xylen, RMG_MPI_IRECV, -1, 1, 0, pct.grid_comm, 13, &rreqs[13]);

        // And the xz-plane edges
        RMG_MPI_trade_f(xzpsms_r_f, ACTIVE_THREADS * xzlen, RMG_MPI_IRECV, 1, 0, -1, pct.grid_comm, 14, &rreqs[14]);
        RMG_MPI_trade_f(xzpsps_r_f, ACTIVE_THREADS * xzlen, RMG_MPI_IRECV, 1, 0, 1, pct.grid_comm, 15, &rreqs[15]);
        RMG_MPI_trade_f(xzmsms_r_f, ACTIVE_THREADS * xzlen, RMG_MPI_IRECV, -1, 0, -1, pct.grid_comm, 16, &rreqs[16]);
        RMG_MPI_trade_f(xzmsps_r_f, ACTIVE_THREADS * xzlen, RMG_MPI_IRECV, -1, 0, 1, pct.grid_comm, 17, &rreqs[17]);

    }

    // Collect data for the corner nodes
    node_idx = 0;
    for(ix = -1;ix <= 1;ix+=2) 
    {
        c1 = (ix == -1) ? (dimx)*incx : 1;

        for(iy = -1;iy <= 1;iy+=2) 
        {
            c2 = (iy == -1) ? (dimy)*incy : 1;

            for(iz = -1;iz <= 1;iz+=2) 
            {
                c3 = (iz == -1) ? (dimz) : 1;

                // Pack the send arrays
                idx1 = node_idx*corner_node_stride + tid;
                m0_s_f[idx1] = f[c1 + c2 + c3];
                node_idx++;

            }
        }
    }


    // Send them
    thread_barrier_wait();
    if(tid == 0) {

        // Corners
        RMG_MPI_trade_f( &m0_s_f[0*corner_node_stride], ACTIVE_THREADS, RMG_MPI_ISEND, 1, 1, 1, pct.grid_comm, 18, &sreqs[18]);
        RMG_MPI_trade_f( &m0_s_f[1*corner_node_stride], ACTIVE_THREADS, RMG_MPI_ISEND, 1, 1, -1, pct.grid_comm, 19, &sreqs[19]);
        RMG_MPI_trade_f( &m0_s_f[2*corner_node_stride], ACTIVE_THREADS, RMG_MPI_ISEND, 1, -1, 1, pct.grid_comm, 20, &sreqs[20]);
        RMG_MPI_trade_f( &m0_s_f[3*corner_node_stride], ACTIVE_THREADS, RMG_MPI_ISEND, 1, -1, -1, pct.grid_comm, 21, &sreqs[21]);
        RMG_MPI_trade_f( &m0_s_f[4*corner_node_stride], ACTIVE_THREADS, RMG_MPI_ISEND, -1, 1, 1, pct.grid_comm, 22, &sreqs[22]);
        RMG_MPI_trade_f( &m0_s_f[5*corner_node_stride], ACTIVE_THREADS, RMG_MPI_ISEND, -1, 1, -1, pct.grid_comm, 23, &sreqs[23]);
        RMG_MPI_trade_f( &m0_s_f[6*corner_node_stride], ACTIVE_THREADS, RMG_MPI_ISEND, -1, -1, 1, pct.grid_comm, 24, &sreqs[24]);
        RMG_MPI_trade_f( &m0_s_f[7*corner_node_stride], ACTIVE_THREADS, RMG_MPI_ISEND, -1, -1, -1, pct.grid_comm, 25, &sreqs[25]);

    }


    /* Collect the positive z-plane and negative z-planes */
    c1 = dimz-1;
    idx = tid * dimx * dimy;
    for (ix = 1; ix < dimx + 1; ix++)
    {
        ixs2 = ix * incx;
        for (iy = 1; iy < dimy + 1; iy++)
        {
            iys2 = ixs2 + iy * incy;
            for (iz = 1; iz < 2; iz++)
            {

                index = iys2 + iz;

                frdz1_f[idx] = f[index];
                frdz2_f[idx] = f[index+c1];
                idx++;

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */


    /* Collect the north and south planes */
    c1 = (dimy-1)*incy;
    idx = tid * dimx * dimz;
    for (ix = 1; ix < dimx + 1; ix++)
    {
        ixs2 = ix * incx;
        for (iy = 1; iy < 2; iy++)
        {

            iys2 = ixs2 + iy * incy;
            for (iz = 1; iz < dimz + 1; iz++)
            {

                index = iys2 + iz;

                frdy1_f[idx] = f[index];
                frdy2_f[idx] = f[index + c1];
                idx++;

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */


    /* Collect the east and west planes */
    c1 = (dimx-1)*incx;
    idx = tid * dimy * dimz;
    for (ix = 1; ix < 2; ix++)
    {
        ixs2 = ix * incx;
        for (iy = 1; iy < dimy + 1; iy++)
        {
            iys2 = ixs2 + iy * incy;
            for (iz = 1; iz < dimz + 1; iz++)
            {

                index = iys2 + iz;

                frdx1_f[idx] = f[index];
                frdx2_f[idx] = f[index + c1];
                idx++;

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */


    thread_barrier_wait();
    if(tid == 0) {

        // Send z planes
        RMG_MPI_trade_f(frdz1_f, ACTIVE_THREADS * zlen, RMG_MPI_ISEND, 0, 0, -1, pct.grid_comm, 0, &sreqs[0]);
        RMG_MPI_trade_f(frdz2_f, ACTIVE_THREADS * zlen, RMG_MPI_ISEND, 0, 0, 1, pct.grid_comm, 1, &sreqs[1]);

        // Send the north and south planes
        RMG_MPI_trade_f(frdy1_f, ACTIVE_THREADS * ylen, RMG_MPI_ISEND, 0, -1, 0, pct.grid_comm, 2, &sreqs[2]);
        RMG_MPI_trade_f(frdy2_f, ACTIVE_THREADS * ylen, RMG_MPI_ISEND, 0, 1, 0, pct.grid_comm, 3, &sreqs[3]);

        // Send the east and west planes
        RMG_MPI_trade_f(frdx1_f, ACTIVE_THREADS * xlen, RMG_MPI_ISEND, -1, 0, 0, pct.grid_comm, 4, &sreqs[4]);
        RMG_MPI_trade_f(frdx2_f, ACTIVE_THREADS * xlen, RMG_MPI_ISEND, 1, 0, 0, pct.grid_comm, 5, &sreqs[5]);

    }


    // Pack yz-plane edges
    c1 = (dimy-1)*incy;
    c2 = dimz-1;
    idx = tid * dimx;
    for (ix = 1; ix < dimx + 1; ix++)
    {
        ixs2 = ix * incx;
        for(iy = 1;iy < 2;iy++) {
            iys2 = ixs2 + iy*incy;
            for(iz = 1;iz < 2;iz++) {
                index = iys2 + iz;
                yzpsms_s_f[idx] = f[index + c2]; 
                yzpsps_s_f[idx] = f[index];
                yzmsms_s_f[idx] = f[index + c2 + c1];
                yzmsps_s_f[idx] = f[index + c1];
                idx++;
            }
        }
    }


    // Pack xy plane edges
    c1 = (dimy-1)*incy;
    c2 = (dimx-1)*incx;
    idx = tid * dimz;
    for (ix = 1; ix < 2; ix++)
    {
        ixs2 = ix * incx;
        for(iy = 1;iy < 2;iy++) {
            iys2 = ixs2 + iy*incy;
            for(iz = 1;iz < dimz + 1;iz++) {
                index = iys2 + iz;
                xypsms_s_f[idx] = f[index + c1]; 
                xypsps_s_f[idx] = f[index];
                xymsms_s_f[idx] = f[index + c2 + c1];
                xymsps_s_f[idx] = f[index + c2];
                idx++;
            }
        }
    }


    // Pack xz plane edges
    c1 = (dimz-1);
    c2 = (dimx-1)*incx;
    idx = tid * dimy;
    for (ix = 1; ix < 2; ix++)
    {
        ixs2 = ix * incx;
        for(iy = 1;iy < dimy + 1;iy++) {
            iys2 = ixs2 + iy*incy;
            for(iz = 1;iz < 2;iz++) {
                index = iys2 + iz;
                xzpsms_s_f[idx] = f[index + c1]; 
                xzpsps_s_f[idx] = f[index];
                xzmsms_s_f[idx] = f[index + c2 + c1];
                xzmsps_s_f[idx] = f[index + c2];
                idx++;
            }
        }
    }


    // Send them
    thread_barrier_wait();
    if(tid == 0) {

        // Send yz-plane edges
        RMG_MPI_trade_f(yzpsms_s_f, ACTIVE_THREADS * yzlen, RMG_MPI_ISEND, 0, -1, 1, pct.grid_comm, 6, &sreqs[6]);
        RMG_MPI_trade_f(yzpsps_s_f, ACTIVE_THREADS * yzlen, RMG_MPI_ISEND, 0, -1, -1, pct.grid_comm, 7, &sreqs[7]);
        RMG_MPI_trade_f(yzmsms_s_f, ACTIVE_THREADS * yzlen, RMG_MPI_ISEND, 0, 1, 1, pct.grid_comm, 8, &sreqs[8]);
        RMG_MPI_trade_f(yzmsps_s_f, ACTIVE_THREADS * yzlen, RMG_MPI_ISEND, 0, 1, -1, pct.grid_comm, 9, &sreqs[9]);

        // Send xy-plane edges
        RMG_MPI_trade_f(xypsms_s_f, ACTIVE_THREADS * xylen, RMG_MPI_ISEND, -1, 1, 0, pct.grid_comm, 10, &sreqs[10]);
        RMG_MPI_trade_f(xypsps_s_f, ACTIVE_THREADS * xylen, RMG_MPI_ISEND, -1, -1, 0, pct.grid_comm, 11, &sreqs[11]);
        RMG_MPI_trade_f(xymsms_s_f, ACTIVE_THREADS * xylen, RMG_MPI_ISEND, 1, 1, 0, pct.grid_comm, 12, &sreqs[12]);
        RMG_MPI_trade_f(xymsps_s_f, ACTIVE_THREADS * xylen, RMG_MPI_ISEND, 1, -1, 0, pct.grid_comm, 13, &sreqs[13]);

        // Send xz-plane edges
        RMG_MPI_trade_f(xzpsms_s_f, ACTIVE_THREADS * xzlen, RMG_MPI_ISEND, -1, 0, 1, pct.grid_comm, 14, &sreqs[14]);
        RMG_MPI_trade_f(xzpsps_s_f, ACTIVE_THREADS * xzlen, RMG_MPI_ISEND, -1, 0, -1, pct.grid_comm, 15, &sreqs[15]);
        RMG_MPI_trade_f(xzmsms_s_f, ACTIVE_THREADS * xzlen, RMG_MPI_ISEND, 1, 0, 1, pct.grid_comm, 16, &sreqs[16]);
        RMG_MPI_trade_f(xzmsps_s_f, ACTIVE_THREADS * xzlen, RMG_MPI_ISEND, 1, 0, -1, pct.grid_comm, 17, &sreqs[17]);

    }

    // Wait for all the recvs to finish
    if(tid == 0) {
        retval = MPI_Waitall(26, rreqs, MPI_STATUSES_IGNORE);
        if(retval != MPI_SUCCESS)  error_handler("Error in MPI_Waitall.\n");
    }
    thread_barrier_wait();


    // Unpack yz-plane edges
    c1 = (dimy + 1) * incy;
    c2 = (dimz + 1);
    idx = tid * dimx;
    for (ix = 0; ix < dimx; ix++)
    {
        ixs2 = (ix + 1) * incx;
        for(iy = 0;iy < 1;iy++) {
            iys2 = ixs2 + iy*incy;
            for(iz = 0;iz < 1;iz++) {
                index = iys2 + iz;
                f[index + c1] = yzpsms_r_f[idx];
                f[index + c2 + c1] = yzpsps_r_f[idx];
                f[index] =  yzmsms_r_f[idx];
                f[index + c2] = yzmsps_r_f[idx];
                idx++;
            }
        }
    }



    // Unpack xy-plane edges
    c1 = (dimy + 1) * incy;
    c2 = (dimx + 1) * incx;
    idx = tid * dimz;
    for (ix = 0; ix < 1; ix++)
    {
        ixs2 = ix * incx;
        for(iy = 0;iy < 1;iy++) {
            iys2 = ixs2 + iy*incy;
            for(iz = 0;iz < dimz;iz++) {
                index = iys2 + iz + 1;
                f[index + c2] = xypsms_r_f[idx];
                f[index + c2 + c1] = xypsps_r_f[idx];
                f[index] =  xymsms_r_f[idx];
                f[index + c1] = xymsps_r_f[idx];
                idx++;
            }
        }
    }


    // Unpack xz-plane edges
    c1 = (dimz + 1);
    c2 = (dimx + 1) * incx;
    idx = tid * dimy;
    for (ix = 0; ix < 1; ix++)
    {
        ixs2 = ix * incx;
        for(iy = 0;iy < dimy;iy++) {
            iys2 = ixs2 + (iy + 1)*incy;
            for(iz = 0;iz < 1;iz++) {
                index = iys2 + iz;
                f[index + c2] = xzpsms_r_f[idx];
                f[index + c2 + c1] = xzpsps_r_f[idx];
                f[index] =  xzmsms_r_f[idx];
                f[index + c1] = xzmsps_r_f[idx];
                idx++;
            }
        }
    }


    /* Unpack z-planes */
    c1 = dimz + 1;
    idx = tid * dimx * dimy;
    for (ix = 0; ix < dimx; ix++)
    {
        ixs2 = (ix + 1) * incx;
        for (iy = 0; iy < dimy; iy++)
        {
            iys2 = ixs2 + (iy + 1) * incy;
            for (iz = 0; iz < 1; iz++)
            {

                index = iys2 + iz;
                f[index] = frdz1n_f[idx];
                f[index + c1] = frdz2n_f[idx];
                idx++;

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */


    /* Unpack the north and south planes */
    c1 = (dimy + 1) * incy;
    idx = tid * dimx * dimz;
    for (ix = 0; ix < dimx; ix++)
    {
        ixs2 = (ix + 1) * incx;
        for (iy = 0; iy < 1; iy++)
        {
            iys2 = ixs2 + iy * incy;
            for (iz = 0; iz < dimz; iz++)
            {
                index = iys2 + iz + 1;
                f[index] = frdy1n_f[idx];
                f[index + c1] = frdy2n_f[idx];
                idx++;

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */


    /* Unpack the east and west planes */
    c1 = (dimx + 1) * incx;
    idx = tid * dimy * dimz;
    for (ix = 0; ix < 1; ix++)
    {
        ixs2 = ix * incx;
        for (iy = 0; iy < dimy; iy++)
        {
            iys2 = ixs2 + (iy + 1) * incy;
            for (iz = 0; iz < dimz; iz++)
            {

                index = iys2 + iz + 1;
                f[index] = frdx1n_f[idx];
                f[index + c1] = frdx2n_f[idx];
                idx++;

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */


    // Unpack the corners
    node_idx = 0;
    for(ix = -1;ix <= 1;ix+=2)
    {
        c1 = (ix == -1) ? 0 : (dimx+1)*incx;

        for(iy = -1;iy <= 1;iy+=2)
        {
            c2 = (iy == -1) ? 0 : (dimy+1)*incy;

            for(iz = -1;iz <= 1;iz+=2)
            {
                c3 = (iz == -1) ? 0 : (dimz+1);

                // unpack the recv arrays
                idx1 = node_idx*corner_node_stride + tid;
                f[c1 + c2 + c3] = m0_r_f[idx1];
                node_idx++;

            }
        }
    }


    // Finally wait for all the sends to finish
    if(tid == 0) {
        retval = MPI_Waitall(26, sreqs, MPI_STATUSES_IGNORE);
        if(retval != MPI_SUCCESS)  error_handler("Error in MPI_Waitall.\n");
    }
    thread_barrier_wait();

#if MD_TIMERS
    time2 = my_crtc ();
    rmg_timings (IMAGE_TIME, (time2 - time1));
#endif

} // end trade_images1_async


// Used by coarse level multigrid routines which assume that the central data is already present
// and that the f array is sized [dimx+2][dimy+2][dimz+2]
//
void trade_images1_central_async_f (rmg_float_t * f, int dimx, int dimy, int dimz)
{
    int ix, iy, iz, incx, incy, index;
    int ixs2, iys2, c1, idx;
    int xlen, ylen, zlen;
    int tid=0, retval;
    int ACTIVE_THREADS = 1;

#if MD_TIMERS
    rmg_double_t time1, time2;
    time1 = my_crtc ();
#endif

    tid = get_thread_tid();
    if(tid < 0) tid = 0;
    if(is_loop_over_states()) ACTIVE_THREADS = ct.THREADS_PER_NODE;


    incx = (dimy + 2) * (dimz + 2);
    incy = dimz + 2;

    zlen = dimx * dimy;
    ylen = dimx * dimz;
    xlen = dimy * dimz;

    // Thread 0 posts all of the receives
    if(tid == 0) {

        // The z planes
        RMG_MPI_trade_f(frdz2n_f, ACTIVE_THREADS * zlen, RMG_MPI_IRECV, 0, 0, 1, pct.grid_comm, 0, &rreqs[0]);
        RMG_MPI_trade_f(frdz1n_f, ACTIVE_THREADS * zlen, RMG_MPI_IRECV, 0, 0, -1, pct.grid_comm, 1, &rreqs[1]);

        // The y planes
        RMG_MPI_trade_f(frdy2n_f, ACTIVE_THREADS * ylen, RMG_MPI_IRECV, 0, 1, 0, pct.grid_comm, 2, &rreqs[2]);
        RMG_MPI_trade_f(frdy1n_f, ACTIVE_THREADS * ylen, RMG_MPI_IRECV, 0, -1, 0, pct.grid_comm, 3, &rreqs[3]);

        // The x planes
        RMG_MPI_trade_f(frdx2n_f, ACTIVE_THREADS * xlen, RMG_MPI_IRECV, 1, 0, 0, pct.grid_comm, 4, &rreqs[4]);
        RMG_MPI_trade_f(frdx1n_f, ACTIVE_THREADS * xlen, RMG_MPI_IRECV, -1, 0, 0, pct.grid_comm, 5, &rreqs[5]);

    }


    /* Collect the positive z-plane and negative z-planes */
    c1 = dimz-1;
    idx = tid * dimx * dimy;
    for (ix = 1; ix < dimx + 1; ix++)
    {
        ixs2 = ix * incx;
        for (iy = 1; iy < dimy + 1; iy++)
        {
            iys2 = ixs2 + iy * incy;
            for (iz = 1; iz < 2; iz++)
            {

                index = iys2 + iz;

                frdz1_f[idx] = f[index];
                frdz2_f[idx] = f[index+c1];
                idx++;

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */


    /* Collect the north and south planes */
    c1 = (dimy-1)*incy;
    idx = tid * dimx * dimz;
    for (ix = 1; ix < dimx + 1; ix++)
    {
        ixs2 = ix * incx;
        for (iy = 1; iy < 2; iy++)
        {

            iys2 = ixs2 + iy * incy;
            for (iz = 1; iz < dimz + 1; iz++)
            {

                index = iys2 + iz;

                frdy1_f[idx] = f[index];
                frdy2_f[idx] = f[index + c1];
                idx++;

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */


    /* Collect the east and west planes */
    c1 = (dimx-1)*incx;
    idx = tid * dimy * dimz;
    for (ix = 1; ix < 2; ix++)
    {
        ixs2 = ix * incx;
        for (iy = 1; iy < dimy + 1; iy++)
        {
            iys2 = ixs2 + iy * incy;
            for (iz = 1; iz < dimz + 1; iz++)
            {

                index = iys2 + iz;

                frdx1_f[idx] = f[index];
                frdx2_f[idx] = f[index + c1];
                idx++;

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */


    thread_barrier_wait();
    if(tid == 0) {

        // Send z planes
        RMG_MPI_trade_f(frdz1_f, ACTIVE_THREADS * zlen, RMG_MPI_ISEND, 0, 0, -1, pct.grid_comm, 0, &sreqs[0]);
        RMG_MPI_trade_f(frdz2_f, ACTIVE_THREADS * zlen, RMG_MPI_ISEND, 0, 0, 1, pct.grid_comm, 1, &sreqs[1]);

        // Send the north and south planes
        RMG_MPI_trade_f(frdy1_f, ACTIVE_THREADS * ylen, RMG_MPI_ISEND, 0, -1, 0, pct.grid_comm, 2, &sreqs[2]);
        RMG_MPI_trade_f(frdy2_f, ACTIVE_THREADS * ylen, RMG_MPI_ISEND, 0, 1, 0, pct.grid_comm, 3, &sreqs[3]);

        // Send the east and west planes
        RMG_MPI_trade_f(frdx1_f, ACTIVE_THREADS * xlen, RMG_MPI_ISEND, -1, 0, 0, pct.grid_comm, 4, &sreqs[4]);
        RMG_MPI_trade_f(frdx2_f, ACTIVE_THREADS * xlen, RMG_MPI_ISEND, 1, 0, 0, pct.grid_comm, 5, &sreqs[5]);

    }

    // Wait for all the recvs to finish
    if(tid == 0) {
        retval = MPI_Waitall(6, rreqs, MPI_STATUSES_IGNORE);
        if(retval != MPI_SUCCESS)  error_handler("Error in MPI_Waitall.\n");
    }
    thread_barrier_wait();


    /* Unpack z-planes */
    c1 = dimz + 1;
    idx = tid * dimx * dimy;
    for (ix = 0; ix < dimx; ix++)
    {
        ixs2 = (ix + 1) * incx;
        for (iy = 0; iy < dimy; iy++)
        {
            iys2 = ixs2 + (iy + 1) * incy;
            for (iz = 0; iz < 1; iz++)
            {

                index = iys2 + iz;
                f[index] = frdz1n_f[idx];
                f[index + c1] = frdz2n_f[idx];
                idx++;

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */


    /* Unpack the north and south planes */
    c1 = (dimy + 1) * incy;
    idx = tid * dimx * dimz;
    for (ix = 0; ix < dimx; ix++)
    {
        ixs2 = (ix + 1) * incx;
        for (iy = 0; iy < 1; iy++)
        {
            iys2 = ixs2 + iy * incy;
            for (iz = 0; iz < dimz; iz++)
            {
                index = iys2 + iz + 1;
                f[index] = frdy1n_f[idx];
                f[index + c1] = frdy2n_f[idx];
                idx++;

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */


    /* Unpack the east and west planes */
    c1 = (dimx + 1) * incx;
    idx = tid * dimy * dimz;
    for (ix = 0; ix < 1; ix++)
    {
        ixs2 = ix * incx;
        for (iy = 0; iy < dimy; iy++)
        {
            iys2 = ixs2 + (iy + 1) * incy;
            for (iz = 0; iz < dimz; iz++)
            {

                index = iys2 + iz + 1;
                f[index] = frdx1n_f[idx];
                f[index + c1] = frdx2n_f[idx];
                idx++;

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */


    // Finally wait for all the sends to finish
    if(tid == 0) {
        retval = MPI_Waitall(6, sreqs, MPI_STATUSES_IGNORE);
        if(retval != MPI_SUCCESS)  error_handler("Error in MPI_Waitall.\n");
    }
    thread_barrier_wait();

#if MD_TIMERS
    time2 = my_crtc ();
    rmg_timings (IMAGE_TIME, (time2 - time1));
#endif

} // end trade_images1_central_async_f


// Used by coarse level multigrid routines which assume that the central data is already present
// and that the f array is sized [dimx+2][dimy+2][dimz+2]
//
void trade_images1_central_async (rmg_double_t * f, int dimx, int dimy, int dimz)
{
    int ix, iy, iz, incx, incy, index;
    int ixs2, iys2, c1, idx;
    int xlen, ylen, zlen;
    int tid=0, retval;
    int ACTIVE_THREADS = 1;

#if MD_TIMERS
    rmg_double_t time1, time2;
    time1 = my_crtc ();
#endif

    tid = get_thread_tid();
    if(tid < 0) tid = 0;
    if(is_loop_over_states()) ACTIVE_THREADS = ct.THREADS_PER_NODE;


    incx = (dimy + 2) * (dimz + 2);
    incy = dimz + 2;

    zlen = dimx * dimy;
    ylen = dimx * dimz;
    xlen = dimy * dimz;

    // Thread 0 posts all of the receives
    if(tid == 0) {

        // The z planes
        RMG_MPI_trade(frdz2n, ACTIVE_THREADS * zlen, RMG_MPI_IRECV, 0, 0, 1, pct.grid_comm, 0, &rreqs[0]);
        RMG_MPI_trade(frdz1n, ACTIVE_THREADS * zlen, RMG_MPI_IRECV, 0, 0, -1, pct.grid_comm, 1, &rreqs[1]);

        // The y planes
        RMG_MPI_trade(frdy2n, ACTIVE_THREADS * ylen, RMG_MPI_IRECV, 0, 1, 0, pct.grid_comm, 2, &rreqs[2]);
        RMG_MPI_trade(frdy1n, ACTIVE_THREADS * ylen, RMG_MPI_IRECV, 0, -1, 0, pct.grid_comm, 3, &rreqs[3]);

        // The x planes
        RMG_MPI_trade(frdx2n, ACTIVE_THREADS * xlen, RMG_MPI_IRECV, 1, 0, 0, pct.grid_comm, 4, &rreqs[4]);
        RMG_MPI_trade(frdx1n, ACTIVE_THREADS * xlen, RMG_MPI_IRECV, -1, 0, 0, pct.grid_comm, 5, &rreqs[5]);

    }


    /* Collect the positive z-plane and negative z-planes */
    c1 = dimz-1;
    idx = tid * dimx * dimy;
    for (ix = 1; ix < dimx + 1; ix++)
    {
        ixs2 = ix * incx;
        for (iy = 1; iy < dimy + 1; iy++)
        {
            iys2 = ixs2 + iy * incy;
            for (iz = 1; iz < 2; iz++)
            {

                index = iys2 + iz;

                frdz1[idx] = f[index];
                frdz2[idx] = f[index+c1];
                idx++;

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */


    /* Collect the north and south planes */
    c1 = (dimy-1)*incy;
    idx = tid * dimx * dimz;
    for (ix = 1; ix < dimx + 1; ix++)
    {
        ixs2 = ix * incx;
        for (iy = 1; iy < 2; iy++)
        {

            iys2 = ixs2 + iy * incy;
            for (iz = 1; iz < dimz + 1; iz++)
            {

                index = iys2 + iz;

                frdy1[idx] = f[index];
                frdy2[idx] = f[index + c1];
                idx++;

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */


    /* Collect the east and west planes */
    c1 = (dimx-1)*incx;
    idx = tid * dimy * dimz;
    for (ix = 1; ix < 2; ix++)
    {
        ixs2 = ix * incx;
        for (iy = 1; iy < dimy + 1; iy++)
        {
            iys2 = ixs2 + iy * incy;
            for (iz = 1; iz < dimz + 1; iz++)
            {

                index = iys2 + iz;

                frdx1[idx] = f[index];
                frdx2[idx] = f[index + c1];
                idx++;

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */


    thread_barrier_wait();
    if(tid == 0) {

        // Send z planes
        RMG_MPI_trade(frdz1, ACTIVE_THREADS * zlen, RMG_MPI_ISEND, 0, 0, -1, pct.grid_comm, 0, &sreqs[0]);
        RMG_MPI_trade(frdz2, ACTIVE_THREADS * zlen, RMG_MPI_ISEND, 0, 0, 1, pct.grid_comm, 1, &sreqs[1]);

        // Send the north and south planes
        RMG_MPI_trade(frdy1, ACTIVE_THREADS * ylen, RMG_MPI_ISEND, 0, -1, 0, pct.grid_comm, 2, &sreqs[2]);
        RMG_MPI_trade(frdy2, ACTIVE_THREADS * ylen, RMG_MPI_ISEND, 0, 1, 0, pct.grid_comm, 3, &sreqs[3]);

        // Send the east and west planes
        RMG_MPI_trade(frdx1, ACTIVE_THREADS * xlen, RMG_MPI_ISEND, -1, 0, 0, pct.grid_comm, 4, &sreqs[4]);
        RMG_MPI_trade(frdx2, ACTIVE_THREADS * xlen, RMG_MPI_ISEND, 1, 0, 0, pct.grid_comm, 5, &sreqs[5]);

    }

    // Wait for all the recvs to finish
    if(tid == 0) {
        retval = MPI_Waitall(6, rreqs, MPI_STATUSES_IGNORE);
        if(retval != MPI_SUCCESS)  error_handler("Error in MPI_Waitall.\n");
    }
    thread_barrier_wait();


    /* Unpack z-planes */
    c1 = dimz + 1;
    idx = tid * dimx * dimy;
    for (ix = 0; ix < dimx; ix++)
    {
        ixs2 = (ix + 1) * incx;
        for (iy = 0; iy < dimy; iy++)
        {
            iys2 = ixs2 + (iy + 1) * incy;
            for (iz = 0; iz < 1; iz++)
            {

                index = iys2 + iz;
                f[index] = frdz1n[idx];
                f[index + c1] = frdz2n[idx];
                idx++;

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */


    /* Unpack the north and south planes */
    c1 = (dimy + 1) * incy;
    idx = tid * dimx * dimz;
    for (ix = 0; ix < dimx; ix++)
    {
        ixs2 = (ix + 1) * incx;
        for (iy = 0; iy < 1; iy++)
        {
            iys2 = ixs2 + iy * incy;
            for (iz = 0; iz < dimz; iz++)
            {
                index = iys2 + iz + 1;
                f[index] = frdy1n[idx];
                f[index + c1] = frdy2n[idx];
                idx++;

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */


    /* Unpack the east and west planes */
    c1 = (dimx + 1) * incx;
    idx = tid * dimy * dimz;
    for (ix = 0; ix < 1; ix++)
    {
        ixs2 = ix * incx;
        for (iy = 0; iy < dimy; iy++)
        {
            iys2 = ixs2 + (iy + 1) * incy;
            for (iz = 0; iz < dimz; iz++)
            {

                index = iys2 + iz + 1;
                f[index] = frdx1n[idx];
                f[index + c1] = frdx2n[idx];
                idx++;

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */


    // Finally wait for all the sends to finish
    if(tid == 0) {
        retval = MPI_Waitall(6, sreqs, MPI_STATUSES_IGNORE);
        if(retval != MPI_SUCCESS)  error_handler("Error in MPI_Waitall.\n");
    }
    thread_barrier_wait();

#if MD_TIMERS
    time2 = my_crtc ();
    rmg_timings (IMAGE_TIME, (time2 - time1));
#endif

} // end trade_images1_central_async

#endif
