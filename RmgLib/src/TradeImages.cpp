/*
 *
 * Copyright (c) 1995, Emil Briggs
 * Copyright (C) 1998  Emil Briggs, Charles Brabec, Mark Wensell, 
 *                     Dan Sullivan, Chris Rapcewicz, Jerzy Bernholc
 * Copyright (C) 2001  Emil Briggs, Wenchang Lu,
 *                     Marco Buongiorno Nardelli,Charles Brabec, 
 *                     Mark Wensell,Dan Sullivan, Chris Rapcewicz,
 *                     Jerzy Bernholc
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the <organization> nor the
 *       names of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.

 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * 
*/

#include "TradeImages.h"
#include "RmgTimer.h"
#include "MpiQueue.h"
#include <cmath>
#include <complex>
#include <mutex>
#include <condition_variable>
#include <boost/lockfree/queue.hpp>


#define ASYNC_MODE      0
#define SYNC_MODE       1


// Force instantiation of float, double and std::complex versions.
template void TradeImages::trade_images<float>(float*, int, int, int, int);
template void TradeImages::trade_images<double>(double*, int, int, int, int);
template void TradeImages::trade_images<std::complex<double> >(std::complex<double>*, int, int, int, int);
template void TradeImages::trade_images<std::complex<float> >(std::complex<float>*, int, int, int, int);
template void TradeImages::trade_imagesx<float>(float*, float*, int, int, int, int, int);
template void TradeImages::trade_imagesx<double>(double*, double*, int, int, int, int, int);
template void TradeImages::trade_imagesx<std::complex<float> >(std::complex <float>*, std::complex <float>*, int, int, int, int, int);
template void TradeImages::trade_imagesx<std::complex<double> >(std::complex <double>*, std::complex <double>*, int, int, int, int, int);



/*
 * INPUTS
 * f[dimx*dimy*dimz] - raw array without images. pack_ptos should not be called, this is 
 *                     handled inside the function now
 
 * w[(dimx+2*images) * (dimy+2*images) * (dimz+2*images)]
 *    This array is completely filled, i.e. the original data is filled in and then 
 *    the image data are added*/


// Constructor
TradeImages::TradeImages(BaseGrid *BG, size_t elem_len, bool new_queue_mode, MpiQueue *newQM)
{

    BaseThread *T = BaseThread::getBaseThread(0);
    this->G = BG;
    this->queue_mode = new_queue_mode;
    this->queue = newQM;

    int grid_xp, grid_yp, grid_zp, grid_max1, grid_max2;

    grid_xp = this->G->get_PX0_GRID(this->G->get_default_FG_RATIO()) + 2*MAX_TRADE_IMAGES;
    grid_yp = this->G->get_PY0_GRID(this->G->get_default_FG_RATIO()) + 2*MAX_TRADE_IMAGES;
    grid_zp = this->G->get_PZ0_GRID(this->G->get_default_FG_RATIO()) + 2*MAX_TRADE_IMAGES;
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


     TradeImages::max_alloc = 6 * MAX_TRADE_IMAGES * grid_max1 * grid_max2 * T->get_threads_per_node();
     int retval = MPI_Alloc_mem(sizeof(std::complex<double>) * this->max_alloc , MPI_INFO_NULL, &swbuf1x);
     if(retval != MPI_SUCCESS) {
         rmg_error_handler (__FILE__, __LINE__, "Error in MPI_Alloc_mem.\n");
     }
     retval = MPI_Alloc_mem(sizeof(std::complex<double>)* this->max_alloc , MPI_INFO_NULL, &swbuf2x);
     if(retval != MPI_SUCCESS) {
         rmg_error_handler (__FILE__, __LINE__, "Error in MPI_Alloc_mem.\n");
     }

     this->local_mode = (this->G->get_PE_X() == 1) && (this->G->get_PE_Y() == 1) && (this->G->get_PE_Z() == 1);
     TradeImages::mode = ASYNC_MODE;
     TradeImages::timer_mode = false;
     TradeImages::init_trade_imagesx_async(elem_len);        

     return;

}

// Destructor to free mem allocated via MPI_Alloc
// Not correct but since this object lives for program duration not a big deal.
TradeImages::~TradeImages(void)
{
    MPI_Free_mem(swbuf1x);
    MPI_Free_mem(swbuf2x);
    MPI_Free_mem(frdx1[0]);
    MPI_Free_mem(frdx2[0]);
    MPI_Free_mem(frdy1[0]);
    MPI_Free_mem(frdy2[0]);
    MPI_Free_mem(frdz1[0]);
    MPI_Free_mem(frdz2[0]);
    MPI_Free_mem(frdx1n[0]);
    MPI_Free_mem(frdx2n[0]);
    MPI_Free_mem(frdy1n[0]);
    MPI_Free_mem(frdy2n[0]);
    MPI_Free_mem(frdz1n[0]);
    MPI_Free_mem(frdz2n[0]);
    delete [] frdx1;
    delete [] frdx2;
    delete [] frdy1;
    delete [] frdy2;
    delete [] frdz1;
    delete [] frdz2;
    delete [] frdx1n;
    delete [] frdx2n;
    delete [] frdy1n;
    delete [] frdy2n;
    delete [] frdz1n;
    delete [] frdz2n;

    MPI_Free_mem(yzpsms_r);
    MPI_Free_mem(xypsms_r);
    MPI_Free_mem(xzpsms_r);
    MPI_Free_mem(m0_r);
    MPI_Free_mem(m0_s);
}


void TradeImages::set_timer_mode(bool verbose)
{
    TradeImages::timer_mode = verbose;
}
void TradeImages::set_queue_mode(bool mode)
{
    BaseThread *T = BaseThread::getBaseThread(0);
    int threads = T->get_threads_per_node();

    // It makes no sense to use queue mode with only a single thread.
    if(threads <= 1)
    {
        TradeImages::queue_mode = false;
        return;
    }
    TradeImages::queue_mode = mode;
}
void TradeImages::set_synchronous_mode(void)
{
    TradeImages::mode = SYNC_MODE;
}
void TradeImages::set_asynchronous_mode(void)
{
    TradeImages::mode = ASYNC_MODE;
}
void TradeImages::set_MPI_comm(MPI_Comm comm)
{
    TradeImages::comm = comm;
}
MPI_Comm TradeImages::get_MPI_comm(void)
{
    return TradeImages::comm;
}


template <typename RmgType>
void TradeImages::trade_imagesx (RmgType *f, RmgType *w, int dimx, int dimy, int dimz, int images, int type)
{
    RmgTimer *RT=NULL;
    if(this->timer_mode) RT = new RmgTimer("Trade images: trade_imagesx");

    if(this->local_mode && (type == CENTRAL_TRADE))
    {
        TradeImages::trade_imagesx_central_local (f, w, dimx, dimy, dimz, images);
        if(this->timer_mode) delete RT;
        return;
    }

    int ix, iy, iz, incx, incy, incx0, incy0, index, tim;
    int ixs, iys, ixs2, iys2, c1, c2, alloc;
    int xlen, ylen, zlen, stop;
    int *nb_ids, factor;
    MPI_Status mrstatus;
    RmgType *frdx1, *frdx2, *frdy1, *frdy2, *frdz1, *frdz2;
    RmgType *frdx1n, *frdx2n, *frdy1n, *frdy2n, *frdz1n, *frdz2n;
    RmgType *swbuf1x_f, *swbuf2x_f;
    BaseThread *T = BaseThread::getBaseThread(0);

    factor = sizeof(RmgType);

    int testsize = (dimx + 2)*(dimx + 2) + (dimy + 2)*(dimy + 2) * (dimz + 2)*(dimz + 2);

    if(TradeImages::mode == ASYNC_MODE) {
        if(type == CENTRAL_TRADE) {
            if(this->queue_mode && T->is_loop_over_states() && (testsize >= 1092))
            {
                TradeImages::trade_imagesx_central_async_managed (f, w, dimx, dimy, dimz, images);
            }
            else
            {
                TradeImages::trade_imagesx_central_async (f, w, dimx, dimy, dimz, images);
            }
        }
        else {
            if(this->queue_mode && T->is_loop_over_states() && (testsize >= 1092))
            {
                //TradeImages::trade_imagesx_async_managed (f, w, dimx, dimy, dimz, images);
                printf("Variant not programmed yet.\n");fflush(NULL);exit(0);
            }
            else
            {
                TradeImages::trade_imagesx_async (f, w, dimx, dimy, dimz, images);
            }
        }
        if(this->timer_mode) delete RT;
        return;
    }

    int ACTIVE_THREADS = 1;
    int tid = T->get_thread_tid();
    if(tid < 0) tid = 0;
    if(T->is_loop_over_states()) ACTIVE_THREADS = T->barrier->barrier_count();

    swbuf1x_f = (RmgType *)TradeImages::swbuf1x;
    swbuf2x_f = (RmgType *)TradeImages::swbuf2x;

    tim = 2 * images;

    incx = (dimy + tim) * (dimz + tim);
    incy = dimz + tim;
    incx0 = dimy * dimz;
    incy0 = dimz;

    zlen = dimx * dimy * images;
    ylen = dimx * images * (dimz + tim);
    xlen = images * (dimy + tim) * (dimz + tim);


    alloc = 2 * (xlen + ylen + zlen) * T->get_threads_per_node();
    // Verify that the required memory will fit into the statically allocated storage
    if(alloc > TradeImages::max_alloc)
        rmg_error_handler (__FILE__, __LINE__, "Not enough memory. This should never happen.");


    nb_ids = this->G->get_neighbors();

    /* Load up w with the basic stuff */
    for (ix = 0; ix < dimx; ix++)
    {
        ixs = ix * incx0;
        ixs2 = (ix + images) * incx;

        for (iy = 0; iy < dimy; iy++)
        {
            iys = ixs + iy * incy0;
            iys2 = ixs2 + (iy + images) * incy;

            for(int idx = 0;idx < dimz;idx++)
                  w[iys2 + images + idx] = f[iys + idx];

        }                       /* end for */

    }                           /* end for */


    /* Collect the positive z-plane and negative z-planes */
    frdz1 = (RmgType *)&swbuf1x_f[tid * zlen];
    frdz2 = (RmgType *)&swbuf1x_f[tid * zlen + ACTIVE_THREADS * zlen];
    frdz2n = (RmgType *)&swbuf2x_f[tid * zlen];
    frdz1n = (RmgType *)&swbuf2x_f[tid * zlen + ACTIVE_THREADS * zlen];
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


    stop = zlen * ACTIVE_THREADS;
    T->thread_barrier_wait(true);
    if(tid == 0) {
        MPI_Sendrecv (&swbuf1x_f[0], factor*stop, MPI_BYTE, nb_ids[NB_D], (1<<16),
                  &swbuf2x_f[0], factor*stop, MPI_BYTE, nb_ids[NB_U], (1<<16), TradeImages::comm, &mrstatus);

        MPI_Sendrecv (&swbuf1x_f[stop], factor*stop, MPI_BYTE, nb_ids[NB_U], (2<<16),
                  &swbuf2x_f[stop], factor*stop, MPI_BYTE, nb_ids[NB_D], (2<<16), TradeImages::comm, &mrstatus);
    }
    T->thread_barrier_wait(true);

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
    frdy1 = (RmgType *)&swbuf1x_f[tid * ylen];
    frdy2 = (RmgType *)&swbuf1x_f[tid * ylen + ACTIVE_THREADS * ylen];
    frdy2n = (RmgType *)&swbuf2x_f[tid * ylen];
    frdy1n = (RmgType *)&swbuf2x_f[tid * ylen + ACTIVE_THREADS * ylen];
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


    stop = ylen * ACTIVE_THREADS;
    T->thread_barrier_wait(true);
    if(tid == 0) {
        MPI_Sendrecv (&swbuf1x_f[0], factor*stop, MPI_BYTE, nb_ids[NB_S], (3<<16),
                  &swbuf2x_f[0], factor*stop, MPI_BYTE, nb_ids[NB_N], (3<<16), TradeImages::comm, &mrstatus);

        MPI_Sendrecv (&swbuf1x_f[stop], factor*stop, MPI_BYTE, nb_ids[NB_N], (4<<16),
                  &swbuf2x_f[stop], factor*stop, MPI_BYTE, nb_ids[NB_S], (4<<16), TradeImages::comm, &mrstatus);
    }
    T->thread_barrier_wait(true);

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
    frdx1 = (RmgType *)&swbuf1x_f[tid * xlen];
    frdx2 = (RmgType *)&swbuf1x_f[tid * xlen + ACTIVE_THREADS * xlen];
    frdx2n = (RmgType *)&swbuf2x_f[tid * xlen];
    frdx1n = (RmgType *)&swbuf2x_f[tid * xlen + ACTIVE_THREADS * xlen];
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


    stop = xlen * ACTIVE_THREADS;
    T->thread_barrier_wait(true);
    if(tid == 0) {

        MPI_Sendrecv (&swbuf1x_f[0], factor*stop, MPI_BYTE, nb_ids[NB_W], (5<<16),
                   &swbuf2x_f[0], factor*stop, MPI_BYTE, nb_ids[NB_E], (5<<16), TradeImages::comm, &mrstatus);

        MPI_Sendrecv (&swbuf1x_f[stop], factor*stop, MPI_BYTE, nb_ids[NB_E], (6<<16),
                  &swbuf2x_f[stop], factor*stop, MPI_BYTE, nb_ids[NB_W], (6<<16), TradeImages::comm, &mrstatus);

    }
    T->thread_barrier_wait(true);


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

    if(this->timer_mode) delete RT;

} // end trade_imagesx

// Local trade images when the object is only defined on one MPI process
template <typename RmgType>
void TradeImages::trade_imagesx_central_local (RmgType * f, RmgType * w, int dimx, int dimy, int dimz, int images)
{

    int tim = 2 * images;
    int incx = (dimy + tim) * (dimz + tim);
    int incy = dimz + tim;
    int incx0 = dimy * dimz;
    int incy0 = dimz;

    /* Load up w with the basic stuff */
    for (int ix = 0; ix < dimx; ix++)
    {
        int ixs = ix * incx0;
        int ixs2 = (ix + images) * incx;

        for (int iy = 0; iy < dimy; iy++)
        {
            int iys = ixs + iy * incy0;
            int iys2 = ixs2 + (iy + images) * incy;

            for(int idx = 0;idx < dimz;idx++)
                  w[iys2 + images + idx] = f[iys + idx];

        }                       /* end for */

    }                           /* end for */


    /* positive z-plane and negative z-planes */
    int c1 = dimz + images;
    for (int ix = 0; ix < dimx; ix++)
    {
        int ixs2 = (ix + images) * incx;
        for (int iy = 0; iy < dimy; iy++)
        {
            int iys2 = ixs2 + (iy + images) * incy;
            for (int iz = 0; iz < images; iz++)
            {
                int index = iys2 + iz;
                w[index + c1] = w[index + images];
                w[index] = w[index + dimz];

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */

    /* north and south planes */
    c1 = images * incy;
    int c2 = dimy * incy;
    int c3 = (dimy + images) * incy;
    for (int ix = 0; ix < dimx; ix++)
    {
        int ixs2 = (ix + images) * incx;
        for (int iy = 0; iy < images; iy++)
        {

            int iys2 = ixs2 + iy * incy;
            for (int iz = 0; iz < dimz + tim; iz++)
            {

                int index = iys2 + iz;
                w[index + c3] = w[index + c1];
                w[index] = w[index + c2];

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */


    c1 = images * incx;
    c2 = dimx * incx;
    c3 = (dimx + images) * incx;
    for (int ix = 0; ix < images; ix++)
    {
        int ixs2 = ix * incx;
        for (int iy = 0; iy < dimy + tim; iy++)
        {
            int iys2 = ixs2 + iy * incy;
            for (int iz = 0; iz < dimz + tim; iz++)
            {

                int index = iys2 + iz;
                w[index + c3]  = w[index + c1];
                w[index]  = w[index + c2];

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */


}

template <typename RmgType>
void TradeImages::trade_images (RmgType * mat, int dimx, int dimy, int dimz, int type)
{
    RmgTimer *RT=NULL;
    if(this->timer_mode) RT = new RmgTimer("Trade images: trade_images");

    if(this->local_mode)
    {
        TradeImages::trade_images_local(mat, dimx, dimy, dimz, type);
        if(this->timer_mode) delete RT;
        return;
    }

    int i, j;
    int incx, incy, incz;
    int xmax, ymax, zmax;
    int ix, iz, iys2;
    int alloc, alloc1, toffset;
    int ACTIVE_THREADS = 1, factor;
    int *nb_ids;
    BaseThread *T = BaseThread::getBaseThread(0);

    factor = sizeof(RmgType);

    MPI_Status mstatus;
    int idx, stop;
    RmgType *nmat1, *nmat2;
    RmgType *swbuf1x_f, *swbuf2x_f;


    swbuf1x_f = (RmgType *)TradeImages::swbuf1x;
    swbuf2x_f = (RmgType *)TradeImages::swbuf2x;

    nb_ids = this->G->get_neighbors();

    int testsize = (dimx + 2)*(dimx + 2) + (dimy + 2)*(dimy + 2) * (dimz + 2)*(dimz + 2);
    if(TradeImages::mode == ASYNC_MODE) {
        if(type == CENTRAL_TRADE)
        {
            if(this->queue_mode && T->is_loop_over_states() && (testsize > 1092))
            {
                TradeImages::trade_images1_central_async_managed (mat, dimx, dimy, dimz);
            }
            else
            {
                TradeImages::trade_images1_central_async (mat, dimx, dimy, dimz);
            }
            if(this->timer_mode) delete RT;
            return;
        }
        else if(this->queue_mode && T->is_loop_over_states() && (testsize > 1092))
        {
            TradeImages::trade_images_async_managed (mat, dimx, dimy, dimz);
            if(this->timer_mode) delete RT;
            return;
        }

    }

    alloc = 4 * (dimx + 2) * (dimy + 2);
    alloc1 = 4 * (dimy + 2) * (dimz + 2);
    if(alloc1 > alloc)
        alloc = alloc1;
    alloc = alloc * T->get_threads_per_node();
    if(alloc > TradeImages::max_alloc)
        rmg_error_handler (__FILE__, __LINE__, "Not enough memory. This should never happen.");


    int tid = T->get_thread_tid();
    if(tid < 0) tid = 0;
    if(T->is_loop_over_states()) ACTIVE_THREADS = T->barrier->barrier_count();

/* precalc some boundaries */
    incx = (dimy + 2) * (dimz + 2);
    incy = (dimz + 2);
    incz = 1;
    xmax = dimx * incx;
    ymax = dimy * incy;
    zmax = dimz * incz;

// For they hybrid case we use a single array for all threads which lets us combine multiple
// Mpi requests into one
    nmat1 = &swbuf1x_f[dimx * dimy * tid];
    nmat2 = &swbuf2x_f[dimx * dimy * tid];

/*
 * First, do Up-Down Trade (+/- z)
 *  (trading dimx * dimy array of data, twice)
 */


    idx = 0;
    for (i = incx; i <= incx * dimx; i += incx)
    {
        for (j = 1; j <= dimy; j++)
        {
            nmat1[idx] = mat[i + j * incy + zmax];
            idx++;
        }
    }

    // We only want one thread to do the MPI call here.
    stop = dimx * dimy * ACTIVE_THREADS;
    T->thread_barrier_wait(true);
    if(tid == 0) {
        MPI_Sendrecv (swbuf1x_f, factor*stop, MPI_BYTE, nb_ids[NB_U], (1<<16), swbuf2x_f, factor*stop,
		  MPI_BYTE, nb_ids[NB_D], (1<<16), TradeImages::comm, &mstatus);
    } 
    T->thread_barrier_wait(true);

    idx = 0;
    for (i = incx; i <= incx * dimx; i += incx)
    {
        for (j = 1; j <= dimy; j++)
        {

            mat[i + j * incy] = nmat2[idx];
            idx++;
        }
    }


    idx = 0;
    for (i = incx; i <= incx * dimx; i += incx)
    {
        for (j = 1; j <= dimy; j++)
        {
            nmat1[idx] = mat[i + j * incy + incz];
            idx++;
        }
    }


    stop = dimx * dimy * ACTIVE_THREADS;
    T->thread_barrier_wait(true);
   if(tid == 0) {
        MPI_Sendrecv (swbuf1x_f, factor*stop, MPI_BYTE, nb_ids[NB_D], (2<<16), swbuf2x_f, factor*stop,
                      MPI_BYTE, nb_ids[NB_U], (2<<16), TradeImages::comm, &mstatus);
    }
    T->thread_barrier_wait(true);


    idx = 0;
    for (i = incx; i <= incx * dimx; i += incx)
    {
        for (j = 1; j <= dimy; j++)
        {
            mat[i + j * incy + zmax + incz] = nmat2[idx];
            idx++;
        }
    }

/*
 * next, do North-South Trade (+/- y)
 *  (trading dimx * (dimz+2) array of data, twice)
 */

    stop =  dimx * (dimz+2);
    toffset = tid*dimx*(dimz+2);
    idx = 0;
    for (ix = 0; ix < dimx; ix++) {
        iys2 = (ix + 1) * incx;
        for (iz = 0; iz < dimz + 2; iz++) {
            swbuf1x_f[idx + toffset] = mat[iys2 + ymax + iz];
            idx++;
        }
    }

    T->thread_barrier_wait(true);
    if(tid == 0) {
        MPI_Sendrecv (swbuf1x_f, factor*stop * ACTIVE_THREADS, MPI_BYTE, nb_ids[NB_N], (3<<16), swbuf2x_f, factor*stop * ACTIVE_THREADS,
                  MPI_BYTE, nb_ids[NB_S], (3<<16), TradeImages::comm, &mstatus);
    }
    T->thread_barrier_wait(true);

    idx = 0;
    for (ix = 0; ix < dimx; ix++) {
        iys2 = (ix + 1) * incx;
        for (iz = 0; iz < dimz + 2; iz++) {
            mat[iys2 + iz] = swbuf2x_f[idx + toffset];
            idx++;
        }
    }

    idx = 0;
    for (ix = 0; ix < dimx; ix++) {
        iys2 = (ix + 1) * incx;
        for (iz = 0; iz < dimz + 2; iz++) {
            swbuf1x_f[idx + toffset] = mat[iys2 + incy + iz];
            idx++;
        }
    }


    T->thread_barrier_wait(true);
    if(tid == 0) {
        MPI_Sendrecv (swbuf1x_f, factor*stop * ACTIVE_THREADS, MPI_BYTE, nb_ids[NB_S], (4<<16), swbuf2x_f, factor*stop * ACTIVE_THREADS,
                  MPI_BYTE, nb_ids[NB_N], (4<<16), TradeImages::comm, &mstatus);
    }
    T->thread_barrier_wait(true);


    idx = 0;
    for (ix = 0; ix < dimx; ix++) {
        iys2 = (ix + 1) * incx;
        for (iz = 0; iz < dimz + 2; iz++) {
            mat[iys2 + ymax + incy + iz] = swbuf2x_f[idx + toffset];
            idx++;
        }
    }


/*
 * Finally, do East - West Trade (+/- x)
 *  (trading (dimy+2) * (dimz+2) array of data, twice)
 */

    stop = (dimy + 2) * (dimz + 2);

    for(int idx = 0;idx < stop;idx++)
        swbuf1x_f[tid * stop + idx] = mat[xmax + idx];



    T->thread_barrier_wait(true);
    if(tid == 0) {
        MPI_Sendrecv (swbuf1x_f, factor*stop * ACTIVE_THREADS, MPI_BYTE, nb_ids[NB_E], (5<<16), swbuf2x_f, factor*stop * ACTIVE_THREADS,
                  MPI_BYTE, nb_ids[NB_W], (5<<16), TradeImages::comm, &mstatus);
    }
    T->thread_barrier_wait(true);

    for(int idx = 0;idx < stop;idx++)
        mat[idx] = swbuf2x_f[tid * stop + idx];

    for(int idx = 0;idx < stop;idx++)
       swbuf1x_f[tid * stop + idx] = mat[incx + idx];

    T->thread_barrier_wait(true);
    if(tid == 0) {
        MPI_Sendrecv (swbuf1x_f, factor*stop * ACTIVE_THREADS, MPI_BYTE, nb_ids[NB_W], (6<<16), swbuf2x_f, factor*stop * ACTIVE_THREADS,
                  MPI_BYTE, nb_ids[NB_E], (6<<16), TradeImages::comm, &mstatus);
    }
    T->thread_barrier_wait(true);
    for(int idx = 0;idx < stop;idx++)
        mat[xmax + incx + idx] = swbuf2x_f[tid * stop + idx];


    /* For clusters set the boundaries to zero -- this is wrong for the hartree
     * potential but we'll fix it up later. */
//    if ((ct.boundaryflag == CLUSTER) || (ct.boundaryflag == SURFACE))
//        set_bc (mat, dimx, dimy, dimz, 1, 0.0);

    if(this->timer_mode) delete RT;
}



template <typename RmgType>
void TradeImages::trade_images_local (RmgType * mat, int dimx, int dimy, int dimz, int type)
{

    /* precalc some boundaries */
    int incx = (dimy + 2) * (dimz + 2);
    int incy = (dimz + 2);
    int incz = 1;
    int xmax = dimx * incx;
    int ymax = dimy * incy;
    int zmax = dimz * incz;

    for (int i = incx; i <= incx * dimx; i += incx)
    {
        for (int j = 1; j <= dimy; j++)
        {
            mat[i + j * incy] = mat[i + j * incy + zmax];
            mat[i + j * incy + zmax + incz] = mat[i + j * incy + incz];
        }
    }

    for (int ix = 0; ix < dimx; ix++) {
        int iys2 = (ix + 1) * incx;
        for (int iz = 0; iz < dimz + 2; iz++) {
            mat[iys2 + ymax + incy + iz] = mat[iys2 + incy + iz];
            mat[iys2 + iz] = mat[iys2 + ymax + iz];
        }
    }

    int stop = (dimy + 2) * (dimz + 2);
    for(int idx = 0;idx < stop;idx++)
       mat[idx] = mat[xmax + idx];

    for(int idx = 0;idx < stop;idx++)
       mat[xmax + incx + idx] = mat[incx + idx];


}

int TRADE_GRID_EDGES;
int GRID_MAX1;
int GRID_MAX2;



#define MAX_IMG2 (MAX_TRADE_IMAGES*MAX_TRADE_IMAGES)
#define MAX_IMG3 (MAX_TRADE_IMAGES*MAX_TRADE_IMAGES*MAX_TRADE_IMAGES)



// This function is used to setup the MPI request
template <typename RmgType>
void TradeImages::RMG_MPI_trade(RmgType *buf, int count, int type, int pe_x_offset, int pe_y_offset, int pe_z_offset, MPI_Comm comm, int tag, int istate, MPI_Request *req)
{
    int ntag, target, factor;

    factor = sizeof(RmgType);

    // Incement offsets so they can act as array indices into send and recv lists
    pe_x_offset++;
    pe_y_offset++;
    pe_z_offset++;

    // The actual tag passed to mpi consists of the state index shifted left 5 bits and
    // the passed tag parameter. Most MPI implementations use 32 bits for the tag so this
    // gives 2^27 bits for the state index. Crays use 2^21 bits for the tag which gives
    // 2^16 states while the MPI spec only requires 2^16 bits which gives 2^11 states.
    ntag = (istate<<5) + tag;

    target = TradeImages::target_node[pe_x_offset][pe_y_offset][pe_z_offset];

    if(type == RMG_MPI_IRECV) {
        MPI_Irecv(buf, factor*count, MPI_BYTE, target,
                   ntag, comm, req);
    }
    else {
        MPI_Isend(buf, factor*count, MPI_BYTE, target,
                       ntag, comm, req);
    }

}

// This function is used to queue an MPI trade request when running in managed mode
template <typename RmgType>
void TradeImages::RMG_MPI_queue_trade(RmgType *buf, int count, int type, int pe_x_offset, int pe_y_offset, int pe_z_offset, MPI_Comm comm, int tag, int istate, mpi_queue_item_t &qitem)
{

    qitem.comm = comm;
    qitem.is_unpacked = false;
    qitem.is_completed->store(false);
    BaseThread *T = BaseThread::getBaseThread(0);
    int tid = T->get_thread_tid();
    // The actual tag passed to mpi consists of the state index shifted left 5 bits and
    // the passed tag parameter. Most MPI implementations use 32 bits for the tag so this
    // gives 2^27 bits for the state index. Crays use 2^21 bits for the tag which gives
    // 2^16 states while the MPI spec only requires 2^16 bits which gives 2^11 states.
    qitem.mpi_tag = (istate<<5) + tag;
    qitem.target = TradeImages::target_node[pe_x_offset+1][pe_y_offset+1][pe_z_offset+1];
    qitem.target_index = (pe_x_offset+1)*9 + (pe_y_offset+1)*3 + pe_z_offset+1;
    qitem.type = type;
    qitem.buf = (void *)buf;
    qitem.buflen = sizeof(RmgType)*count;

    // Push it onto the queue
    this->queue->push(tid, qitem);

}

template <typename RmgType>
void TradeImages::RMG_MPI_queue_allreduce(RmgType *buf, int count, MPI_Datatype datatype, MPI_Comm comm, mpi_queue_item_t &qitem)
{

    // Tag is based on tid in the lower 8 bits which gives us up to 256 threads
    qitem.comm = comm;
    qitem.is_completed = false;
    qitem.is_unpacked = false;
    qitem.type = RMG_MPI_SUM;
    qitem.datatype = datatype;
    qitem.buf = (void *)buf;
    qitem.buflen = sizeof(RmgType)*count;

    // Push it onto the queue
//    this->queue->push(qitem);

}

void TradeImages::allocate_buffers(double ** &P, int nthreads, int length_per_thread, size_t elem_len)
{
    int factor = elem_len / sizeof(double);
    P = new double * [nthreads];
    int retval = MPI_Alloc_mem(elem_len * nthreads * length_per_thread , MPI_INFO_NULL, &P[0]);
    if(retval != MPI_SUCCESS) {
        rmg_error_handler (__FILE__, __LINE__, "Error in MPI_Alloc_mem.\n");
    }
    for(int tid = 1;tid < nthreads;tid++) P[tid] = P[tid-1] + factor * length_per_thread;
#if 0
    for(int tid = 1;tid < nthreads;tid++)
    {
        int retval = MPI_Alloc_mem(elem_len * length_per_thread , MPI_INFO_NULL, &P[tid]);
        if(retval != MPI_SUCCESS) {
            rmg_error_handler (__FILE__, __LINE__, "Error in MPI_Alloc_mem.\n");
        }
    }
#endif
}

// Allocates memory via MPI_Alloc_mem for use on systems with RDMA capability
void TradeImages::init_trade_imagesx_async(size_t elem_len) 
{
    int THREADS_PER_NODE;
    int ix, iy, iz;
    int pe_x, pe_y, pe_z;
    int t_pe_x, t_pe_y, t_pe_z;
    int grid_xp, grid_yp, grid_zp;
    BaseThread *T = BaseThread::getBaseThread(0);

    THREADS_PER_NODE = T->get_threads_per_node();

    //printf("Using Async trade_images with max images = %d.\n", MAX_TRADE_IMAGES);

    grid_xp = this->G->get_PX0_GRID(this->G->get_default_FG_RATIO()) + 2*MAX_TRADE_IMAGES;
    grid_yp = this->G->get_PY0_GRID(this->G->get_default_FG_RATIO()) + 2*MAX_TRADE_IMAGES;
    grid_zp = this->G->get_PZ0_GRID(this->G->get_default_FG_RATIO()) + 2*MAX_TRADE_IMAGES;
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
    this->G->pe2xyz(this->G->get_rank(), &pe_x, &pe_y, &pe_z);
    for(ix = -1;ix <= 1;ix++) {

        t_pe_x = (pe_x + ix + this->G->get_PE_X()) % this->G->get_PE_X();

        for(iy = -1;iy <= 1;iy++) {

            t_pe_y = (pe_y + iy + this->G->get_PE_Y()) % this->G->get_PE_Y();

            for(iz = -1;iz <= 1;iz++) {

                t_pe_z = (pe_z + iz + this->G->get_PE_Z()) % this->G->get_PE_Z();
                TradeImages::target_node[ix+1][iy+1][iz+1] = t_pe_x*this->G->get_PE_Y()*this->G->get_PE_Z() + t_pe_y*this->G->get_PE_Z() + t_pe_z;

            }
        }
    } // end for

    // Allocate memory buffers using MPI_Alloc_mem
    TradeImages::allocate_buffers(frdx1, THREADS_PER_NODE, MAX_TRADE_IMAGES * GRID_MAX1 * GRID_MAX2, elem_len);
    TradeImages::allocate_buffers(frdx2, THREADS_PER_NODE, MAX_TRADE_IMAGES * GRID_MAX1 * GRID_MAX2, elem_len);
    TradeImages::allocate_buffers(frdy1, THREADS_PER_NODE, MAX_TRADE_IMAGES * GRID_MAX1 * GRID_MAX2, elem_len);
    TradeImages::allocate_buffers(frdy2, THREADS_PER_NODE, MAX_TRADE_IMAGES * GRID_MAX1 * GRID_MAX2, elem_len);
    TradeImages::allocate_buffers(frdz1, THREADS_PER_NODE, MAX_TRADE_IMAGES * GRID_MAX1 * GRID_MAX2, elem_len);
    TradeImages::allocate_buffers(frdz2, THREADS_PER_NODE, MAX_TRADE_IMAGES * GRID_MAX1 * GRID_MAX2, elem_len);
    TradeImages::allocate_buffers(frdx1n, THREADS_PER_NODE, MAX_TRADE_IMAGES * GRID_MAX1 * GRID_MAX2, elem_len);
    TradeImages::allocate_buffers(frdx2n, THREADS_PER_NODE, MAX_TRADE_IMAGES * GRID_MAX1 * GRID_MAX2, elem_len);
    TradeImages::allocate_buffers(frdy1n, THREADS_PER_NODE, MAX_TRADE_IMAGES * GRID_MAX1 * GRID_MAX2, elem_len);
    TradeImages::allocate_buffers(frdy2n, THREADS_PER_NODE, MAX_TRADE_IMAGES * GRID_MAX1 * GRID_MAX2, elem_len);
    TradeImages::allocate_buffers(frdz1n, THREADS_PER_NODE, MAX_TRADE_IMAGES * GRID_MAX1 * GRID_MAX2, elem_len);
    TradeImages::allocate_buffers(frdz2n, THREADS_PER_NODE, MAX_TRADE_IMAGES * GRID_MAX1 * GRID_MAX2, elem_len);


    //retval = MPI_Alloc_mem(elem_len * 8 * MAX_IMG2 * THREADS_PER_NODE * TRADE_GRID_EDGES , MPI_INFO_NULL, &yzpsms_r);
    TradeImages::allocate_buffers(yzpsms_r, THREADS_PER_NODE, MAX_IMG2 * TRADE_GRID_EDGES, elem_len);
    TradeImages::allocate_buffers(yzpsps_r, THREADS_PER_NODE, MAX_IMG2 * TRADE_GRID_EDGES, elem_len);
    TradeImages::allocate_buffers(yzmsms_r, THREADS_PER_NODE, MAX_IMG2 * TRADE_GRID_EDGES, elem_len);
    TradeImages::allocate_buffers(yzmsps_r, THREADS_PER_NODE, MAX_IMG2 * TRADE_GRID_EDGES, elem_len);
    TradeImages::allocate_buffers(yzpsms_s, THREADS_PER_NODE, MAX_IMG2 * TRADE_GRID_EDGES, elem_len);
    TradeImages::allocate_buffers(yzpsps_s, THREADS_PER_NODE, MAX_IMG2 * TRADE_GRID_EDGES, elem_len);
    TradeImages::allocate_buffers(yzmsms_s, THREADS_PER_NODE, MAX_IMG2 * TRADE_GRID_EDGES, elem_len);
    TradeImages::allocate_buffers(yzmsps_s, THREADS_PER_NODE, MAX_IMG2 * TRADE_GRID_EDGES, elem_len);

    //retval = MPI_Alloc_mem(elem_len * 8 * MAX_IMG2 * THREADS_PER_NODE * TRADE_GRID_EDGES , MPI_INFO_NULL, &xypsms_r);
    TradeImages::allocate_buffers(xypsms_r, THREADS_PER_NODE, MAX_IMG2 * TRADE_GRID_EDGES, elem_len);
    TradeImages::allocate_buffers(xypsps_r, THREADS_PER_NODE, MAX_IMG2 * TRADE_GRID_EDGES, elem_len);
    TradeImages::allocate_buffers(xymsms_r, THREADS_PER_NODE, MAX_IMG2 * TRADE_GRID_EDGES, elem_len);
    TradeImages::allocate_buffers(xymsps_r, THREADS_PER_NODE, MAX_IMG2 * TRADE_GRID_EDGES, elem_len);
    TradeImages::allocate_buffers(xypsms_s, THREADS_PER_NODE, MAX_IMG2 * TRADE_GRID_EDGES, elem_len);
    TradeImages::allocate_buffers(xypsps_s, THREADS_PER_NODE, MAX_IMG2 * TRADE_GRID_EDGES, elem_len);
    TradeImages::allocate_buffers(xymsms_s, THREADS_PER_NODE, MAX_IMG2 * TRADE_GRID_EDGES, elem_len);
    TradeImages::allocate_buffers(xymsps_s, THREADS_PER_NODE, MAX_IMG2 * TRADE_GRID_EDGES, elem_len);

    //retval = MPI_Alloc_mem(elem_len * 8 * MAX_IMG2 * THREADS_PER_NODE * TRADE_GRID_EDGES , MPI_INFO_NULL, &xzpsms_r);
    TradeImages::allocate_buffers(xzpsms_r, THREADS_PER_NODE, MAX_IMG2 * TRADE_GRID_EDGES, elem_len);
    TradeImages::allocate_buffers(xzpsps_r, THREADS_PER_NODE, MAX_IMG2 * TRADE_GRID_EDGES, elem_len);
    TradeImages::allocate_buffers(xzmsms_r, THREADS_PER_NODE, MAX_IMG2 * TRADE_GRID_EDGES, elem_len);
    TradeImages::allocate_buffers(xzmsps_r, THREADS_PER_NODE, MAX_IMG2 * TRADE_GRID_EDGES, elem_len);
    TradeImages::allocate_buffers(xzpsms_s, THREADS_PER_NODE, MAX_IMG2 * TRADE_GRID_EDGES, elem_len);
    TradeImages::allocate_buffers(xzpsps_s, THREADS_PER_NODE, MAX_IMG2 * TRADE_GRID_EDGES, elem_len);
    TradeImages::allocate_buffers(xzmsms_s, THREADS_PER_NODE, MAX_IMG2 * TRADE_GRID_EDGES, elem_len);
    TradeImages::allocate_buffers(xzmsps_s, THREADS_PER_NODE, MAX_IMG2 * TRADE_GRID_EDGES, elem_len);


    //retval = MPI_Alloc_mem(elem_len * 8 * MAX_IMG3 * THREADS_PER_NODE , MPI_INFO_NULL, &m0_r);
    TradeImages::allocate_buffers(m0_r, THREADS_PER_NODE, 8*MAX_IMG3, elem_len);

    //retval = MPI_Alloc_mem(elem_len * 8 * MAX_IMG3 * THREADS_PER_NODE , MPI_INFO_NULL, &m0_s);
    TradeImages::allocate_buffers(m0_s, THREADS_PER_NODE, 8*MAX_IMG3, elem_len);

}

template <typename RmgType>
void TradeImages::trade_imagesx_async (RmgType * f, RmgType * w, int dimx, int dimy, int dimz, int images)
{
    BaseThread *T = BaseThread::getBaseThread(0);
    int tid = T->get_thread_tid();
    if(tid < 0) tid = 0;
    int ACTIVE_THREADS = 1;
    if(T->is_loop_over_states()) ACTIVE_THREADS = T->barrier->barrier_count();

    int istate = T->get_thread_basetag();
    MPI_Comm grid_comm = T->get_unique_comm(istate);

    int ix, iy, iz, ix1, iy1, iz1, incx, incy, incx0, incy0, index, tim;
    int ixs, iys, ixs2, iys2, c1, c2, c3, idx, idx1, img3;
    int xlen, ylen, zlen, yzlen, xylen, xzlen;
    int corner_node_stride, node_idx, retval;

    RmgType *frdx1_f, *frdx2_f, *frdy1_f, *frdy2_f, *frdz1_f, *frdz2_f;
    RmgType *frdx1n_f, *frdx2n_f, *frdy1n_f, *frdy2n_f, *frdz1n_f, *frdz2n_f;
    RmgType *yzpsms_r_f, *yzpsps_r_f, *yzmsms_r_f, *yzmsps_r_f;
    RmgType *yzpsms_s_f, *yzpsps_s_f, *yzmsms_s_f, *yzmsps_s_f;
    RmgType *xzpsms_r_f, *xzpsps_r_f, *xzmsms_r_f, *xzmsps_r_f;
    RmgType *xzpsms_s_f, *xzpsps_s_f, *xzmsms_s_f, *xzmsps_s_f;
    RmgType *xypsms_r_f, *xypsps_r_f, *xymsms_r_f, *xymsps_r_f;
    RmgType *xypsms_s_f, *xypsps_s_f, *xymsms_s_f, *xymsps_s_f;
    RmgType *m0_s_f, *m0_r_f;

    frdx1_f = (RmgType *)TradeImages::frdx1[0];
    frdx2_f = (RmgType *)TradeImages::frdx2[0];
    frdy1_f = (RmgType *)TradeImages::frdy1[0];
    frdy2_f = (RmgType *)TradeImages::frdy2[0];
    frdz1_f = (RmgType *)TradeImages::frdz1[0];
    frdz2_f = (RmgType *)TradeImages::frdz2[0];
    frdx1n_f = (RmgType *)TradeImages::frdx1n[0];
    frdx2n_f = (RmgType *)TradeImages::frdx2n[0];
    frdy1n_f = (RmgType *)TradeImages::frdy1n[0];
    frdy2n_f = (RmgType *)TradeImages::frdy2n[0];
    frdz1n_f = (RmgType *)TradeImages::frdz1n[0];
    frdz2n_f = (RmgType *)TradeImages::frdz2n[0];
    yzpsms_r_f = (RmgType *)TradeImages::yzpsms_r[0];
    yzpsps_r_f = (RmgType *)TradeImages::yzpsps_r[0];
    yzmsms_r_f = (RmgType *)TradeImages::yzmsms_r[0];
    yzmsps_r_f = (RmgType *)TradeImages::yzmsps_r[0];
    yzpsms_s_f = (RmgType *)TradeImages::yzpsms_s[0];
    yzpsps_s_f = (RmgType *)TradeImages::yzpsps_s[0];
    yzmsms_s_f = (RmgType *)TradeImages::yzmsms_s[0];
    yzmsps_s_f = (RmgType *)TradeImages::yzmsps_s[0];
    xzpsms_r_f = (RmgType *)TradeImages::xzpsms_r[0];
    xzpsps_r_f = (RmgType *)TradeImages::xzpsps_r[0];
    xzmsms_r_f = (RmgType *)TradeImages::xzmsms_r[0];
    xzmsps_r_f = (RmgType *)TradeImages::xzmsps_r[0];
    xzpsms_s_f = (RmgType *)TradeImages::xzpsms_s[0];
    xzpsps_s_f = (RmgType *)TradeImages::xzpsps_s[0];
    xzmsms_s_f = (RmgType *)TradeImages::xzmsms_s[0];
    xzmsps_s_f = (RmgType *)TradeImages::xzmsps_s[0];
    xypsms_r_f = (RmgType *)TradeImages::xypsms_r[0];
    xypsps_r_f = (RmgType *)TradeImages::xypsps_r[0];
    xymsms_r_f = (RmgType *)TradeImages::xymsms_r[0];
    xymsps_r_f = (RmgType *)TradeImages::xymsps_r[0];
    xypsms_s_f = (RmgType *)TradeImages::xypsms_s[0];
    xypsps_s_f = (RmgType *)TradeImages::xypsps_s[0];
    xymsms_s_f = (RmgType *)TradeImages::xymsms_s[0];
    xymsps_s_f = (RmgType *)TradeImages::xymsps_s[0];
    m0_s_f = (RmgType *)TradeImages::m0_s[0];
    m0_r_f = (RmgType *)TradeImages::m0_r[0];

    if(images > MAX_TRADE_IMAGES) {
       rmg_error_handler (__FILE__, __LINE__, "Images count too high in trade_imagesx_async. Modify and recompile may be required.\n");
    }

    corner_node_stride = ACTIVE_THREADS * MAX_IMG3; 
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
      TradeImages::RMG_MPI_trade( &m0_r_f[0*corner_node_stride], ACTIVE_THREADS * img3, RMG_MPI_IRECV, -1, -1, -1, grid_comm, 18, istate, &TradeImages::rreqs[18]);
      TradeImages::RMG_MPI_trade( &m0_r_f[1*corner_node_stride], ACTIVE_THREADS * img3, RMG_MPI_IRECV, -1, -1, 1, grid_comm, 19, istate, &TradeImages::rreqs[19]);
      TradeImages::RMG_MPI_trade( &m0_r_f[2*corner_node_stride], ACTIVE_THREADS * img3, RMG_MPI_IRECV, -1, 1, -1, grid_comm, 20, istate, &TradeImages::rreqs[20]);
      TradeImages::RMG_MPI_trade( &m0_r_f[3*corner_node_stride], ACTIVE_THREADS * img3, RMG_MPI_IRECV, -1, 1, 1, grid_comm, 21, istate, &TradeImages::rreqs[21]);
      TradeImages::RMG_MPI_trade( &m0_r_f[4*corner_node_stride], ACTIVE_THREADS * img3, RMG_MPI_IRECV, 1, -1, -1, grid_comm, 22, istate, &TradeImages::rreqs[22]);
      TradeImages::RMG_MPI_trade( &m0_r_f[5*corner_node_stride], ACTIVE_THREADS * img3, RMG_MPI_IRECV, 1, -1, 1, grid_comm, 23, istate, &TradeImages::rreqs[23]);
      TradeImages::RMG_MPI_trade( &m0_r_f[6*corner_node_stride], ACTIVE_THREADS * img3, RMG_MPI_IRECV, 1, 1, -1, grid_comm, 24, istate, &TradeImages::rreqs[24]);
      TradeImages::RMG_MPI_trade( &m0_r_f[7*corner_node_stride], ACTIVE_THREADS * img3, RMG_MPI_IRECV, 1, 1, 1, grid_comm, 25, istate, &TradeImages::rreqs[25]);

        // The z planes
      TradeImages::RMG_MPI_trade(frdz2n_f, ACTIVE_THREADS * zlen, RMG_MPI_IRECV, 0, 0, 1, grid_comm, 0, istate, &TradeImages::rreqs[0]);
      TradeImages::RMG_MPI_trade(frdz1n_f, ACTIVE_THREADS * zlen, RMG_MPI_IRECV, 0, 0, -1, grid_comm, 1, istate, &TradeImages::rreqs[1]);

        // The y planes
      TradeImages::RMG_MPI_trade(frdy2n_f, ACTIVE_THREADS * ylen, RMG_MPI_IRECV, 0, 1, 0, grid_comm, 2, istate, &TradeImages::rreqs[2]);
      TradeImages::RMG_MPI_trade(frdy1n_f, ACTIVE_THREADS * ylen, RMG_MPI_IRECV, 0, -1, 0, grid_comm, 3, istate, &TradeImages::rreqs[3]);

        // The x planes
      TradeImages::RMG_MPI_trade(frdx2n_f, ACTIVE_THREADS * xlen, RMG_MPI_IRECV, 1, 0, 0, grid_comm, 4, istate, &TradeImages::rreqs[4]);
      TradeImages::RMG_MPI_trade(frdx1n_f, ACTIVE_THREADS * xlen, RMG_MPI_IRECV, -1, 0, 0, grid_comm, 5, istate, &TradeImages::rreqs[5]);

        // And the yz-plane edges
      TradeImages::RMG_MPI_trade(yzpsms_r_f, ACTIVE_THREADS * yzlen, RMG_MPI_IRECV, 0, 1, -1, grid_comm, 6, istate, &TradeImages::rreqs[6]);
      TradeImages::RMG_MPI_trade(yzpsps_r_f, ACTIVE_THREADS * yzlen, RMG_MPI_IRECV, 0, 1, 1, grid_comm, 7, istate, &TradeImages::rreqs[7]);
      TradeImages::RMG_MPI_trade(yzmsms_r_f, ACTIVE_THREADS * yzlen, RMG_MPI_IRECV, 0, -1, -1, grid_comm, 8, istate, &TradeImages::rreqs[8]);
      TradeImages::RMG_MPI_trade(yzmsps_r_f, ACTIVE_THREADS * yzlen, RMG_MPI_IRECV, 0, -1, 1, grid_comm, 9, istate, &TradeImages::rreqs[9]);

        // And the xy-plane edges
      TradeImages::RMG_MPI_trade(xypsms_r_f, ACTIVE_THREADS * xylen, RMG_MPI_IRECV, 1, -1, 0, grid_comm, 10, istate, &TradeImages::rreqs[10]);
      TradeImages::RMG_MPI_trade(xypsps_r_f, ACTIVE_THREADS * xylen, RMG_MPI_IRECV, 1, 1, 0, grid_comm, 11, istate, &TradeImages::rreqs[11]);
      TradeImages::RMG_MPI_trade(xymsms_r_f, ACTIVE_THREADS * xylen, RMG_MPI_IRECV, -1, -1, 0, grid_comm, 12, istate, &TradeImages::rreqs[12]);
      TradeImages::RMG_MPI_trade(xymsps_r_f, ACTIVE_THREADS * xylen, RMG_MPI_IRECV, -1, 1, 0, grid_comm, 13, istate, &TradeImages::rreqs[13]);

        // And the xz-plane edges
      TradeImages::RMG_MPI_trade(xzpsms_r_f, ACTIVE_THREADS * xzlen, RMG_MPI_IRECV, 1, 0, -1, grid_comm, 14, istate, &TradeImages::rreqs[14]);
      TradeImages::RMG_MPI_trade(xzpsps_r_f, ACTIVE_THREADS * xzlen, RMG_MPI_IRECV, 1, 0, 1, grid_comm, 15, istate, &TradeImages::rreqs[15]);
      TradeImages::RMG_MPI_trade(xzmsms_r_f, ACTIVE_THREADS * xzlen, RMG_MPI_IRECV, -1, 0, -1, grid_comm, 16, istate, &TradeImages::rreqs[16]);
      TradeImages::RMG_MPI_trade(xzmsps_r_f, ACTIVE_THREADS * xzlen, RMG_MPI_IRECV, -1, 0, 1, grid_comm, 17, istate, &TradeImages::rreqs[17]);

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
    T->thread_barrier_wait(true);
    if(tid == 0) {

        // Corners
      TradeImages::RMG_MPI_trade( &m0_s_f[0*corner_node_stride], ACTIVE_THREADS * img3, RMG_MPI_ISEND, 1, 1, 1, grid_comm, 18, istate, &TradeImages::sreqs[18]);
      TradeImages::RMG_MPI_trade( &m0_s_f[1*corner_node_stride], ACTIVE_THREADS * img3, RMG_MPI_ISEND, 1, 1, -1, grid_comm, 19, istate, &TradeImages::sreqs[19]);
      TradeImages::RMG_MPI_trade( &m0_s_f[2*corner_node_stride], ACTIVE_THREADS * img3, RMG_MPI_ISEND, 1, -1, 1, grid_comm, 20, istate, &TradeImages::sreqs[20]);
      TradeImages::RMG_MPI_trade( &m0_s_f[3*corner_node_stride], ACTIVE_THREADS * img3, RMG_MPI_ISEND, 1, -1, -1, grid_comm, 21, istate, &TradeImages::sreqs[21]);
      TradeImages::RMG_MPI_trade( &m0_s_f[4*corner_node_stride], ACTIVE_THREADS * img3, RMG_MPI_ISEND, -1, 1, 1, grid_comm, 22, istate, &TradeImages::sreqs[22]);
      TradeImages::RMG_MPI_trade( &m0_s_f[5*corner_node_stride], ACTIVE_THREADS * img3, RMG_MPI_ISEND, -1, 1, -1, grid_comm, 23, istate, &TradeImages::sreqs[23]);
      TradeImages::RMG_MPI_trade( &m0_s_f[6*corner_node_stride], ACTIVE_THREADS * img3, RMG_MPI_ISEND, -1, -1, 1, grid_comm, 24, istate, &TradeImages::sreqs[24]);
      TradeImages::RMG_MPI_trade( &m0_s_f[7*corner_node_stride], ACTIVE_THREADS * img3, RMG_MPI_ISEND, -1, -1, -1, grid_comm, 25, istate, &TradeImages::sreqs[25]);

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


    T->thread_barrier_wait(true);
    if(tid == 0) {

        // Send z planes
      TradeImages::RMG_MPI_trade(frdz1_f, ACTIVE_THREADS * zlen, RMG_MPI_ISEND, 0, 0, -1, grid_comm, 0, istate, &TradeImages::sreqs[0]);
      TradeImages::RMG_MPI_trade(frdz2_f, ACTIVE_THREADS * zlen, RMG_MPI_ISEND, 0, 0, 1, grid_comm, 1, istate, &TradeImages::sreqs[1]);

        // Send the north and south planes
      TradeImages::RMG_MPI_trade(frdy1_f, ACTIVE_THREADS * ylen, RMG_MPI_ISEND, 0, -1, 0, grid_comm, 2, istate, &TradeImages::sreqs[2]);
      TradeImages::RMG_MPI_trade(frdy2_f, ACTIVE_THREADS * ylen, RMG_MPI_ISEND, 0, 1, 0, grid_comm, 3, istate, &TradeImages::sreqs[3]);

        // Send the east and west planes
      TradeImages::RMG_MPI_trade(frdx1_f, ACTIVE_THREADS * xlen, RMG_MPI_ISEND, -1, 0, 0, grid_comm, 4, istate, &TradeImages::sreqs[4]);
      TradeImages::RMG_MPI_trade(frdx2_f, ACTIVE_THREADS * xlen, RMG_MPI_ISEND, 1, 0, 0, grid_comm, 5, istate, &TradeImages::sreqs[5]);

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
    T->thread_barrier_wait(true);
    if(tid == 0) {

        // Send yz-plane edges
      TradeImages::RMG_MPI_trade(yzpsms_s_f, ACTIVE_THREADS * yzlen, RMG_MPI_ISEND, 0, -1, 1, grid_comm, 6, istate, &TradeImages::sreqs[6]);
      TradeImages::RMG_MPI_trade(yzpsps_s_f, ACTIVE_THREADS * yzlen, RMG_MPI_ISEND, 0, -1, -1, grid_comm, 7, istate, &TradeImages::sreqs[7]);
      TradeImages::RMG_MPI_trade(yzmsms_s_f, ACTIVE_THREADS * yzlen, RMG_MPI_ISEND, 0, 1, 1, grid_comm, 8, istate, &TradeImages::sreqs[8]);
      TradeImages::RMG_MPI_trade(yzmsps_s_f, ACTIVE_THREADS * yzlen, RMG_MPI_ISEND, 0, 1, -1, grid_comm, 9, istate, &TradeImages::sreqs[9]);

        // Send xy-plane edges
      TradeImages::RMG_MPI_trade(xypsms_s_f, ACTIVE_THREADS * xylen, RMG_MPI_ISEND, -1, 1, 0, grid_comm, 10, istate, &TradeImages::sreqs[10]);
      TradeImages::RMG_MPI_trade(xypsps_s_f, ACTIVE_THREADS * xylen, RMG_MPI_ISEND, -1, -1, 0, grid_comm, 11, istate, &TradeImages::sreqs[11]);
      TradeImages::RMG_MPI_trade(xymsms_s_f, ACTIVE_THREADS * xylen, RMG_MPI_ISEND, 1, 1, 0, grid_comm, 12, istate, &TradeImages::sreqs[12]);
      TradeImages::RMG_MPI_trade(xymsps_s_f, ACTIVE_THREADS * xylen, RMG_MPI_ISEND, 1, -1, 0, grid_comm, 13, istate, &TradeImages::sreqs[13]);

        // Send xz-plane edges
      TradeImages::RMG_MPI_trade(xzpsms_s_f, ACTIVE_THREADS * xzlen, RMG_MPI_ISEND, -1, 0, 1, grid_comm, 14, istate, &TradeImages::sreqs[14]);
      TradeImages::RMG_MPI_trade(xzpsps_s_f, ACTIVE_THREADS * xzlen, RMG_MPI_ISEND, -1, 0, -1, grid_comm, 15, istate, &TradeImages::sreqs[15]);
      TradeImages::RMG_MPI_trade(xzmsms_s_f, ACTIVE_THREADS * xzlen, RMG_MPI_ISEND, 1, 0, 1, grid_comm, 16, istate, &TradeImages::sreqs[16]);
      TradeImages::RMG_MPI_trade(xzmsps_s_f, ACTIVE_THREADS * xzlen, RMG_MPI_ISEND, 1, 0, -1, grid_comm, 17, istate, &TradeImages::sreqs[17]);

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

            for(int idx = 0;idx < dimz;idx++)
                w[iys2 + images + idx] = f[iys + idx];


        }                       /* end for */

    }                           /* end for */


    // Wait for all the recvs to finish
    if(tid == 0) {
        retval = MPI_Waitall(26, TradeImages::rreqs, MPI_STATUSES_IGNORE);
        if(retval != MPI_SUCCESS) rmg_error_handler (__FILE__, __LINE__, "Error in MPI_Waitall.\n");
    }
    T->thread_barrier_wait(true);


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
        retval = MPI_Waitall(26, TradeImages::sreqs, MPI_STATUSES_IGNORE);
        if(retval != MPI_SUCCESS) rmg_error_handler (__FILE__, __LINE__, "Error in MPI_Waitall.\n");
    }
    T->thread_barrier_wait(true);

} // end trade_imagesx_async



// Asynchronous image trades for central finite difference operators
template <typename RmgType>
void TradeImages::trade_imagesx_central_async (RmgType * f, RmgType * w, int dimx, int dimy, int dimz, int images)
{
    BaseThread *T = BaseThread::getBaseThread(0);
    int tid = T->get_thread_tid();
    if(tid < 0) tid = 0;
    int ACTIVE_THREADS = 1;
    if(T->is_loop_over_states()) ACTIVE_THREADS = T->barrier->barrier_count();

    int istate = T->get_thread_basetag();
    MPI_Comm grid_comm = T->get_unique_comm(istate);

    int ix, iy, iz, incx, incy, incx0, incy0, index, tim;
    int ixs, iys, ixs2, iys2, c1, idx;
    int xlen, ylen, zlen;
    int retval;

    RmgType *frdx1_f, *frdx2_f, *frdy1_f, *frdy2_f, *frdz1_f, *frdz2_f;
    RmgType *frdx1n_f, *frdx2n_f, *frdy1n_f, *frdy2n_f, *frdz1n_f, *frdz2n_f;

    frdx1_f = (RmgType *)TradeImages::frdx1[0];
    frdx2_f = (RmgType *)TradeImages::frdx2[0];
    frdy1_f = (RmgType *)TradeImages::frdy1[0];
    frdy2_f = (RmgType *)TradeImages::frdy2[0];
    frdz1_f = (RmgType *)TradeImages::frdz1[0];
    frdz2_f = (RmgType *)TradeImages::frdz2[0];
    frdx1n_f = (RmgType *)TradeImages::frdx1n[0];
    frdx2n_f = (RmgType *)TradeImages::frdx2n[0];
    frdy1n_f = (RmgType *)TradeImages::frdy1n[0];
    frdy2n_f = (RmgType *)TradeImages::frdy2n[0];
    frdz1n_f = (RmgType *)TradeImages::frdz1n[0];
    frdz2n_f = (RmgType *)TradeImages::frdz2n[0];


    if(images > MAX_TRADE_IMAGES) {
       rmg_error_handler (__FILE__, __LINE__, "Images count too high in trade_imagesx_async. Modify and recompile may be required.\n");
    }


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
      TradeImages::RMG_MPI_trade(frdz2n_f, ACTIVE_THREADS * zlen, RMG_MPI_IRECV, 0, 0, 1, grid_comm, 0, istate, &TradeImages::rreqs[0]);
      TradeImages::RMG_MPI_trade(frdz1n_f, ACTIVE_THREADS * zlen, RMG_MPI_IRECV, 0, 0, -1, grid_comm, 1, istate, &TradeImages::rreqs[1]);

        // The y planes
      TradeImages::RMG_MPI_trade(frdy2n_f, ACTIVE_THREADS * ylen, RMG_MPI_IRECV, 0, 1, 0, grid_comm, 2, istate, &TradeImages::rreqs[2]);
      TradeImages::RMG_MPI_trade(frdy1n_f, ACTIVE_THREADS * ylen, RMG_MPI_IRECV, 0, -1, 0, grid_comm, 3, istate, &TradeImages::rreqs[3]);

        // The x planes
      TradeImages::RMG_MPI_trade(frdx2n_f, ACTIVE_THREADS * xlen, RMG_MPI_IRECV, 1, 0, 0, grid_comm, 4, istate, &TradeImages::rreqs[4]);
      TradeImages::RMG_MPI_trade(frdx1n_f, ACTIVE_THREADS * xlen, RMG_MPI_IRECV, -1, 0, 0, grid_comm, 5, istate, &TradeImages::rreqs[5]);

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


    T->thread_barrier_wait(true);
    if(tid == 0) {

        // Send z planes
      TradeImages::RMG_MPI_trade(frdz1_f, ACTIVE_THREADS * zlen, RMG_MPI_ISEND, 0, 0, -1, grid_comm, 0, istate, &TradeImages::sreqs[0]);
      TradeImages::RMG_MPI_trade(frdz2_f, ACTIVE_THREADS * zlen, RMG_MPI_ISEND, 0, 0, 1, grid_comm, 1, istate, &TradeImages::sreqs[1]);

        // Send the north and south planes
      TradeImages::RMG_MPI_trade(frdy1_f, ACTIVE_THREADS * ylen, RMG_MPI_ISEND, 0, -1, 0, grid_comm, 2, istate, &TradeImages::sreqs[2]);
      TradeImages::RMG_MPI_trade(frdy2_f, ACTIVE_THREADS * ylen, RMG_MPI_ISEND, 0, 1, 0, grid_comm, 3, istate, &TradeImages::sreqs[3]);

        // Send the east and west planes
      TradeImages::RMG_MPI_trade(frdx1_f, ACTIVE_THREADS * xlen, RMG_MPI_ISEND, -1, 0, 0, grid_comm, 4, istate, &TradeImages::sreqs[4]);
      TradeImages::RMG_MPI_trade(frdx2_f, ACTIVE_THREADS * xlen, RMG_MPI_ISEND, 1, 0, 0, grid_comm, 5, istate, &TradeImages::sreqs[5]);

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

            for(int idx = 0;idx < dimz;idx++)
                w[iys2 + images + idx] = f[iys + idx];


        }                       /* end for */

    }                           /* end for */


    // Wait for all the recvs to finish
    if(tid == 0) {
        retval = MPI_Waitall(6, TradeImages::rreqs, MPI_STATUSES_IGNORE);
        if(retval != MPI_SUCCESS) rmg_error_handler (__FILE__, __LINE__, "Error in MPI_Waitall.\n");
    }
    T->thread_barrier_wait(true);


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
        retval = MPI_Waitall(6, TradeImages::sreqs, MPI_STATUSES_IGNORE);
        if(retval != MPI_SUCCESS) rmg_error_handler (__FILE__, __LINE__, "Error in MPI_Waitall.\n");
    }
    T->thread_barrier_wait(true);

} // end trade_imagesx_central_async


// Asynchronous image trades for central finite difference operators using the queue manager
template <typename RmgType>
void TradeImages::trade_imagesx_central_async_managed (RmgType * f, RmgType * w, int dimx, int dimy, int dimz, int images)
{

    if(images > MAX_TRADE_IMAGES) {
       rmg_error_handler (__FILE__, __LINE__, "Images count too high in trade_imagesx_async. Modify and recompile may be required.\n");
    }

    BaseThread *T = BaseThread::getBaseThread(0);
    int tid = T->get_thread_tid();
    if(tid < 0) tid = 0;
    int istate = T->get_thread_basetag();
    MPI_Comm grid_comm = T->get_unique_comm(istate);

    int ix, iy, iz, index;
    int ixs, iys, ixs2, iys2, c1, idx;

    std::atomic_bool is_completed_r[6];
    std::atomic_bool is_completed_s[6];
    std::mutex send_mutex;
    std::condition_variable send_cv;

    mpi_queue_item_t qitems_r[6];
    mpi_queue_item_t qitems_s[6];

    RmgType *frdx1_f, *frdx2_f, *frdy1_f, *frdy2_f, *frdz1_f, *frdz2_f;
    RmgType *frdx1n_f, *frdx2n_f, *frdy1n_f, *frdy2n_f, *frdz1n_f, *frdz2n_f;

    int tim = 2 * images;
    int incx = (dimy + tim) * (dimz + tim);
    int incy = dimz + tim;
    int incx0 = dimy * dimz;
    int incy0 = dimz;
    int zlen = dimx * dimy * images;
    int ylen = dimx * dimz * images;
    int xlen = dimy * dimz * images;

    for(int it = 0;it < 6;it++)
    {
        qitems_r[it].is_completed = &is_completed_r[it];
        qitems_s[it].is_completed = &is_completed_s[it];
    }


    frdx1_f = (RmgType *)TradeImages::frdx1[tid];
    frdx2_f = (RmgType *)TradeImages::frdx2[tid];
    frdy1_f = (RmgType *)TradeImages::frdy1[tid];
    frdy2_f = (RmgType *)TradeImages::frdy2[tid];
    frdz1_f = (RmgType *)TradeImages::frdz1[tid];
    frdz2_f = (RmgType *)TradeImages::frdz2[tid];
    frdx1n_f = (RmgType *)TradeImages::frdx1n[tid];
    frdx2n_f = (RmgType *)TradeImages::frdx2n[tid];
    frdy1n_f = (RmgType *)TradeImages::frdy1n[tid];
    frdy2n_f = (RmgType *)TradeImages::frdy2n[tid];
    frdz1n_f = (RmgType *)TradeImages::frdz1n[tid];
    frdz2n_f = (RmgType *)TradeImages::frdz2n[tid];



    // The z planes
    TradeImages::RMG_MPI_queue_trade(frdz2n_f, zlen, RMG_MPI_IRECV, 0, 0, 1, grid_comm, 0, istate, qitems_r[0]);
    TradeImages::RMG_MPI_queue_trade(frdz1n_f, zlen, RMG_MPI_IRECV, 0, 0, -1, grid_comm, 1, istate, qitems_r[1]);

    // The y planes
    TradeImages::RMG_MPI_queue_trade(frdy2n_f, ylen, RMG_MPI_IRECV, 0, 1, 0, grid_comm, 2, istate, qitems_r[2]);
    TradeImages::RMG_MPI_queue_trade(frdy1n_f, ylen, RMG_MPI_IRECV, 0, -1, 0, grid_comm, 3, istate, qitems_r[3]);

    // The x planes
    TradeImages::RMG_MPI_queue_trade(frdx2n_f, xlen, RMG_MPI_IRECV, 1, 0, 0, grid_comm, 4, istate, qitems_r[4]);
    TradeImages::RMG_MPI_queue_trade(frdx1n_f, xlen, RMG_MPI_IRECV, -1, 0, 0, grid_comm, 5, istate, qitems_r[5]);


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

                frdz1_f[idx] = f[index];
                frdz2_f[idx] = f[index+c1];
                idx++;

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */


    // Send z planes
    TradeImages::RMG_MPI_queue_trade(frdz1_f, zlen, RMG_MPI_ISEND, 0, 0, -1, grid_comm, 0, istate, qitems_s[0]);
    TradeImages::RMG_MPI_queue_trade(frdz2_f, zlen, RMG_MPI_ISEND, 0, 0, 1, grid_comm, 1, istate, qitems_s[1]);


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

                frdy1_f[idx] = f[index];
                frdy2_f[idx] = f[index + c1];
                idx++;

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */


    // Send the north and south planes
    TradeImages::RMG_MPI_queue_trade(frdy1_f, ylen, RMG_MPI_ISEND, 0, -1, 0, grid_comm, 2, istate, qitems_s[2]);
    TradeImages::RMG_MPI_queue_trade(frdy2_f, ylen, RMG_MPI_ISEND, 0, 1, 0, grid_comm, 3, istate, qitems_s[3]);


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

                frdx1_f[idx] = f[index];
                frdx2_f[idx] = f[index + c1];
                idx++;

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */


    // Send the east and west planes
    TradeImages::RMG_MPI_queue_trade(frdx1_f, xlen, RMG_MPI_ISEND, -1, 0, 0, grid_comm, 4, istate, qitems_s[4]);
    TradeImages::RMG_MPI_queue_trade(frdx2_f, xlen, RMG_MPI_ISEND, 1, 0, 0, grid_comm, 5, istate, qitems_s[5]);


    // Load up w with the interior points
    for (ix = 0; ix < dimx; ix++)
    {
        ixs = ix * incx0;
        ixs2 = (ix + images) * incx;

        for (iy = 0; iy < dimy; iy++)
        {
            iys = ixs + iy * incy0;
            iys2 = ixs2 + (iy + images) * incy;

            for(int idx = 0;idx < dimz;idx++)
                w[iys2 + images + idx] = f[iys + idx];


        }                       /* end for */

    }                           /* end for */


    // Loop and unpack items as they become available
    int items_completed = 0;
    while(items_completed < 6)
    {

        if(!qitems_r[0].is_unpacked)
        {
            if(is_completed_r[0].load(std::memory_order_acquire))
            {
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
                            w[index + c1] = frdz2n_f[idx];
                            idx++;
                        }                   /* end for */
                    }                       /* end for */
                }                           /* end for */
                items_completed++;
                qitems_r[0].is_unpacked = true;
            }
        }

        if(!qitems_r[1].is_unpacked)
        {
            if(is_completed_r[1].load(std::memory_order_acquire))
            {
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
                            w[index] = frdz1n_f[idx];
                            idx++;
                        }                   /* end for */
                    }                       /* end for */
                }                           /* end for */
                items_completed++;
                qitems_r[1].is_unpacked = true;
            }
        }

// ***********************************************************

        if(!qitems_r[2].is_unpacked)
        { 
            if(is_completed_r[2].load(std::memory_order_acquire))
            {
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
                            w[index + c1] = frdy2n_f[idx];
                            idx++;
                        }                   /* end for */
                    }                       /* end for */
                }                           /* end for */
                items_completed++;
                qitems_r[2].is_unpacked = true;
            }
        }

        if(!qitems_r[3].is_unpacked)
        {
            if(is_completed_r[3].load(std::memory_order_acquire))
            {
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
                            w[index] = frdy1n_f[idx];
                            idx++;
                        }                   /* end for */
                    }                       /* end for */
                }                           /* end for */
                items_completed++;
                qitems_r[3].is_unpacked = true;
            }
        }

// ***********************************************************

        if(!qitems_r[4].is_unpacked)
        {
            if(is_completed_r[4].load(std::memory_order_acquire))
            {
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
                            w[index + c1] = frdx2n_f[idx];
                            idx++;
                        }                   /* end for */
                    }                       /* end for */
                }                           /* end for */
                items_completed++;
                qitems_r[4].is_unpacked = true;
            }
        }
 
        if(!qitems_r[5].is_unpacked)
        {
            if(is_completed_r[5].load(std::memory_order_acquire))
            {
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
                            w[index] = frdx1n_f[idx];
                            idx++;
                        }                   /* end for */
                    }                       /* end for */
                }
                items_completed++;
                qitems_r[5].is_unpacked = true;
            }
        } 

    } // end while

   //this->queue->spinwaitall(is_completed_r, 6);
   this->queue->spinwaitall(is_completed_s, 6);

} // end trade_imagesx_central_async_managed


// Used by coarse level multigrid routines which assume that the central data is already present
// and that the f array is sized [dimx+2][dimy+2][dimz+2]
//
template <typename RmgType>
void TradeImages::trade_images1_async (RmgType * f, int dimx, int dimy, int dimz)
{
    BaseThread *T = BaseThread::getBaseThread(0);
    int tid = T->get_thread_tid();
    if(tid < 0) tid = 0;
    int ACTIVE_THREADS = 1;
    if(T->is_loop_over_states()) ACTIVE_THREADS = T->barrier->barrier_count();

    int istate = T->get_thread_basetag();
    MPI_Comm grid_comm = T->get_unique_comm(istate);

    int ix, iy, iz, incx, incy, index;
    int ixs2, iys2, c1, c2, c3, idx, idx1;
    int xlen, ylen, zlen, yzlen, xylen, xzlen;
    int corner_node_stride, node_idx, retval;

    RmgType *frdx1_f, *frdx2_f, *frdy1_f, *frdy2_f, *frdz1_f, *frdz2_f;
    RmgType *frdx1n_f, *frdx2n_f, *frdy1n_f, *frdy2n_f, *frdz1n_f, *frdz2n_f;
    RmgType *yzpsms_r_f, *yzpsps_r_f, *yzmsms_r_f, *yzmsps_r_f;
    RmgType *yzpsms_s_f, *yzpsps_s_f, *yzmsms_s_f, *yzmsps_s_f;
    RmgType *xzpsms_r_f, *xzpsps_r_f, *xzmsms_r_f, *xzmsps_r_f;
    RmgType *xzpsms_s_f, *xzpsps_s_f, *xzmsms_s_f, *xzmsps_s_f;
    RmgType *xypsms_r_f, *xypsps_r_f, *xymsms_r_f, *xymsps_r_f;
    RmgType *xypsms_s_f, *xypsps_s_f, *xymsms_s_f, *xymsps_s_f;
    RmgType *m0_s_f, *m0_r_f;

    frdx1_f = (RmgType *)TradeImages::frdx1[0];
    frdx2_f = (RmgType *)TradeImages::frdx2[0];
    frdy1_f = (RmgType *)TradeImages::frdy1[0];
    frdy2_f = (RmgType *)TradeImages::frdy2[0];
    frdz1_f = (RmgType *)TradeImages::frdz1[0];
    frdz2_f = (RmgType *)TradeImages::frdz2[0];
    frdx1n_f = (RmgType *)TradeImages::frdx1n[0];
    frdx2n_f = (RmgType *)TradeImages::frdx2n[0];
    frdy1n_f = (RmgType *)TradeImages::frdy1n[0];
    frdy2n_f = (RmgType *)TradeImages::frdy2n[0];
    frdz1n_f = (RmgType *)TradeImages::frdz1n[0];
    frdz2n_f = (RmgType *)TradeImages::frdz2n[0];
    yzpsms_r_f = (RmgType *)TradeImages::yzpsms_r[0];
    yzpsps_r_f = (RmgType *)TradeImages::yzpsps_r[0];
    yzmsms_r_f = (RmgType *)TradeImages::yzmsms_r[0];
    yzmsps_r_f = (RmgType *)TradeImages::yzmsps_r[0];
    yzpsms_s_f = (RmgType *)TradeImages::yzpsms_s[0];
    yzpsps_s_f = (RmgType *)TradeImages::yzpsps_s[0];
    yzmsms_s_f = (RmgType *)TradeImages::yzmsms_s[0];
    yzmsps_s_f = (RmgType *)TradeImages::yzmsps_s[0];
    xzpsms_r_f = (RmgType *)TradeImages::xzpsms_r[0];
    xzpsps_r_f = (RmgType *)TradeImages::xzpsps_r[0];
    xzmsms_r_f = (RmgType *)TradeImages::xzmsms_r[0];
    xzmsps_r_f = (RmgType *)TradeImages::xzmsps_r[0];
    xzpsms_s_f = (RmgType *)TradeImages::xzpsms_s[0];
    xzpsps_s_f = (RmgType *)TradeImages::xzpsps_s[0];
    xzmsms_s_f = (RmgType *)TradeImages::xzmsms_s[0];
    xzmsps_s_f = (RmgType *)TradeImages::xzmsps_s[0];
    xypsms_r_f = (RmgType *)TradeImages::xypsms_r[0];
    xypsps_r_f = (RmgType *)TradeImages::xypsps_r[0];
    xymsms_r_f = (RmgType *)TradeImages::xymsms_r[0];
    xymsps_r_f = (RmgType *)TradeImages::xymsps_r[0];
    xypsms_s_f = (RmgType *)TradeImages::xypsms_s[0];
    xypsps_s_f = (RmgType *)TradeImages::xypsps_s[0];
    xymsms_s_f = (RmgType *)TradeImages::xymsms_s[0];
    xymsps_s_f = (RmgType *)TradeImages::xymsps_s[0];
    m0_s_f = (RmgType *)TradeImages::m0_s[0];
    m0_r_f = (RmgType *)TradeImages::m0_r[0];

    corner_node_stride = ACTIVE_THREADS * MAX_IMG3;

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
      TradeImages::RMG_MPI_trade( &m0_r_f[0*corner_node_stride], ACTIVE_THREADS, RMG_MPI_IRECV, -1, -1, -1, grid_comm, 18, istate, &TradeImages::rreqs[18]);
      TradeImages::RMG_MPI_trade( &m0_r_f[1*corner_node_stride], ACTIVE_THREADS, RMG_MPI_IRECV, -1, -1, 1, grid_comm, 19, istate, &TradeImages::rreqs[19]);
      TradeImages::RMG_MPI_trade( &m0_r_f[2*corner_node_stride], ACTIVE_THREADS, RMG_MPI_IRECV, -1, 1, -1, grid_comm, 20, istate, &TradeImages::rreqs[20]);
      TradeImages::RMG_MPI_trade( &m0_r_f[3*corner_node_stride], ACTIVE_THREADS, RMG_MPI_IRECV, -1, 1, 1, grid_comm, 21, istate, &TradeImages::rreqs[21]);
      TradeImages::RMG_MPI_trade( &m0_r_f[4*corner_node_stride], ACTIVE_THREADS, RMG_MPI_IRECV, 1, -1, -1, grid_comm, 22, istate, &TradeImages::rreqs[22]);
      TradeImages::RMG_MPI_trade( &m0_r_f[5*corner_node_stride], ACTIVE_THREADS, RMG_MPI_IRECV, 1, -1, 1, grid_comm, 23, istate, &TradeImages::rreqs[23]);
      TradeImages::RMG_MPI_trade( &m0_r_f[6*corner_node_stride], ACTIVE_THREADS, RMG_MPI_IRECV, 1, 1, -1, grid_comm, 24, istate, &TradeImages::rreqs[24]);
      TradeImages::RMG_MPI_trade( &m0_r_f[7*corner_node_stride], ACTIVE_THREADS, RMG_MPI_IRECV, 1, 1, 1, grid_comm, 25, istate, &TradeImages::rreqs[25]);

        // The z planes
      TradeImages::RMG_MPI_trade(frdz2n_f, ACTIVE_THREADS * zlen, RMG_MPI_IRECV, 0, 0, 1, grid_comm, 0, istate, &TradeImages::rreqs[0]);
      TradeImages::RMG_MPI_trade(frdz1n_f, ACTIVE_THREADS * zlen, RMG_MPI_IRECV, 0, 0, -1, grid_comm, 1, istate, &TradeImages::rreqs[1]);

        // The y planes
      TradeImages::RMG_MPI_trade(frdy2n_f, ACTIVE_THREADS * ylen, RMG_MPI_IRECV, 0, 1, 0, grid_comm, 2, istate, &TradeImages::rreqs[2]);
      TradeImages::RMG_MPI_trade(frdy1n_f, ACTIVE_THREADS * ylen, RMG_MPI_IRECV, 0, -1, 0, grid_comm, 3, istate, &TradeImages::rreqs[3]);

        // The x planes
      TradeImages::RMG_MPI_trade(frdx2n_f, ACTIVE_THREADS * xlen, RMG_MPI_IRECV, 1, 0, 0, grid_comm, 4, istate, &TradeImages::rreqs[4]);
      TradeImages::RMG_MPI_trade(frdx1n_f, ACTIVE_THREADS * xlen, RMG_MPI_IRECV, -1, 0, 0, grid_comm, 5, istate, &TradeImages::rreqs[5]);

        // And the yz-plane edges
      TradeImages::RMG_MPI_trade(yzpsms_r_f, ACTIVE_THREADS * yzlen, RMG_MPI_IRECV, 0, 1, -1, grid_comm, 6, istate, &TradeImages::rreqs[6]);
      TradeImages::RMG_MPI_trade(yzpsps_r_f, ACTIVE_THREADS * yzlen, RMG_MPI_IRECV, 0, 1, 1, grid_comm, 7, istate, &TradeImages::rreqs[7]);
      TradeImages::RMG_MPI_trade(yzmsms_r_f, ACTIVE_THREADS * yzlen, RMG_MPI_IRECV, 0, -1, -1, grid_comm, 8, istate, &TradeImages::rreqs[8]);
      TradeImages::RMG_MPI_trade(yzmsps_r_f, ACTIVE_THREADS * yzlen, RMG_MPI_IRECV, 0, -1, 1, grid_comm, 9, istate, &TradeImages::rreqs[9]);

        // And the xy-plane edges
      TradeImages::RMG_MPI_trade(xypsms_r_f, ACTIVE_THREADS * xylen, RMG_MPI_IRECV, 1, -1, 0, grid_comm, 10, istate, &TradeImages::rreqs[10]);
      TradeImages::RMG_MPI_trade(xypsps_r_f, ACTIVE_THREADS * xylen, RMG_MPI_IRECV, 1, 1, 0, grid_comm, 11, istate, &TradeImages::rreqs[11]);
      TradeImages::RMG_MPI_trade(xymsms_r_f, ACTIVE_THREADS * xylen, RMG_MPI_IRECV, -1, -1, 0, grid_comm, 12, istate, &TradeImages::rreqs[12]);
      TradeImages::RMG_MPI_trade(xymsps_r_f, ACTIVE_THREADS * xylen, RMG_MPI_IRECV, -1, 1, 0, grid_comm, 13, istate, &TradeImages::rreqs[13]);

        // And the xz-plane edges
      TradeImages::RMG_MPI_trade(xzpsms_r_f, ACTIVE_THREADS * xzlen, RMG_MPI_IRECV, 1, 0, -1, grid_comm, 14, istate, &TradeImages::rreqs[14]);
      TradeImages::RMG_MPI_trade(xzpsps_r_f, ACTIVE_THREADS * xzlen, RMG_MPI_IRECV, 1, 0, 1, grid_comm, 15, istate, &TradeImages::rreqs[15]);
      TradeImages::RMG_MPI_trade(xzmsms_r_f, ACTIVE_THREADS * xzlen, RMG_MPI_IRECV, -1, 0, -1, grid_comm, 16, istate, &TradeImages::rreqs[16]);
      TradeImages::RMG_MPI_trade(xzmsps_r_f, ACTIVE_THREADS * xzlen, RMG_MPI_IRECV, -1, 0, 1, grid_comm, 17, istate, &TradeImages::rreqs[17]);

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
    T->thread_barrier_wait(true);
    if(tid == 0) {

        // Corners
      TradeImages::RMG_MPI_trade( &m0_s_f[0*corner_node_stride], ACTIVE_THREADS, RMG_MPI_ISEND, 1, 1, 1, grid_comm, 18, istate, &TradeImages::sreqs[18]);
      TradeImages::RMG_MPI_trade( &m0_s_f[1*corner_node_stride], ACTIVE_THREADS, RMG_MPI_ISEND, 1, 1, -1, grid_comm, 19, istate, &TradeImages::sreqs[19]);
      TradeImages::RMG_MPI_trade( &m0_s_f[2*corner_node_stride], ACTIVE_THREADS, RMG_MPI_ISEND, 1, -1, 1, grid_comm, 20, istate, &TradeImages::sreqs[20]);
      TradeImages::RMG_MPI_trade( &m0_s_f[3*corner_node_stride], ACTIVE_THREADS, RMG_MPI_ISEND, 1, -1, -1, grid_comm, 21, istate, &TradeImages::sreqs[21]);
      TradeImages::RMG_MPI_trade( &m0_s_f[4*corner_node_stride], ACTIVE_THREADS, RMG_MPI_ISEND, -1, 1, 1, grid_comm, 22, istate, &TradeImages::sreqs[22]);
      TradeImages::RMG_MPI_trade( &m0_s_f[5*corner_node_stride], ACTIVE_THREADS, RMG_MPI_ISEND, -1, 1, -1, grid_comm, 23, istate, &TradeImages::sreqs[23]);
      TradeImages::RMG_MPI_trade( &m0_s_f[6*corner_node_stride], ACTIVE_THREADS, RMG_MPI_ISEND, -1, -1, 1, grid_comm, 24, istate, &TradeImages::sreqs[24]);
      TradeImages::RMG_MPI_trade( &m0_s_f[7*corner_node_stride], ACTIVE_THREADS, RMG_MPI_ISEND, -1, -1, -1, grid_comm, 25, istate, &TradeImages::sreqs[25]);

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


    T->thread_barrier_wait(true);
    if(tid == 0) {

        // Send z planes
      TradeImages::RMG_MPI_trade(frdz1_f, ACTIVE_THREADS * zlen, RMG_MPI_ISEND, 0, 0, -1, grid_comm, 0, istate, &TradeImages::sreqs[0]);
      TradeImages::RMG_MPI_trade(frdz2_f, ACTIVE_THREADS * zlen, RMG_MPI_ISEND, 0, 0, 1, grid_comm, 1, istate, &TradeImages::sreqs[1]);

        // Send the north and south planes
      TradeImages::RMG_MPI_trade(frdy1_f, ACTIVE_THREADS * ylen, RMG_MPI_ISEND, 0, -1, 0, grid_comm, 2, istate, &TradeImages::sreqs[2]);
      TradeImages::RMG_MPI_trade(frdy2_f, ACTIVE_THREADS * ylen, RMG_MPI_ISEND, 0, 1, 0, grid_comm, 3, istate, &TradeImages::sreqs[3]);

        // Send the east and west planes
      TradeImages::RMG_MPI_trade(frdx1_f, ACTIVE_THREADS * xlen, RMG_MPI_ISEND, -1, 0, 0, grid_comm, 4, istate, &TradeImages::sreqs[4]);
      TradeImages::RMG_MPI_trade(frdx2_f, ACTIVE_THREADS * xlen, RMG_MPI_ISEND, 1, 0, 0, grid_comm, 5, istate, &TradeImages::sreqs[5]);

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
    T->thread_barrier_wait(true);
    if(tid == 0) {

        // Send yz-plane edges
      TradeImages::RMG_MPI_trade(yzpsms_s_f, ACTIVE_THREADS * yzlen, RMG_MPI_ISEND, 0, -1, 1, grid_comm, 6, istate, &TradeImages::sreqs[6]);
      TradeImages::RMG_MPI_trade(yzpsps_s_f, ACTIVE_THREADS * yzlen, RMG_MPI_ISEND, 0, -1, -1, grid_comm, 7, istate, &TradeImages::sreqs[7]);
      TradeImages::RMG_MPI_trade(yzmsms_s_f, ACTIVE_THREADS * yzlen, RMG_MPI_ISEND, 0, 1, 1, grid_comm, 8, istate, &TradeImages::sreqs[8]);
      TradeImages::RMG_MPI_trade(yzmsps_s_f, ACTIVE_THREADS * yzlen, RMG_MPI_ISEND, 0, 1, -1, grid_comm, 9, istate, &TradeImages::sreqs[9]);

        // Send xy-plane edges
      TradeImages::RMG_MPI_trade(xypsms_s_f, ACTIVE_THREADS * xylen, RMG_MPI_ISEND, -1, 1, 0, grid_comm, 10, istate, &TradeImages::sreqs[10]);
      TradeImages::RMG_MPI_trade(xypsps_s_f, ACTIVE_THREADS * xylen, RMG_MPI_ISEND, -1, -1, 0, grid_comm, 11, istate, &TradeImages::sreqs[11]);
      TradeImages::RMG_MPI_trade(xymsms_s_f, ACTIVE_THREADS * xylen, RMG_MPI_ISEND, 1, 1, 0, grid_comm, 12, istate, &TradeImages::sreqs[12]);
      TradeImages::RMG_MPI_trade(xymsps_s_f, ACTIVE_THREADS * xylen, RMG_MPI_ISEND, 1, -1, 0, grid_comm, 13, istate, &TradeImages::sreqs[13]);

        // Send xz-plane edges
      TradeImages::RMG_MPI_trade(xzpsms_s_f, ACTIVE_THREADS * xzlen, RMG_MPI_ISEND, -1, 0, 1, grid_comm, 14, istate, &TradeImages::sreqs[14]);
      TradeImages::RMG_MPI_trade(xzpsps_s_f, ACTIVE_THREADS * xzlen, RMG_MPI_ISEND, -1, 0, -1, grid_comm, 15, istate, &TradeImages::sreqs[15]);
      TradeImages::RMG_MPI_trade(xzmsms_s_f, ACTIVE_THREADS * xzlen, RMG_MPI_ISEND, 1, 0, 1, grid_comm, 16, istate, &TradeImages::sreqs[16]);
      TradeImages::RMG_MPI_trade(xzmsps_s_f, ACTIVE_THREADS * xzlen, RMG_MPI_ISEND, 1, 0, -1, grid_comm, 17, istate, &TradeImages::sreqs[17]);

    }

    // Wait for all the recvs to finish
    if(tid == 0) {
        retval = MPI_Waitall(26, TradeImages::rreqs, MPI_STATUSES_IGNORE);
        if(retval != MPI_SUCCESS) rmg_error_handler (__FILE__, __LINE__, "Error in MPI_Waitall.\n");
    }
    T->thread_barrier_wait(true);


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
        retval = MPI_Waitall(26, TradeImages::sreqs, MPI_STATUSES_IGNORE);
        if(retval != MPI_SUCCESS) rmg_error_handler (__FILE__, __LINE__, "Error in MPI_Waitall.\n");
    }
    T->thread_barrier_wait(true);

} // end trade_images1_async


// Used by coarse level multigrid routines which assume that the central data is already present
// and that the f array is sized [dimx+2][dimy+2][dimz+2]
//
template <typename RmgType>
void TradeImages::trade_images1_central_async (RmgType * f, int dimx, int dimy, int dimz)
{
    BaseThread *T = BaseThread::getBaseThread(0);
    int tid = T->get_thread_tid();
    if(tid < 0) tid = 0;
    int ACTIVE_THREADS = 1;
    if(T->is_loop_over_states()) ACTIVE_THREADS = T->barrier->barrier_count();

    int istate = T->get_thread_basetag();
    MPI_Comm grid_comm = T->get_unique_comm(istate);

    int ix, iy, iz, incx, incy, index;
    int ixs2, iys2, c1, idx;
    int xlen, ylen, zlen;
    int retval;

    RmgType *frdx1_f, *frdx2_f, *frdy1_f, *frdy2_f, *frdz1_f, *frdz2_f;
    RmgType *frdx1n_f, *frdx2n_f, *frdy1n_f, *frdy2n_f, *frdz1n_f, *frdz2n_f;

    frdx1_f = (RmgType *)TradeImages::frdx1[0];
    frdx2_f = (RmgType *)TradeImages::frdx2[0];
    frdy1_f = (RmgType *)TradeImages::frdy1[0];
    frdy2_f = (RmgType *)TradeImages::frdy2[0];
    frdz1_f = (RmgType *)TradeImages::frdz1[0];
    frdz2_f = (RmgType *)TradeImages::frdz2[0];
    frdx1n_f = (RmgType *)TradeImages::frdx1n[0];
    frdx2n_f = (RmgType *)TradeImages::frdx2n[0];
    frdy1n_f = (RmgType *)TradeImages::frdy1n[0];
    frdy2n_f = (RmgType *)TradeImages::frdy2n[0];
    frdz1n_f = (RmgType *)TradeImages::frdz1n[0];
    frdz2n_f = (RmgType *)TradeImages::frdz2n[0];


    incx = (dimy + 2) * (dimz + 2);
    incy = dimz + 2;

    zlen = dimx * dimy;
    ylen = dimx * dimz;
    xlen = dimy * dimz;


    // Thread 0 posts all of the receives
    if(tid == 0) {

        // The z planes
      TradeImages::RMG_MPI_trade(frdz2n_f, ACTIVE_THREADS * zlen, RMG_MPI_IRECV, 0, 0, 1, grid_comm, 0, istate, &TradeImages::rreqs[0]);
      TradeImages::RMG_MPI_trade(frdz1n_f, ACTIVE_THREADS * zlen, RMG_MPI_IRECV, 0, 0, -1, grid_comm, 1, istate, &TradeImages::rreqs[1]);

        // The y planes
      TradeImages::RMG_MPI_trade(frdy2n_f, ACTIVE_THREADS * ylen, RMG_MPI_IRECV, 0, 1, 0, grid_comm, 2, istate, &TradeImages::rreqs[2]);
      TradeImages::RMG_MPI_trade(frdy1n_f, ACTIVE_THREADS * ylen, RMG_MPI_IRECV, 0, -1, 0, grid_comm, 3, istate, &TradeImages::rreqs[3]);

        // The x planes
      TradeImages::RMG_MPI_trade(frdx2n_f, ACTIVE_THREADS * xlen, RMG_MPI_IRECV, 1, 0, 0, grid_comm, 4, istate, &TradeImages::rreqs[4]);
      TradeImages::RMG_MPI_trade(frdx1n_f, ACTIVE_THREADS * xlen, RMG_MPI_IRECV, -1, 0, 0, grid_comm, 5, istate, &TradeImages::rreqs[5]);

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


    T->thread_barrier_wait(true);
    if(tid == 0) {

        // Send z planes
      TradeImages::RMG_MPI_trade(frdz1_f, ACTIVE_THREADS * zlen, RMG_MPI_ISEND, 0, 0, -1, grid_comm, 0, istate, &TradeImages::sreqs[0]);
      TradeImages::RMG_MPI_trade(frdz2_f, ACTIVE_THREADS * zlen, RMG_MPI_ISEND, 0, 0, 1, grid_comm, 1, istate, &TradeImages::sreqs[1]);

        // Send the north and south planes
      TradeImages::RMG_MPI_trade(frdy1_f, ACTIVE_THREADS * ylen, RMG_MPI_ISEND, 0, -1, 0, grid_comm, 2, istate, &TradeImages::sreqs[2]);
      TradeImages::RMG_MPI_trade(frdy2_f, ACTIVE_THREADS * ylen, RMG_MPI_ISEND, 0, 1, 0, grid_comm, 3, istate, &TradeImages::sreqs[3]);

        // Send the east and west planes
      TradeImages::RMG_MPI_trade(frdx1_f, ACTIVE_THREADS * xlen, RMG_MPI_ISEND, -1, 0, 0, grid_comm, 4, istate, &TradeImages::sreqs[4]);
      TradeImages::RMG_MPI_trade(frdx2_f, ACTIVE_THREADS * xlen, RMG_MPI_ISEND, 1, 0, 0, grid_comm, 5, istate, &TradeImages::sreqs[5]);

    }

    // Wait for all the recvs to finish
    if(tid == 0) {
        retval = MPI_Waitall(6, TradeImages::rreqs, MPI_STATUSES_IGNORE);
        if(retval != MPI_SUCCESS) rmg_error_handler (__FILE__, __LINE__, "Error in MPI_Waitall.\n");
    }
    T->thread_barrier_wait(true);


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
        retval = MPI_Waitall(6, TradeImages::sreqs, MPI_STATUSES_IGNORE);
        if(retval != MPI_SUCCESS) rmg_error_handler (__FILE__, __LINE__, "Error in MPI_Waitall.\n");
    }
    T->thread_barrier_wait(true);

} // end trade_images1_central_async


// Used by coarse level multigrid routines which assume that the central data is already present
// and that the f array is sized [dimx+2][dimy+2][dimz+2]
//
template <typename RmgType>
void TradeImages::trade_images1_central_async_managed (RmgType * f, int dimx, int dimy, int dimz)
{
    BaseThread *T = BaseThread::getBaseThread(0);
    int tid = T->get_thread_tid();
    if(tid < 0) tid = 0;
    int istate = T->get_thread_basetag();
    MPI_Comm grid_comm = T->get_unique_comm(istate);


    int ix, iy, iz, incx, incy, index;
    int ixs2, iys2, c1, idx;
    int xlen, ylen, zlen;

    std::atomic_bool is_completed_r[6];
    std::atomic_bool is_completed_s[6];

    incx = (dimy + 2) * (dimz + 2);
    incy = dimz + 2;

    zlen = dimx * dimy;
    ylen = dimx * dimz;
    xlen = dimy * dimz;

    mpi_queue_item_t qitems_r[6];
    mpi_queue_item_t qitems_s[6];

    for(int it = 0;it < 6;it++)
    {
        qitems_r[it].is_completed = &is_completed_r[it];
        qitems_s[it].is_completed = &is_completed_s[it];
    }

    RmgType *frdx1_f = (RmgType *)TradeImages::frdx1[tid];
    RmgType *frdx2_f = (RmgType *)TradeImages::frdx2[tid];
    RmgType *frdy1_f = (RmgType *)TradeImages::frdy1[tid];
    RmgType *frdy2_f = (RmgType *)TradeImages::frdy2[tid];
    RmgType *frdz1_f = (RmgType *)TradeImages::frdz1[tid];
    RmgType *frdz2_f = (RmgType *)TradeImages::frdz2[tid];
    RmgType *frdx1n_f = (RmgType *)TradeImages::frdx1n[tid];
    RmgType *frdx2n_f = (RmgType *)TradeImages::frdx2n[tid];
    RmgType *frdy1n_f = (RmgType *)TradeImages::frdy1n[tid];
    RmgType *frdy2n_f = (RmgType *)TradeImages::frdy2n[tid];
    RmgType *frdz1n_f = (RmgType *)TradeImages::frdz1n[tid];
    RmgType *frdz2n_f = (RmgType *)TradeImages::frdz2n[tid];


    // The z planes
    TradeImages::RMG_MPI_queue_trade(frdz2n_f, zlen, RMG_MPI_IRECV, 0, 0, 1, grid_comm, 0, istate, qitems_r[0]);
    TradeImages::RMG_MPI_queue_trade(frdz1n_f, zlen, RMG_MPI_IRECV, 0, 0, -1, grid_comm, 1, istate, qitems_r[1]);

    // The y planes
    TradeImages::RMG_MPI_queue_trade(frdy2n_f, ylen, RMG_MPI_IRECV, 0, 1, 0, grid_comm, 2, istate, qitems_r[2]);
    TradeImages::RMG_MPI_queue_trade(frdy1n_f, ylen, RMG_MPI_IRECV, 0, -1, 0, grid_comm, 3, istate, qitems_r[3]);

    // The x planes
    TradeImages::RMG_MPI_queue_trade(frdx2n_f, xlen, RMG_MPI_IRECV, 1, 0, 0, grid_comm, 4, istate, qitems_r[4]);
    TradeImages::RMG_MPI_queue_trade(frdx1n_f, xlen, RMG_MPI_IRECV, -1, 0, 0, grid_comm, 5, istate, qitems_r[5]);



    /* Collect the positive z-plane and negative z-planes */
    c1 = dimz-1;
    idx = 0;
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

    // Send z planes
    TradeImages::RMG_MPI_queue_trade(frdz1_f, zlen, RMG_MPI_ISEND, 0, 0, -1, grid_comm, 0, istate, qitems_s[0]);
    TradeImages::RMG_MPI_queue_trade(frdz2_f, zlen, RMG_MPI_ISEND, 0, 0, 1, grid_comm, 1, istate, qitems_s[1]);

    /* Collect the north and south planes */
    c1 = (dimy-1)*incy;
    idx = 0;
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


    // Send the north and south planes
    TradeImages::RMG_MPI_queue_trade(frdy1_f, ylen, RMG_MPI_ISEND, 0, -1, 0, grid_comm, 2, istate, qitems_s[2]);
    TradeImages::RMG_MPI_queue_trade(frdy2_f, ylen, RMG_MPI_ISEND, 0, 1, 0, grid_comm, 3, istate, qitems_s[3]);


    /* Collect the east and west planes */
    c1 = (dimx-1)*incx;
    idx = 0;
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

    // Send the east and west planes
    TradeImages::RMG_MPI_queue_trade(frdx1_f, xlen, RMG_MPI_ISEND, -1, 0, 0, grid_comm, 4, istate, qitems_s[4]);
    TradeImages::RMG_MPI_queue_trade(frdx2_f, xlen, RMG_MPI_ISEND, 1, 0, 0, grid_comm, 5, istate, qitems_s[5]);

    // Loop and unpack items as they become available
    int items_completed = 0;
    while(items_completed < 1)
    {
        if(!qitems_r[0].is_unpacked)
        {
            if(is_completed_r[0].load(std::memory_order_acquire))
            {
                /* Unpack z-planes */
                c1 = dimz + 1;
                idx = 0;
                for (ix = 0; ix < dimx; ix++)
                {
                    ixs2 = (ix + 1) * incx;
                    for (iy = 0; iy < dimy; iy++)
                    {
                        iys2 = ixs2 + (iy + 1) * incy;
                        for (iz = 0; iz < 1; iz++)
                        {
                            index = iys2 + iz;
                            f[index + c1] = frdz2n_f[idx];
                            idx++;
                        }                   /* end for */
                    }                       /* end for */
                }                           /* end for */
                items_completed++;
                qitems_r[0].is_unpacked = true;
            }
            else
            {
                if(!qitems_r[0].is_unpacked)this->queue->spinwait(10);
            }
        }
    }

    while(items_completed < 2)
    {
        if(!qitems_r[1].is_unpacked)
        {
            if(is_completed_r[1].load(std::memory_order_acquire))
            {
                /* Unpack z-planes */
                c1 = dimz + 1;
                idx = 0;
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
                            idx++;
                        }                   /* end for */
                    }                       /* end for */
                }                           /* end for */
                items_completed++;
                qitems_r[1].is_unpacked = true;
            }
            else
            {
                if(!qitems_r[1].is_unpacked)this->queue->spinwait(10);
            }
        }
    }

    // ***********************************************************

    while(items_completed < 3)
    {
        if(!qitems_r[2].is_unpacked)
        {
            if(is_completed_r[2].load(std::memory_order_acquire))
            {
                /* Unpack the north and south planes */
                c1 = (dimy + 1) * incy;
                idx = 0;
                for (ix = 0; ix < dimx; ix++)
                {
                    ixs2 = (ix + 1) * incx;
                    for (iy = 0; iy < 1; iy++)
                    {
                        iys2 = ixs2 + iy * incy;
                        for (iz = 0; iz < dimz; iz++)
                        {
                            index = iys2 + iz + 1;
                            f[index + c1] = frdy2n_f[idx];
                            idx++;
                        }                   /* end for */
                    }                       /* end for */
                }                           /* end for */
                items_completed++;
                qitems_r[2].is_unpacked = true;
            }
            else
            {
                if(!qitems_r[2].is_unpacked)this->queue->spinwait(10);
            }
        }
    }

    while(items_completed < 4)
    {
        if(!qitems_r[3].is_unpacked)
        {
            if(is_completed_r[3].load(std::memory_order_acquire))
            {
                /* Unpack the north and south planes */
                c1 = (dimy + 1) * incy;
                idx = 0;
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
                            idx++;
                        }                   /* end for */
                    }                       /* end for */
                }                           /* end for */
                items_completed++;
                qitems_r[3].is_unpacked = true;
            }
            else
            {
                if(!qitems_r[3].is_unpacked)this->queue->spinwait(10);
            }
        }
    }

    // ***********************************************************

    while(items_completed < 5)
    {
        if(!qitems_r[4].is_unpacked)
        {
            if(is_completed_r[4].load(std::memory_order_acquire))
            {
                /* Unpack the east and west planes */
                c1 = (dimx + 1) * incx;
                idx = 0;
                for (ix = 0; ix < 1; ix++)
                {
                    ixs2 = ix * incx;
                    for (iy = 0; iy < dimy; iy++)
                    {
                        iys2 = ixs2 + (iy + 1) * incy;
                        for (iz = 0; iz < dimz; iz++)
                        {
                            index = iys2 + iz + 1;
                            f[index + c1] = frdx2n_f[idx];
                            idx++;
                        }                   /* end for */
                    }                       /* end for */
                }                           /* end for */
                items_completed++;
                qitems_r[4].is_unpacked = true;
            }
            else
            {
                if(!qitems_r[4].is_unpacked)this->queue->spinwait(10);
            }
        }
    }

    while(items_completed < 6)
    {
        if(!qitems_r[5].is_unpacked)
        {
            if(is_completed_r[5].load(std::memory_order_acquire))
            {
                /* Unpack the east and west planes */
                c1 = (dimx + 1) * incx;
                idx = 0;
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
                            idx++;
                        }                   /* end for */
                    }                       /* end for */
                }                           /* end for */
                items_completed++;
                qitems_r[5].is_unpacked = true;
            }
            else
            {
                if(!qitems_r[5].is_unpacked)this->queue->spinwait(10);
            }
        }

    } // end while

   //this->queue->spinwaitall(is_completed_r, 6);
   this->queue->spinwaitall(is_completed_s, 6);

} // end trade_images1_central_async_managed


// Used by coarse level multigrid routines which assume that the central data is already present
// and that the f array is sized [dimx+2][dimy+2][dimz+2]
//
template <typename RmgType>
void TradeImages::trade_images1_async_managed (RmgType * f, int dimx, int dimy, int dimz)
{
    BaseThread *T = BaseThread::getBaseThread(0);
    int tid = T->get_thread_tid();
    if(tid < 0) tid = 0;
    int istate = T->get_thread_basetag();
    MPI_Comm grid_comm = T->get_unique_comm(istate);


    int ix, iy, iz, index;
    int ixs2, iys2, c1, c2, c3, idx;
    int node_idx;

    int incx = (dimy + 2) * (dimz + 2);
    int incy = dimz + 2;

    int zlen = dimx * dimy;
    int ylen = dimx * dimz;
    int xlen = dimy * dimz;

    int yzlen = dimx;
    int xylen = dimz;
    int xzlen = dimy;

    std::atomic_bool is_completed_r[26];
    std::atomic_bool is_completed_s[26];

    mpi_queue_item_t qitems_r[26];
    mpi_queue_item_t qitems_s[26];

    for(int it = 0;it < 26;it++)
    {
        qitems_r[it].is_completed = &is_completed_r[it];
        qitems_s[it].is_completed = &is_completed_s[it];
    }

    RmgType *frdx1_f = (RmgType *)TradeImages::frdx1[tid];
    RmgType *frdx2_f = (RmgType *)TradeImages::frdx2[tid];
    RmgType *frdy1_f = (RmgType *)TradeImages::frdy1[tid];
    RmgType *frdy2_f = (RmgType *)TradeImages::frdy2[tid];
    RmgType *frdz1_f = (RmgType *)TradeImages::frdz1[tid];
    RmgType *frdz2_f = (RmgType *)TradeImages::frdz2[tid];
    RmgType *frdx1n_f = (RmgType *)TradeImages::frdx1n[tid];
    RmgType *frdx2n_f = (RmgType *)TradeImages::frdx2n[tid];
    RmgType *frdy1n_f = (RmgType *)TradeImages::frdy1n[tid];
    RmgType *frdy2n_f = (RmgType *)TradeImages::frdy2n[tid];
    RmgType *frdz1n_f = (RmgType *)TradeImages::frdz1n[tid];
    RmgType *frdz2n_f = (RmgType *)TradeImages::frdz2n[tid];

    RmgType *yzpsms_r_f = (RmgType *)TradeImages::yzpsms_r[tid];
    RmgType *yzpsps_r_f = (RmgType *)TradeImages::yzpsps_r[tid];
    RmgType *yzmsms_r_f = (RmgType *)TradeImages::yzmsms_r[tid];
    RmgType *yzmsps_r_f = (RmgType *)TradeImages::yzmsps_r[tid];
    RmgType *yzpsms_s_f = (RmgType *)TradeImages::yzpsms_s[tid];
    RmgType *yzpsps_s_f = (RmgType *)TradeImages::yzpsps_s[tid];
    RmgType *yzmsms_s_f = (RmgType *)TradeImages::yzmsms_s[tid];
    RmgType *yzmsps_s_f = (RmgType *)TradeImages::yzmsps_s[tid];
    RmgType *xzpsms_r_f = (RmgType *)TradeImages::xzpsms_r[tid];
    RmgType *xzpsps_r_f = (RmgType *)TradeImages::xzpsps_r[tid];
    RmgType *xzmsms_r_f = (RmgType *)TradeImages::xzmsms_r[tid];
    RmgType *xzmsps_r_f = (RmgType *)TradeImages::xzmsps_r[tid];
    RmgType *xzpsms_s_f = (RmgType *)TradeImages::xzpsms_s[tid];
    RmgType *xzpsps_s_f = (RmgType *)TradeImages::xzpsps_s[tid];
    RmgType *xzmsms_s_f = (RmgType *)TradeImages::xzmsms_s[tid];
    RmgType *xzmsps_s_f = (RmgType *)TradeImages::xzmsps_s[tid];
    RmgType *xypsms_r_f = (RmgType *)TradeImages::xypsms_r[tid];
    RmgType *xypsps_r_f = (RmgType *)TradeImages::xypsps_r[tid];
    RmgType *xymsms_r_f = (RmgType *)TradeImages::xymsms_r[tid];
    RmgType *xymsps_r_f = (RmgType *)TradeImages::xymsps_r[tid];
    RmgType *xypsms_s_f = (RmgType *)TradeImages::xypsms_s[tid];
    RmgType *xypsps_s_f = (RmgType *)TradeImages::xypsps_s[tid];
    RmgType *xymsms_s_f = (RmgType *)TradeImages::xymsms_s[tid];
    RmgType *xymsps_s_f = (RmgType *)TradeImages::xymsps_s[tid];

    RmgType *m0_s_f = (RmgType *)TradeImages::m0_s[tid];
    RmgType *m0_r_f = (RmgType *)TradeImages::m0_r[tid];


    // Corner nodes
    TradeImages::RMG_MPI_queue_trade( &m0_r_f[0], 1, RMG_MPI_IRECV, -1, -1, -1, grid_comm, 18, istate, qitems_r[18]);
    TradeImages::RMG_MPI_queue_trade( &m0_r_f[1], 1, RMG_MPI_IRECV, -1, -1, 1, grid_comm, 19, istate, qitems_r[19]);
    TradeImages::RMG_MPI_queue_trade( &m0_r_f[2], 1, RMG_MPI_IRECV, -1, 1, -1, grid_comm, 20, istate, qitems_r[20]);
    TradeImages::RMG_MPI_queue_trade( &m0_r_f[3], 1, RMG_MPI_IRECV, -1, 1, 1, grid_comm, 21, istate, qitems_r[21]);
    TradeImages::RMG_MPI_queue_trade( &m0_r_f[4], 1, RMG_MPI_IRECV, 1, -1, -1, grid_comm, 22, istate, qitems_r[22]);
    TradeImages::RMG_MPI_queue_trade( &m0_r_f[5], 1, RMG_MPI_IRECV, 1, -1, 1, grid_comm, 23, istate, qitems_r[23]);
    TradeImages::RMG_MPI_queue_trade( &m0_r_f[6], 1, RMG_MPI_IRECV, 1, 1, -1, grid_comm, 24, istate, qitems_r[24]);
    TradeImages::RMG_MPI_queue_trade( &m0_r_f[7], 1, RMG_MPI_IRECV, 1, 1, 1, grid_comm, 25, istate, qitems_r[25]);

    // The z planes
    TradeImages::RMG_MPI_queue_trade(frdz2n_f,  zlen, RMG_MPI_IRECV, 0, 0, 1, grid_comm, 0, istate, qitems_r[0]);
    TradeImages::RMG_MPI_queue_trade(frdz1n_f,  zlen, RMG_MPI_IRECV, 0, 0, -1, grid_comm, 1, istate, qitems_r[1]);

    // The y planes
    TradeImages::RMG_MPI_queue_trade(frdy2n_f,  ylen, RMG_MPI_IRECV, 0, 1, 0, grid_comm, 2, istate, qitems_r[2]);
    TradeImages::RMG_MPI_queue_trade(frdy1n_f,  ylen, RMG_MPI_IRECV, 0, -1, 0, grid_comm, 3, istate, qitems_r[3]);

    // The x planes
    TradeImages::RMG_MPI_queue_trade(frdx2n_f,  xlen, RMG_MPI_IRECV, 1, 0, 0, grid_comm, 4, istate, qitems_r[4]);
    TradeImages::RMG_MPI_queue_trade(frdx1n_f,  xlen, RMG_MPI_IRECV, -1, 0, 0, grid_comm, 5, istate, qitems_r[5]);

    // And the yz-plane edges
    TradeImages::RMG_MPI_queue_trade(yzpsms_r_f,  yzlen, RMG_MPI_IRECV, 0, 1, -1, grid_comm, 6, istate, qitems_r[6]);
    TradeImages::RMG_MPI_queue_trade(yzpsps_r_f,  yzlen, RMG_MPI_IRECV, 0, 1, 1, grid_comm, 7, istate, qitems_r[7]);
    TradeImages::RMG_MPI_queue_trade(yzmsms_r_f,  yzlen, RMG_MPI_IRECV, 0, -1, -1, grid_comm, 8, istate, qitems_r[8]);
    TradeImages::RMG_MPI_queue_trade(yzmsps_r_f,  yzlen, RMG_MPI_IRECV, 0, -1, 1, grid_comm, 9, istate, qitems_r[9]);

    // And the xy-plane edges
    TradeImages::RMG_MPI_queue_trade(xypsms_r_f,  xylen, RMG_MPI_IRECV, 1, -1, 0, grid_comm, 10, istate, qitems_r[10]);
    TradeImages::RMG_MPI_queue_trade(xypsps_r_f,  xylen, RMG_MPI_IRECV, 1, 1, 0, grid_comm, 11, istate, qitems_r[11]);
    TradeImages::RMG_MPI_queue_trade(xymsms_r_f,  xylen, RMG_MPI_IRECV, -1, -1, 0, grid_comm, 12, istate, qitems_r[12]);
    TradeImages::RMG_MPI_queue_trade(xymsps_r_f,  xylen, RMG_MPI_IRECV, -1, 1, 0, grid_comm, 13, istate, qitems_r[13]);

    // And the xz-plane edges
    TradeImages::RMG_MPI_queue_trade(xzpsms_r_f,  xzlen, RMG_MPI_IRECV, 1, 0, -1, grid_comm, 14, istate, qitems_r[14]);
    TradeImages::RMG_MPI_queue_trade(xzpsps_r_f,  xzlen, RMG_MPI_IRECV, 1, 0, 1, grid_comm, 15, istate, qitems_r[15]);
    TradeImages::RMG_MPI_queue_trade(xzmsms_r_f,  xzlen, RMG_MPI_IRECV, -1, 0, -1, grid_comm, 16, istate, qitems_r[16]);
    TradeImages::RMG_MPI_queue_trade(xzmsps_r_f,  xzlen, RMG_MPI_IRECV, -1, 0, 1, grid_comm, 17, istate, qitems_r[17]);



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
                m0_s_f[node_idx] = f[c1 + c2 + c3];
                node_idx++;

            }
        }
    }


    // Corners
    TradeImages::RMG_MPI_queue_trade( &m0_s_f[0], 1, RMG_MPI_ISEND, 1, 1, 1, grid_comm, 18, istate, qitems_s[18]);
    TradeImages::RMG_MPI_queue_trade( &m0_s_f[1], 1, RMG_MPI_ISEND, 1, 1, -1, grid_comm, 19, istate, qitems_s[19]);
    TradeImages::RMG_MPI_queue_trade( &m0_s_f[2], 1, RMG_MPI_ISEND, 1, -1, 1, grid_comm, 20, istate, qitems_s[20]);
    TradeImages::RMG_MPI_queue_trade( &m0_s_f[3], 1, RMG_MPI_ISEND, 1, -1, -1, grid_comm, 21, istate, qitems_s[21]);
    TradeImages::RMG_MPI_queue_trade( &m0_s_f[4], 1, RMG_MPI_ISEND, -1, 1, 1, grid_comm, 22, istate, qitems_s[22]);
    TradeImages::RMG_MPI_queue_trade( &m0_s_f[5], 1, RMG_MPI_ISEND, -1, 1, -1, grid_comm, 23, istate, qitems_s[23]);
    TradeImages::RMG_MPI_queue_trade( &m0_s_f[6], 1, RMG_MPI_ISEND, -1, -1, 1, grid_comm, 24, istate, qitems_s[24]);
    TradeImages::RMG_MPI_queue_trade( &m0_s_f[7], 1, RMG_MPI_ISEND, -1, -1, -1, grid_comm, 25, istate, qitems_s[25]);


    /* Collect the positive z-plane and negative z-planes */
    c1 = dimz-1;
    idx = 0;
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

    // Send z planes
    TradeImages::RMG_MPI_queue_trade(frdz1_f,  zlen, RMG_MPI_ISEND, 0, 0, -1, grid_comm, 0, istate, qitems_s[0]);
    TradeImages::RMG_MPI_queue_trade(frdz2_f,  zlen, RMG_MPI_ISEND, 0, 0, 1, grid_comm, 1, istate, qitems_s[1]);

    /* Collect the north and south planes */
    c1 = (dimy-1)*incy;
    idx = 0;
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

    // Send the north and south planes
    TradeImages::RMG_MPI_queue_trade(frdy1_f,  ylen, RMG_MPI_ISEND, 0, -1, 0, grid_comm, 2, istate, qitems_s[2]);
    TradeImages::RMG_MPI_queue_trade(frdy2_f,  ylen, RMG_MPI_ISEND, 0, 1, 0, grid_comm, 3, istate, qitems_s[3]);

    /* Collect the east and west planes */
    c1 = (dimx-1)*incx;
    idx = 0;
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

    // Send the east and west planes
    TradeImages::RMG_MPI_queue_trade(frdx1_f,  xlen, RMG_MPI_ISEND, -1, 0, 0, grid_comm, 4, istate, qitems_s[4]);
    TradeImages::RMG_MPI_queue_trade(frdx2_f,  xlen, RMG_MPI_ISEND, 1, 0, 0, grid_comm, 5, istate, qitems_s[5]);


    // Pack yz-plane edges
    c1 = (dimy-1)*incy;
    c2 = dimz-1;
    idx = 0;
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

    // Send yz-plane edges
    TradeImages::RMG_MPI_queue_trade(yzpsms_s_f,  yzlen, RMG_MPI_ISEND, 0, -1, 1, grid_comm, 6, istate, qitems_s[6]);
    TradeImages::RMG_MPI_queue_trade(yzpsps_s_f,  yzlen, RMG_MPI_ISEND, 0, -1, -1, grid_comm, 7, istate, qitems_s[7]);
    TradeImages::RMG_MPI_queue_trade(yzmsms_s_f,  yzlen, RMG_MPI_ISEND, 0, 1, 1, grid_comm, 8, istate, qitems_s[8]);
    TradeImages::RMG_MPI_queue_trade(yzmsps_s_f,  yzlen, RMG_MPI_ISEND, 0, 1, -1, grid_comm, 9, istate, qitems_s[9]);

    // Pack xy plane edges
    c1 = (dimy-1)*incy;
    c2 = (dimx-1)*incx;
    idx = 0;
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

    // Send xy-plane edges
    TradeImages::RMG_MPI_queue_trade(xypsms_s_f,  xylen, RMG_MPI_ISEND, -1, 1, 0, grid_comm, 10, istate, qitems_s[10]);
    TradeImages::RMG_MPI_queue_trade(xypsps_s_f,  xylen, RMG_MPI_ISEND, -1, -1, 0, grid_comm, 11, istate, qitems_s[11]);
    TradeImages::RMG_MPI_queue_trade(xymsms_s_f,  xylen, RMG_MPI_ISEND, 1, 1, 0, grid_comm, 12, istate, qitems_s[12]);
    TradeImages::RMG_MPI_queue_trade(xymsps_s_f,  xylen, RMG_MPI_ISEND, 1, -1, 0, grid_comm, 13, istate, qitems_s[13]);

    // Pack xz plane edges
    c1 = (dimz-1);
    c2 = (dimx-1)*incx;
    idx = 0;
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

    // Send xz-plane edges
    TradeImages::RMG_MPI_queue_trade(xzpsms_s_f,  xzlen, RMG_MPI_ISEND, -1, 0, 1, grid_comm, 14, istate, qitems_s[14]);
    TradeImages::RMG_MPI_queue_trade(xzpsps_s_f,  xzlen, RMG_MPI_ISEND, -1, 0, -1, grid_comm, 15, istate, qitems_s[15]);
    TradeImages::RMG_MPI_queue_trade(xzmsms_s_f,  xzlen, RMG_MPI_ISEND, 1, 0, 1, grid_comm, 16, istate, qitems_s[16]);
    TradeImages::RMG_MPI_queue_trade(xzmsps_s_f,  xzlen, RMG_MPI_ISEND, 1, 0, -1, grid_comm, 17, istate, qitems_s[17]);


    // Loop and unpack items as they become available
    int items_completed = 0;
    while(items_completed < 8)
    {
        for(int i = 18;i < 26;i++)
        {
            if(!qitems_r[i].is_unpacked)
            {
                if(is_completed_r[i].load(std::memory_order_acquire))
                {
                    items_completed++;
                    // actually unpack them following the loop but since corners are only 1 item no overlap possible
                    qitems_r[i].is_unpacked = true;
                }
            }
        }
    }

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
                f[c1 + c2 + c3] = m0_r_f[node_idx];
                node_idx++;

            }
        }
    }

    while(items_completed < 20)
    {
        for(int i = 6;i < 18;i++)
        {
            if(!qitems_r[i].is_unpacked)
            {
                if(is_completed_r[i].load(std::memory_order_acquire))
                {
                    items_completed++;
                    // actually unpack them following the loop but edges are so small don't worry about overlap
                    qitems_r[i].is_unpacked = true;
                }
            }
        }
    }


    // Unpack yz-plane edges
    c1 = (dimy + 1) * incy;
    c2 = (dimz + 1);
    idx = 0;
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
    idx = 0;
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
    idx = 0;
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


    items_completed = 0;
    while(items_completed < 6)
    {

        if(!qitems_r[0].is_unpacked && !qitems_r[1].is_unpacked)
        {
            if(is_completed_r[0].load(std::memory_order_acquire) && is_completed_r[1].load(std::memory_order_acquire))
            {
                /* Unpack z-planes */
                c1 = dimz + 1;
                idx = 0;
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
                items_completed+=2;
                qitems_r[0].is_unpacked = true;
                qitems_r[1].is_unpacked = true;
            }
        }

        if(!qitems_r[2].is_unpacked && !qitems_r[3].is_unpacked)
        {
            if(is_completed_r[2].load(std::memory_order_acquire) && is_completed_r[3].load(std::memory_order_acquire))
            {
                /* Unpack the north and south planes */
                c1 = (dimy + 1) * incy;
                idx = 0;
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
                items_completed+=2;
                qitems_r[2].is_unpacked = true;
                qitems_r[3].is_unpacked = true;
            }
        }


        if(!qitems_r[4].is_unpacked && !qitems_r[5].is_unpacked)
        {
            if(is_completed_r[4].load(std::memory_order_acquire) && is_completed_r[5].load(std::memory_order_acquire))
            {
                /* Unpack the east and west planes */
                c1 = (dimx + 1) * incx;
                idx = 0;
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
                items_completed+=2;
                qitems_r[4].is_unpacked = true;
                qitems_r[5].is_unpacked = true;
            }
        }

    } // end while

   //this->queue->spinwaitall(is_completed_r, 26);
   this->queue->spinwaitall(is_completed_s, 26);

} // end trade_images1_async_managed


template <typename RmgType>
void TradeImages::trade_images_async_managed (RmgType * mat, int dimx, int dimy, int dimz)
{

    BaseThread *T = BaseThread::getBaseThread(0);
    int tid = T->get_thread_tid();
    if(tid < 0) tid = 0;
    int istate = T->get_thread_basetag();
    MPI_Comm grid_comm = T->get_unique_comm(istate);

    int idx;


/* precalc some boundaries */
    int incx = (dimy + 2) * (dimz + 2);
    int incy = (dimz + 2);
    int incz = 1;
    int xmax = dimx * incx;
    int ymax = dimy * incy;
    int zmax = dimz * incz;

    std::atomic_bool is_completed_r[6];
    std::atomic_bool is_completed_s[6];

    mpi_queue_item_t qitems_r[6];
    mpi_queue_item_t qitems_s[6];

    for(int it = 0;it < 6;it++)
    {
        qitems_r[it].is_completed = &is_completed_r[it];
        qitems_s[it].is_completed = &is_completed_s[it];
    }

    RmgType *frdy1_f = (RmgType *)TradeImages::frdy1[tid];
    RmgType *frdy2_f = (RmgType *)TradeImages::frdy2[tid];
    RmgType *frdz1_f = (RmgType *)TradeImages::frdz1[tid];
    RmgType *frdz2_f = (RmgType *)TradeImages::frdz2[tid];
    RmgType *frdy1n_f = (RmgType *)TradeImages::frdy1n[tid];
    RmgType *frdy2n_f = (RmgType *)TradeImages::frdy2n[tid];
    RmgType *frdz1n_f = (RmgType *)TradeImages::frdz1n[tid];
    RmgType *frdz2n_f = (RmgType *)TradeImages::frdz2n[tid];

    // The z planes
    TradeImages::RMG_MPI_queue_trade(frdz2n_f, dimx*dimy, RMG_MPI_IRECV, 0, 0, 1, grid_comm, 0, istate, qitems_r[0]);
    TradeImages::RMG_MPI_queue_trade(frdz1n_f, dimx*dimy, RMG_MPI_IRECV, 0, 0, -1, grid_comm, 1, istate, qitems_r[1]);

    // The y planes
    TradeImages::RMG_MPI_queue_trade(frdy2n_f, dimx*(dimz+2), RMG_MPI_IRECV, 0, 1, 0, grid_comm, 2, istate, qitems_r[2]);
    TradeImages::RMG_MPI_queue_trade(frdy1n_f, dimx*(dimz+2), RMG_MPI_IRECV, 0, -1, 0, grid_comm, 3, istate, qitems_r[3]);

    // The x planes
    TradeImages::RMG_MPI_queue_trade(&mat[xmax+incx], (dimy + 2) * (dimz + 2), RMG_MPI_IRECV, 1, 0, 0, grid_comm, 4, istate, qitems_r[4]);
    TradeImages::RMG_MPI_queue_trade(&mat[0], (dimy + 2) * (dimz + 2), RMG_MPI_IRECV, -1, 0, 0, grid_comm, 5, istate, qitems_r[5]);


/*
 * First, do Up-Down Trade (+/- z)
 *  (trading dimx * dimy array of data, twice)
 */

    idx = 0;
    for (int i = incx; i <= incx * dimx; i += incx)
    {
        for (int j = 1; j <= dimy; j++)
        {
            frdz2_f[idx] = mat[i + j * incy + zmax];
            frdz1_f[idx] = mat[i + j * incy + incz];
            idx++;
        }
    }


    // Send z planes
    TradeImages::RMG_MPI_queue_trade(frdz1_f,  dimx*dimy, RMG_MPI_ISEND, 0, 0, -1, grid_comm, 0, istate, qitems_s[0]);
    TradeImages::RMG_MPI_queue_trade(frdz2_f,  dimx*dimy, RMG_MPI_ISEND, 0, 0, 1, grid_comm, 1, istate, qitems_s[1]);


    // Wait for them to finish and then unpack
    int items_completed = 0;
    while(items_completed < 2)
    {
        if(!qitems_r[0].is_unpacked && is_completed_r[0].load(std::memory_order_acquire))
        {
            idx = 0;
            for (int i = incx; i <= incx * dimx; i += incx)
            {
                for (int j = 1; j <= dimy; j++)
                {
                    mat[i + j * incy + zmax + incz] = frdz2n_f[idx];
                    idx++;
                }
            }
            items_completed++;
            qitems_r[0].is_unpacked = true;
        }
        if(!qitems_r[1].is_unpacked && is_completed_r[1].load(std::memory_order_acquire))
        {
            idx = 0;
            for (int i = incx; i <= incx * dimx; i += incx)
            {
                for (int j = 1; j <= dimy; j++)
                {
                    mat[i + j * incy] = frdz1n_f[idx];
                    idx++;
                }
            }
            items_completed++;
            qitems_r[1].is_unpacked = true;
        }
    }


/*
 * next, do North-South Trade (+/- y)
 *  (trading dimx * (dimz+2) array of data, twice)
 */

    idx = 0;
    for (int ix = 0; ix < dimx; ix++) {
        int iys2 = (ix + 1) * incx;
        for (int iz = 0; iz < dimz + 2; iz++) {
            frdy2_f[idx] = mat[iys2 + ymax + iz];
            frdy1_f[idx] = mat[iys2 + incy + iz];
            idx++;
        }
    }

    // Send y planes
    TradeImages::RMG_MPI_queue_trade(frdy1_f, dimx*(dimz+2), RMG_MPI_ISEND, 0, -1, 0, grid_comm, 2, istate, qitems_s[2]);
    TradeImages::RMG_MPI_queue_trade(frdy2_f, dimx*(dimz+2), RMG_MPI_ISEND, 0, 1, 0, grid_comm, 3, istate, qitems_s[3]);

    // Wait for them to finish and then unpack
    items_completed = 0;
    while(items_completed < 2)
    {
        if(!qitems_r[2].is_unpacked)
        {
            if(is_completed_r[2].load(std::memory_order_acquire))
            {
                int idx = 0;
                for (int ix = 0; ix < dimx; ix++) {
                    int iys2 = (ix + 1) * incx;
                    for (int iz = 0; iz < dimz + 2; iz++) {
                        mat[iys2 + ymax + incy + iz] = frdy2n_f[idx];
                        idx++;
                    }
                }
                items_completed++;
                qitems_r[2].is_unpacked = true;
            }
        }
        if(!qitems_r[3].is_unpacked)
        {
            if(is_completed_r[3].load(std::memory_order_acquire))
            {
                int idx = 0;
                for (int ix = 0; ix < dimx; ix++) {
                    int iys2 = (ix + 1) * incx;
                    for (int iz = 0; iz < dimz + 2; iz++) {
                        mat[iys2 + iz] = frdy1n_f[idx];
                        idx++;
                    }
                }
                items_completed++;
                qitems_r[3].is_unpacked = true;
            }
        }
    }


/*
 * Finally, do East - West Trade (+/- x)
 *  (trading (dimy+2) * (dimz+2) array of data, twice)
 */

    // Send x planes in place, no unpacking required
    TradeImages::RMG_MPI_queue_trade(&mat[incx], (dimy+2)*(dimz+2), RMG_MPI_ISEND, -1, 0, 0, grid_comm, 4, istate, qitems_s[4]);
    TradeImages::RMG_MPI_queue_trade(&mat[xmax], (dimy+2)*(dimz+2), RMG_MPI_ISEND, 1, 0, 0, grid_comm, 5, istate, qitems_s[5]);

   this->queue->spinwaitall(is_completed_r, 6);
   this->queue->spinwaitall(is_completed_s, 6);

}
