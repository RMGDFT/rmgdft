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

#ifndef RMG_TradeImages_H
#define RMG_TradeImages_H 1

#include <mpi.h>
#include "BaseGrid.h"
#include "rmg_error.h"
#include "MpiQueue.h"



/* Type of image trading */
#define FULL_TRADE 1
#define CENTRAL_TRADE 2

/* Maximum number of images for finite difference routines */
#define MAX_TRADE_IMAGES 6

#if __cplusplus

#include "BaseThread.h"
#include <boost/lockfree/queue.hpp>



class TradeImages {


private:

    BaseGrid *G;
    MpiQueue *queue;
    void allocate_buffers(double ** &P, int nthreads, int length_per_thread, size_t elem_len);
    bool queue_mode;

    /// Synchronous/asynchronous mode. 0=asnychronous (default) 1=synchronous
    int mode;

    // Timer mode 0=off (default) 1=on
    bool timer_mode;

    // Is array local? (Only 1 MPI process)
    bool local_mode;

    // rank of this node in comm
    int gridpe;

    /// Rank of target node based on offsets from current node. Used by asynchronous comm routines.
    int target_node[3][3][3];

    int max_alloc;

    /// Buffers allocated via MPI_Allocmem
    double *swbuf1x;
    double *swbuf2x;

    double **frdx1, **frdx2, **frdy1, **frdy2, **frdz1, **frdz2;
    double **frdx1n, **frdx2n, **frdy1n, **frdy2n, **frdz1n, **frdz2n;
    double **yzpsms_r, **yzpsps_r, **yzmsms_r, **yzmsps_r;
    double **yzpsms_s, **yzpsps_s, **yzmsms_s, **yzmsps_s;
    double **xzpsms_r, **xzpsps_r, **xzmsms_r, **xzmsps_r;
    double **xzpsms_s, **xzpsps_s, **xzmsms_s, **xzmsps_s;
    double **xypsms_r, **xypsps_r, **xymsms_r, **xymsps_r;
    double **xypsms_s, **xypsps_s, **xymsms_s, **xymsps_s;
    double **m0_s, **m0_r;

    MPI_Request sreqs[26];
    MPI_Request rreqs[26];

    void init_trade_imagesx_async(size_t elem_len);
    template <typename RmgType> void RMG_MPI_trade(RmgType *buf, int count, int type, int pe_x_offset, int pe_y_offset, int pe_z_offset, MPI_Comm comm, int tag, int extra_tag, MPI_Request *req);

    template <typename RmgType> void RMG_MPI_queue_trade(RmgType *buf, int count, int type, int pe_x_offset, int pe_y_offset, int pe_z_offset, MPI_Comm comm, int tag, int extra_tag, mpi_queue_item_t &qitem);

    template <typename RmgType> void RMG_MPI_queue_allreduce(RmgType *buf, int count, MPI_Datatype type, MPI_Comm comm, mpi_queue_item_t &qitem);

    template <typename RmgType> void trade_imagesx_async (RmgType * f, RmgType * w, int dimx, int dimy, int dimz, int images);
    template <typename RmgType> void trade_imagesx_central_async (RmgType * f, RmgType * w, int dimx, int dimy, int dimz, int images);
    template <typename RmgType> void trade_imagesx_central_local (RmgType * f, RmgType * w, int dimx, int dimy, int dimz, int images);
    template <typename RmgType> void trade_imagesx_central_async_managed (RmgType * f, RmgType * w, int dimx, int dimy, int dimz, int images);
    template <typename RmgType> void trade_images1_central_async (RmgType * f, int dimx, int dimy, int dimz);
    template <typename RmgType> void trade_images1_async (RmgType * f, int dimx, int dimy, int dimz);
    template <typename RmgType> void trade_images1_async_managed (RmgType * f, int dimx, int dimy, int dimz);
    template <typename RmgType> void trade_images1_central_async_managed (RmgType * f, int dimx, int dimy, int dimz);
    template <typename RmgType> void trade_images_async_managed (RmgType * f, int dimx, int dimy, int dimz);
    template <typename RmgType> void trade_images_local (RmgType * mat, int dimx, int dimy, int dimz, int type);



public:
    /// MPI communicator to use
    MPI_Comm comm;

    TradeImages(BaseGrid *BG, size_t elem_len, bool new_queue_mode, MpiQueue *newQM);
    ~TradeImages(void);
    void set_synchronous_mode(void);
    void set_asynchronous_mode(void);
    void set_queue_mode(bool mode);
    void set_timer_mode(bool verbose);
    void set_MPI_comm(MPI_Comm comm);
    MPI_Comm get_MPI_comm(void);
    void set_gridpe(int gridpe);
    template <typename RmgType> void trade_imagesx (RmgType *f, RmgType *w, int dimx, int dimy, int dimz, int images, int type);
    template <typename RmgType> void trade_images (RmgType * mat, int dimx, int dimy, int dimz, int type);




};

#endif
#endif
