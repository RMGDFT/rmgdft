/************************** SVN Revision Information **************************
 **    $Id: trade_imagesx.c 2091 2014-02-05 15:48:31Z ebriggs $    **
******************************************************************************/

#include "TradeImages.h"
#include "transition.h"
#include "common_prototypes.h"
#include "BlasWrappers.h"
#include <cmath>
#include <complex>

using namespace std;

#define ASYNC_MODE      0
#define SYNC_MODE       1

/* Type of async request passed to the mpi trade_images manager */
#define RMG_MPI_ISEND 1
#define RMG_MPI_IRECV 2


// Force instantiation of float, double and complex versions.
template void TradeImages::trade_images<rmg_float_t>(rmg_float_t*, int, int, int, int);
template void TradeImages::trade_images<rmg_double_t>(rmg_double_t*, int, int, int, int);
template void TradeImages::trade_imagesx<rmg_float_t>(rmg_float_t*, rmg_float_t*, int, int, int, int, int);
template void TradeImages::trade_imagesx<rmg_double_t>(rmg_double_t*, rmg_double_t*, int, int, int, int, int);
template void TradeImages::trade_imagesx<complex<float> >(complex <float>*, complex <float>*, int, int, int, int, int);
template void TradeImages::trade_imagesx<complex<double> >(complex <double>*, complex <double>*, int, int, int, int, int);

extern "C"
{
rmg_double_t my_crtc (void);
}

/*
 * INPUTS
 * f[dimx*dimy*dimz] - raw array without images. pack_ptos should not be called, this is 
 *                     handled inside the function now
 
 * w[(dimx+2*images) * (dimy+2*images) * (dimz+2*images)]
 *    This array is completely filled, i.e. the original data is filled in and then 
 *    the image data are added*/

rmg_double_t *TradeImages::swbuf1x = NULL;
rmg_double_t *TradeImages::swbuf2x = NULL;
int TradeImages::max_alloc;

// Rank of target node based on offsets from current node
int TradeImages::target_node[3][3][3];
int TradeImages::mode;


rmg_double_t *TradeImages::frdx1, *TradeImages::frdx2, *TradeImages::frdy1, *TradeImages::frdy2, *TradeImages::frdz1, *TradeImages::frdz2;
rmg_double_t *TradeImages::frdx1n, *TradeImages::frdx2n, *TradeImages::frdy1n, *TradeImages::frdy2n, *TradeImages::frdz1n, *TradeImages::frdz2n;
rmg_double_t *TradeImages::yzpsms_r, *TradeImages::yzpsps_r, *TradeImages::yzmsms_r, *TradeImages::yzmsps_r;
rmg_double_t *TradeImages::yzpsms_s, *TradeImages::yzpsps_s, *TradeImages::yzmsms_s, *TradeImages::yzmsps_s;
rmg_double_t *TradeImages::xzpsms_r, *TradeImages::xzpsps_r, *TradeImages::xzmsms_r, *TradeImages::xzmsps_r;
rmg_double_t *TradeImages::xzpsms_s, *TradeImages::xzpsps_s, *TradeImages::xzmsms_s, *TradeImages::xzmsps_s;
rmg_double_t *TradeImages::xypsms_r, *TradeImages::xypsps_r, *TradeImages::xymsms_r, *TradeImages::xymsps_r;
rmg_double_t *TradeImages::xypsms_s, *TradeImages::xypsps_s, *TradeImages::xymsms_s, *TradeImages::xymsps_s;
rmg_double_t *TradeImages::m0_s, *TradeImages::m0_r;

MPI_Request TradeImages::TradeImages::sreqs[26];
MPI_Request TradeImages::TradeImages::rreqs[26];

// Constructor
TradeImages::TradeImages(void)
{

    if((swbuf1x == NULL) && (swbuf2x == NULL)) {
        int grid_xp, grid_yp, grid_zp, grid_max1, grid_max2, retval;
        BaseGrid G;
        BaseThread T(0);

        grid_xp = G.get_FPX0_GRID() + 2*MAX_TRADE_IMAGES;
        grid_yp = G.get_FPY0_GRID() + 2*MAX_TRADE_IMAGES;
        grid_zp = G.get_FPZ0_GRID() + 2*MAX_TRADE_IMAGES;
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
         retval = MPI_Alloc_mem(6 * sizeof(complex<double>) * MAX_TRADE_IMAGES * T.get_threads_per_node() * grid_max1*grid_max2 , MPI_INFO_NULL, &swbuf1x);
         if(retval != MPI_SUCCESS) {
             rmg_error_handler (__FILE__, __LINE__, "Error in MPI_Alloc_mem.\n");
         }
         retval = MPI_Alloc_mem(6 * sizeof(complex<double>) * MAX_TRADE_IMAGES * T.get_threads_per_node() * grid_max1*grid_max2 , MPI_INFO_NULL, &swbuf2x);
         if(retval != MPI_SUCCESS) {
             rmg_error_handler (__FILE__, __LINE__, "Error in MPI_Alloc_mem.\n");
         }


         //my_malloc(swbuf1x, 6 * MAX_TRADE_IMAGES * grid_max1 * grid_max2 * ct.THREADS_PER_NODE, rmg_double_t);
         //my_malloc(swbuf2x, 6 * MAX_TRADE_IMAGES * grid_max1 * grid_max2 * ct.THREADS_PER_NODE, rmg_double_t);
         max_alloc = 6 * MAX_TRADE_IMAGES * grid_max1 * grid_max2 * T.get_threads_per_node();

         TradeImages::init_trade_imagesx_async();        

         return;
    }

}

void TradeImages::set_synchronous_mode(void)
{
    TradeImages::mode = SYNC_MODE;
}
void TradeImages::set_asynchronous_mode(void)
{
    TradeImages::mode = ASYNC_MODE;
}

template <typename RmgType>
void TradeImages::trade_imagesx (RmgType *f, RmgType *w, int dimx, int dimy, int dimz, int images, int type)
{
    int ix, iy, iz, incx, incy, incx0, incy0, index, tim, ione = 1;
    int ixs, iys, ixs2, iys2, c1, c2, alloc;
    int xlen, ylen, zlen, stop, tid;
    int *nb_ids, factor;
    MPI_Status mrstatus;
    RmgType *frdx1, *frdx2, *frdy1, *frdy2, *frdz1, *frdz2;
    RmgType *frdx1n, *frdx2n, *frdy1n, *frdy2n, *frdz1n, *frdz2n;
    RmgType *swbuf1x_f, *swbuf2x_f;
    int ACTIVE_THREADS = 1;
    BaseThread T(0);
    BaseGrid G;

    factor = sizeof(RmgType);


    if(TradeImages::mode == ASYNC_MODE) {
        if(type == CENTRAL_TRADE) {
            TradeImages::trade_imagesx_central_async (f, w, dimx, dimy, dimz, images);
            return;
        }
        else {
            TradeImages::trade_imagesx_async (f, w, dimx, dimy, dimz, images);
            return;
        }
    }

    tid = 0;
    tid = T.get_thread_tid();
    if(tid < 0) tid = 0;
    if(T.is_loop_over_states()) ACTIVE_THREADS = T.get_threads_per_node();

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

    alloc = 2 * (xlen + ylen + zlen) * T.get_threads_per_node();
    // Verify that the required memory will fit into the statically allocated storage
    if(alloc > TradeImages::max_alloc)
        rmg_error_handler (__FILE__, __LINE__, "Not enough memory. This should never happen.");


    nb_ids = G.get_neighbors();

    /* Load up w with the basic stuff */
    for (ix = 0; ix < dimx; ix++)
    {
        ixs = ix * incx0;
        ixs2 = (ix + images) * incx;

        for (iy = 0; iy < dimy; iy++)
        {
            iys = ixs + iy * incy0;
            iys2 = ixs2 + (iy + images) * incy;

            QMD_copy (dimz, &f[iys], ione, &w[iys2 + images], ione);

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
    T.thread_barrier_wait();
    if(tid == 0) {
        MPI_Sendrecv (&swbuf1x_f[0], factor*stop, MPI_BYTE, nb_ids[NB_D], (1>>16),
                  &swbuf2x_f[0], factor*stop, MPI_BYTE, nb_ids[NB_U], (1>>16), transition_get_grid_comm(), &mrstatus);

        MPI_Sendrecv (&swbuf1x_f[stop], factor*stop, MPI_BYTE, nb_ids[NB_U], (2>>16),
                  &swbuf2x_f[stop], factor*stop, MPI_BYTE, nb_ids[NB_D], (2>>16), transition_get_grid_comm(), &mrstatus);
    }
    T.thread_barrier_wait();

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
    T.thread_barrier_wait();
    if(tid == 0) {
        MPI_Sendrecv (&swbuf1x_f[0], factor*stop, MPI_BYTE, nb_ids[NB_S], (3>>16),
                  &swbuf2x_f[0], factor*stop, MPI_BYTE, nb_ids[NB_N], (3>>16), transition_get_grid_comm(), &mrstatus);

        MPI_Sendrecv (&swbuf1x_f[stop], factor*stop, MPI_BYTE, nb_ids[NB_N], (4>>16),
                  &swbuf2x_f[stop], factor*stop, MPI_BYTE, nb_ids[NB_S], (4>>16), transition_get_grid_comm(), &mrstatus);
    }
    T.thread_barrier_wait();

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
    T.thread_barrier_wait();
    if(tid == 0) {

        MPI_Sendrecv (&swbuf1x_f[0], factor*stop, MPI_BYTE, nb_ids[NB_W], (5>>16),
                   &swbuf2x_f[0], factor*stop, MPI_BYTE, nb_ids[NB_E], (5>>16), transition_get_grid_comm(), &mrstatus);

        MPI_Sendrecv (&swbuf1x_f[stop], factor*stop, MPI_BYTE, nb_ids[NB_E], (6>>16),
                  &swbuf2x_f[stop], factor*stop, MPI_BYTE, nb_ids[NB_W], (6>>16), transition_get_grid_comm(), &mrstatus);

    }
    T.thread_barrier_wait();


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


} // end trade_imagesx



template <typename RmgType>
void TradeImages::trade_images (RmgType * mat, int dimx, int dimy, int dimz, int type)
{
    int i, j, ione=1;
    int incx, incy, incz;
    int xmax, ymax, zmax;
    int ix, iz, iys2;
    int alloc, alloc1, toffset;
    int ACTIVE_THREADS = 1, factor;
    int *nb_ids;
    BaseGrid G;
    BaseThread T(0);

    factor = sizeof(RmgType);

    MPI_Status mstatus;
    int idx, stop, tid=0;
    RmgType *nmat1, *nmat2;
    RmgType *swbuf1x_f, *swbuf2x_f;


    swbuf1x_f = (RmgType *)TradeImages::swbuf1x;
    swbuf2x_f = (RmgType *)TradeImages::swbuf2x;

    nb_ids = G.get_neighbors();

    if(TradeImages::mode == ASYNC_MODE) {
        if(type == CENTRAL_TRADE) {
            TradeImages::trade_images1_central_async (mat, dimx, dimy, dimz);
            return;
        }
    }

    alloc = 4 * (dimx + 2) * (dimy + 2);
    alloc1 = 4 * (dimy + 2) * (dimz + 2);
    if(alloc1 > alloc)
        alloc = alloc1;
    alloc = alloc * T.get_threads_per_node();
    if(alloc > TradeImages::max_alloc)
        rmg_error_handler (__FILE__, __LINE__, "Not enough memory. This should never happen.");


    tid = T.get_thread_tid();
    if(tid < 0) tid = 0;
    if(T.is_loop_over_states()) ACTIVE_THREADS = T.get_threads_per_node();


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
    T.thread_barrier_wait();
    if(tid == 0) {
        MPI_Sendrecv (swbuf1x_f, factor*stop, MPI_BYTE, nb_ids[NB_U], (1>>16), swbuf2x_f, factor*stop,
		  MPI_BYTE, nb_ids[NB_D], (1>>16), transition_get_grid_comm(), &mstatus);
    } 
    T.thread_barrier_wait();

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
    T.thread_barrier_wait();
    if(tid == 0) {
        MPI_Sendrecv (swbuf1x_f, factor*stop, MPI_BYTE, nb_ids[NB_D], (2>>16), swbuf2x_f, factor*stop,
                      MPI_BYTE, nb_ids[NB_U], (2>>16), transition_get_grid_comm(), &mstatus);
    }
    T.thread_barrier_wait();


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

    T.thread_barrier_wait();
    if(tid == 0) {
        MPI_Sendrecv (swbuf1x_f, factor*stop * ACTIVE_THREADS, MPI_BYTE, nb_ids[NB_N], (3>>16), swbuf2x_f, factor*stop * ACTIVE_THREADS,
                  MPI_BYTE, nb_ids[NB_S], (3>>16), transition_get_grid_comm(), &mstatus);
    }
    T.thread_barrier_wait();

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


    T.thread_barrier_wait();
    if(tid == 0) {
        MPI_Sendrecv (swbuf1x_f, factor*stop * ACTIVE_THREADS, MPI_BYTE, nb_ids[NB_S], (4>>16), swbuf2x_f, factor*stop * ACTIVE_THREADS,
                  MPI_BYTE, nb_ids[NB_N], (4>>16), transition_get_grid_comm(), &mstatus);
    }
    T.thread_barrier_wait();


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

    QMD_copy (stop, &mat[xmax], ione, &swbuf1x_f[tid * stop], ione);


    T.thread_barrier_wait();
    if(tid == 0) {
        MPI_Sendrecv (swbuf1x_f, factor*stop * ACTIVE_THREADS, MPI_BYTE, nb_ids[NB_E], (5>>16), swbuf2x_f, factor*stop * ACTIVE_THREADS,
                  MPI_BYTE, nb_ids[NB_W], (5>>16), transition_get_grid_comm(), &mstatus);
    }
    T.thread_barrier_wait();

    QMD_copy (stop, &swbuf2x_f[tid * stop], ione, mat, ione);
    QMD_copy (stop, &mat[incx], ione, &swbuf1x_f[tid * stop], ione);

    T.thread_barrier_wait();
    if(tid == 0) {
        MPI_Sendrecv (swbuf1x_f, factor*stop * ACTIVE_THREADS, MPI_BYTE, nb_ids[NB_W], (6>>16), swbuf2x_f, factor*stop * ACTIVE_THREADS,
                  MPI_BYTE, nb_ids[NB_E], (6>>16), transition_get_grid_comm(), &mstatus);
    }
    T.thread_barrier_wait();
    QMD_copy (stop, &swbuf2x_f[tid * stop], ione, &mat[xmax + incx], ione);


    /* For clusters set the boundaries to zero -- this is wrong for the hartree
     * potential but we'll fix it up later. */
//    if ((ct.boundaryflag == CLUSTER) || (ct.boundaryflag == SURFACE))
//        set_bc (mat, dimx, dimy, dimz, 1, 0.0);


}




int TRADE_GRID_EDGES;
int GRID_MAX1;
int GRID_MAX2;



#define MAX_IMG2 (MAX_TRADE_IMAGES*MAX_TRADE_IMAGES)
#define MAX_IMG3 (MAX_TRADE_IMAGES*MAX_TRADE_IMAGES*MAX_TRADE_IMAGES)



// This function is used to setup the MPI request
template <typename RmgType>
void TradeImages::RMG_MPI_trade(RmgType *buf, int count, int type, int pe_x_offset, int pe_y_offset, int pe_z_offset, MPI_Comm comm, int tag, MPI_Request *req)
{
    int tid=0, ntag, target, factor;
    BaseThread T(0);

    factor = sizeof(RmgType);

    // Incement offsets so they can act as array indices into send and recv lists
    pe_x_offset++;
    pe_y_offset++;
    pe_z_offset++;

    // Tag is based on tid in the lower 8 bits which gives us up to 256 threads
    tid = T.get_thread_tid();
    if(tid == -1) tid = 0;


    ntag = (tag<<8) + tid;
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


// Allocates memory via MPI_Alloc_mem for use on systems with RDMA capability
void TradeImages::init_trade_imagesx_async(void) 
{
    int retval, THREADS_PER_NODE;
    int ix, iy, iz;
    int pe_x, pe_y, pe_z;
    int t_pe_x, t_pe_y, t_pe_z;
    int grid_xp, grid_yp, grid_zp;
    BaseThread T(0);
    BaseGrid G;

    THREADS_PER_NODE = T.get_threads_per_node();

    printf("Using Async trade_images with max images = %d.\n", MAX_TRADE_IMAGES);

    grid_xp = G.get_FPX0_GRID() + 2*MAX_TRADE_IMAGES;
    grid_yp = G.get_FPY0_GRID() + 2*MAX_TRADE_IMAGES;
    grid_zp = G.get_FPZ0_GRID() + 2*MAX_TRADE_IMAGES;
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
    pe2xyz(transition_get_gridpe(), &pe_x, &pe_y, &pe_z);
    for(ix = -1;ix <= 1;ix++) {

        t_pe_x = (pe_x + ix + G.PE_X) % G.PE_X;

        for(iy = -1;iy <= 1;iy++) {

            t_pe_y = (pe_y + iy + G.PE_Y) % G.PE_Y;

            for(iz = -1;iz <= 1;iz++) {

                t_pe_z = (pe_z + iz + G.PE_Z) % G.PE_Z;
                TradeImages::target_node[ix+1][iy+1][iz+1] = t_pe_x*G.PE_Y*G.PE_Z + t_pe_y*G.PE_Z + t_pe_z;

            }
        }
    } // end for


    // Allocate memory buffers using MPI_Alloc_mem
    retval = MPI_Alloc_mem(sizeof(complex<double>) * MAX_TRADE_IMAGES * THREADS_PER_NODE * GRID_MAX1 * GRID_MAX2 , MPI_INFO_NULL, &TradeImages::frdx1);
    if(retval != MPI_SUCCESS) {
        rmg_error_handler (__FILE__, __LINE__, "Error in MPI_Alloc_mem.\n");
    }
    retval = MPI_Alloc_mem(sizeof(complex<double>) * MAX_TRADE_IMAGES * THREADS_PER_NODE * GRID_MAX1 * GRID_MAX2 , MPI_INFO_NULL, &TradeImages::frdx2);
    if(retval != MPI_SUCCESS) {
        rmg_error_handler (__FILE__, __LINE__, "Error in MPI_Alloc_mem.\n");
    }
    retval = MPI_Alloc_mem(sizeof(complex<double>) * MAX_TRADE_IMAGES * THREADS_PER_NODE * GRID_MAX1 * GRID_MAX2 , MPI_INFO_NULL, &TradeImages::frdy1);
    if(retval != MPI_SUCCESS) {
        rmg_error_handler (__FILE__, __LINE__, "Error in MPI_Alloc_mem.\n");
    }
    retval = MPI_Alloc_mem(sizeof(complex<double>) * MAX_TRADE_IMAGES * THREADS_PER_NODE * GRID_MAX1 * GRID_MAX2 , MPI_INFO_NULL, &TradeImages::frdy2);
    if(retval != MPI_SUCCESS) {
        rmg_error_handler (__FILE__, __LINE__, "Error in MPI_Alloc_mem.\n");
    }
    retval = MPI_Alloc_mem(sizeof(complex<double>) * MAX_TRADE_IMAGES * THREADS_PER_NODE * GRID_MAX1 * GRID_MAX2 , MPI_INFO_NULL, &TradeImages::frdz1);
    if(retval != MPI_SUCCESS) {
        rmg_error_handler (__FILE__, __LINE__, "Error in MPI_Alloc_mem.\n");
    }
    retval = MPI_Alloc_mem(sizeof(complex<double>) * MAX_TRADE_IMAGES * THREADS_PER_NODE * GRID_MAX1 * GRID_MAX2 , MPI_INFO_NULL, &TradeImages::frdz2);
    if(retval != MPI_SUCCESS) {
        rmg_error_handler (__FILE__, __LINE__, "Error in MPI_Alloc_mem.\n");
    }
    retval = MPI_Alloc_mem(sizeof(complex<double>) * MAX_TRADE_IMAGES * THREADS_PER_NODE * GRID_MAX1 * GRID_MAX2 , MPI_INFO_NULL, &TradeImages::frdx1n);
    if(retval != MPI_SUCCESS) {
        rmg_error_handler (__FILE__, __LINE__, "Error in MPI_Alloc_mem.\n");
    }
    retval = MPI_Alloc_mem(sizeof(complex<double>) * MAX_TRADE_IMAGES * THREADS_PER_NODE * GRID_MAX1 * GRID_MAX2 , MPI_INFO_NULL, &TradeImages::frdx2n);
    if(retval != MPI_SUCCESS) {
        rmg_error_handler (__FILE__, __LINE__, "Error in MPI_Alloc_mem.\n");
    }
    retval = MPI_Alloc_mem(sizeof(complex<double>) * MAX_TRADE_IMAGES * THREADS_PER_NODE * GRID_MAX1 * GRID_MAX2 , MPI_INFO_NULL, &TradeImages::frdy1n);
    if(retval != MPI_SUCCESS) {
        rmg_error_handler (__FILE__, __LINE__, "Error in MPI_Alloc_mem.\n");
    }
    retval = MPI_Alloc_mem(sizeof(complex<double>) * MAX_TRADE_IMAGES * THREADS_PER_NODE * GRID_MAX1 * GRID_MAX2 , MPI_INFO_NULL, &TradeImages::frdy2n);
    if(retval != MPI_SUCCESS) {
        rmg_error_handler (__FILE__, __LINE__, "Error in MPI_Alloc_mem.\n");
    }
    retval = MPI_Alloc_mem(sizeof(complex<double>) * MAX_TRADE_IMAGES * THREADS_PER_NODE * GRID_MAX1 * GRID_MAX2 , MPI_INFO_NULL, &TradeImages::frdz1n);
    if(retval != MPI_SUCCESS) {
        rmg_error_handler (__FILE__, __LINE__, "Error in MPI_Alloc_mem.\n");
    }
    retval = MPI_Alloc_mem(sizeof(complex<double>) * MAX_TRADE_IMAGES * THREADS_PER_NODE * GRID_MAX1 * GRID_MAX2 , MPI_INFO_NULL, &TradeImages::frdz2n);
    if(retval != MPI_SUCCESS) {
        rmg_error_handler (__FILE__, __LINE__, "Error in MPI_Alloc_mem.\n");
    }
    retval = MPI_Alloc_mem(sizeof(complex<double>) * 8 * MAX_IMG2 * THREADS_PER_NODE * TRADE_GRID_EDGES , MPI_INFO_NULL, &TradeImages::yzpsms_r);
    if(retval != MPI_SUCCESS) {
        rmg_error_handler (__FILE__, __LINE__, "Error in MPI_Alloc_mem.\n");
    }
    TradeImages::yzpsps_r = TradeImages::yzpsms_r + MAX_IMG2 * THREADS_PER_NODE * TRADE_GRID_EDGES;
    TradeImages::yzmsms_r = TradeImages::yzpsps_r + MAX_IMG2 * THREADS_PER_NODE * TRADE_GRID_EDGES;
    TradeImages::yzmsps_r = TradeImages::yzmsms_r + MAX_IMG2 * THREADS_PER_NODE * TRADE_GRID_EDGES;
    TradeImages::yzpsms_s = TradeImages::yzmsps_r + MAX_IMG2 * THREADS_PER_NODE * TRADE_GRID_EDGES;
    TradeImages::yzpsps_s = TradeImages::yzpsms_s + MAX_IMG2 * THREADS_PER_NODE * TRADE_GRID_EDGES;
    TradeImages::yzmsms_s = TradeImages::yzpsps_s + MAX_IMG2 * THREADS_PER_NODE * TRADE_GRID_EDGES;
    TradeImages::yzmsps_s = TradeImages::yzmsms_s + MAX_IMG2 * THREADS_PER_NODE * TRADE_GRID_EDGES;

    retval = MPI_Alloc_mem(sizeof(complex<double>) * 8 * MAX_IMG2 * THREADS_PER_NODE * TRADE_GRID_EDGES , MPI_INFO_NULL, &TradeImages::xypsms_r);
    if(retval != MPI_SUCCESS) {
        rmg_error_handler (__FILE__, __LINE__, "Error in MPI_Alloc_mem.\n");
    }
    TradeImages::xypsps_r = TradeImages::xypsms_r + MAX_IMG2 * THREADS_PER_NODE * TRADE_GRID_EDGES;
    TradeImages::xymsms_r = TradeImages::xypsps_r + MAX_IMG2 * THREADS_PER_NODE * TRADE_GRID_EDGES;
    TradeImages::xymsps_r = TradeImages::xymsms_r + MAX_IMG2 * THREADS_PER_NODE * TRADE_GRID_EDGES;
    TradeImages::xypsms_s = TradeImages::xymsps_r + MAX_IMG2 * THREADS_PER_NODE * TRADE_GRID_EDGES;
    TradeImages::xypsps_s = TradeImages::xypsms_s + MAX_IMG2 * THREADS_PER_NODE * TRADE_GRID_EDGES;
    TradeImages::xymsms_s = TradeImages::xypsps_s + MAX_IMG2 * THREADS_PER_NODE * TRADE_GRID_EDGES;
    TradeImages::xymsps_s = TradeImages::xymsms_s + MAX_IMG2 * THREADS_PER_NODE * TRADE_GRID_EDGES;

    retval = MPI_Alloc_mem(sizeof(complex<double>) * 8 * MAX_IMG2 * THREADS_PER_NODE * TRADE_GRID_EDGES , MPI_INFO_NULL, &TradeImages::xzpsms_r);
    if(retval != MPI_SUCCESS) {
        rmg_error_handler (__FILE__, __LINE__, "Error in MPI_Alloc_mem.\n");
    }
    TradeImages::xzpsps_r = TradeImages::xzpsms_r + MAX_IMG2 * THREADS_PER_NODE * TRADE_GRID_EDGES;
    TradeImages::xzmsms_r = TradeImages::xzpsps_r + MAX_IMG2 * THREADS_PER_NODE * TRADE_GRID_EDGES;
    TradeImages::xzmsps_r = TradeImages::xzmsms_r + MAX_IMG2 * THREADS_PER_NODE * TRADE_GRID_EDGES;
    TradeImages::xzpsms_s = TradeImages::xzmsps_r + MAX_IMG2 * THREADS_PER_NODE * TRADE_GRID_EDGES;
    TradeImages::xzpsps_s = TradeImages::xzpsms_s + MAX_IMG2 * THREADS_PER_NODE * TRADE_GRID_EDGES;
    TradeImages::xzmsms_s = TradeImages::xzpsps_s + MAX_IMG2 * THREADS_PER_NODE * TRADE_GRID_EDGES;
    TradeImages::xzmsps_s = TradeImages::xzmsms_s + MAX_IMG2 * THREADS_PER_NODE * TRADE_GRID_EDGES;

    retval = MPI_Alloc_mem(sizeof(complex<double>) * 8 * MAX_IMG3 * THREADS_PER_NODE , MPI_INFO_NULL, &TradeImages::m0_r);
    if(retval != MPI_SUCCESS) {
        rmg_error_handler (__FILE__, __LINE__, "Error in MPI_Alloc_mem.\n");
    }

    retval = MPI_Alloc_mem(sizeof(complex<double>) * 8 * MAX_IMG3 * THREADS_PER_NODE , MPI_INFO_NULL, &TradeImages::m0_s);
    if(retval != MPI_SUCCESS) {
        rmg_error_handler (__FILE__, __LINE__, "Error in MPI_Alloc_mem.\n");
    }
}

template <typename RmgType>
void TradeImages::trade_imagesx_async (RmgType * f, RmgType * w, int dimx, int dimy, int dimz, int images)
{
    int ix, iy, iz, ix1, iy1, iz1, incx, incy, incx0, incy0, index, tim, ione = 1;
    int ixs, iys, ixs2, iys2, c1, c2, c3, idx, idx1, img3;
    int xlen, ylen, zlen, yzlen, xylen, xzlen;
    int tid=0, corner_node_stride, node_idx, retval;

    RmgType *frdx1_f, *frdx2_f, *frdy1_f, *frdy2_f, *frdz1_f, *frdz2_f;
    RmgType *frdx1n_f, *frdx2n_f, *frdy1n_f, *frdy2n_f, *frdz1n_f, *frdz2n_f;
    RmgType *yzpsms_r_f, *yzpsps_r_f, *yzmsms_r_f, *yzmsps_r_f;
    RmgType *yzpsms_s_f, *yzpsps_s_f, *yzmsms_s_f, *yzmsps_s_f;
    RmgType *xzpsms_r_f, *xzpsps_r_f, *xzmsms_r_f, *xzmsps_r_f;
    RmgType *xzpsms_s_f, *xzpsps_s_f, *xzmsms_s_f, *xzmsps_s_f;
    RmgType *xypsms_r_f, *xypsps_r_f, *xymsms_r_f, *xymsps_r_f;
    RmgType *xypsms_s_f, *xypsps_s_f, *xymsms_s_f, *xymsps_s_f;
    RmgType *m0_s_f, *m0_r_f;

    frdx1_f = (RmgType *)TradeImages::frdx1;
    frdx2_f = (RmgType *)TradeImages::frdx2;
    frdy1_f = (RmgType *)TradeImages::frdy1;
    frdy2_f = (RmgType *)TradeImages::frdy2;
    frdz1_f = (RmgType *)TradeImages::frdz1;
    frdz2_f = (RmgType *)TradeImages::frdz2;
    frdx1n_f = (RmgType *)TradeImages::frdx1n;
    frdx2n_f = (RmgType *)TradeImages::frdx2n;
    frdy1n_f = (RmgType *)TradeImages::frdy1n;
    frdy2n_f = (RmgType *)TradeImages::frdy2n;
    frdz1n_f = (RmgType *)TradeImages::frdz1n;
    frdz2n_f = (RmgType *)TradeImages::frdz2n;
    yzpsms_r_f = (RmgType *)TradeImages::yzpsms_r;
    yzpsps_r_f = (RmgType *)TradeImages::yzpsps_r;
    yzmsms_r_f = (RmgType *)TradeImages::yzmsms_r;
    yzmsps_r_f = (RmgType *)TradeImages::yzmsps_r;
    yzpsms_s_f = (RmgType *)TradeImages::yzpsms_s;
    yzpsps_s_f = (RmgType *)TradeImages::yzpsps_s;
    yzmsms_s_f = (RmgType *)TradeImages::yzmsms_s;
    yzmsps_s_f = (RmgType *)TradeImages::yzmsps_s;
    xzpsms_r_f = (RmgType *)TradeImages::xzpsms_r;
    xzpsps_r_f = (RmgType *)TradeImages::xzpsps_r;
    xzmsms_r_f = (RmgType *)TradeImages::xzmsms_r;
    xzmsps_r_f = (RmgType *)TradeImages::xzmsps_r;
    xzpsms_s_f = (RmgType *)TradeImages::xzpsms_s;
    xzpsps_s_f = (RmgType *)TradeImages::xzpsps_s;
    xzmsms_s_f = (RmgType *)TradeImages::xzmsms_s;
    xzmsps_s_f = (RmgType *)TradeImages::xzmsps_s;
    xypsms_r_f = (RmgType *)TradeImages::xypsms_r;
    xypsps_r_f = (RmgType *)TradeImages::xypsps_r;
    xymsms_r_f = (RmgType *)TradeImages::xymsms_r;
    xymsps_r_f = (RmgType *)TradeImages::xymsps_r;
    xypsms_s_f = (RmgType *)TradeImages::xypsms_s;
    xypsps_s_f = (RmgType *)TradeImages::xypsps_s;
    xymsms_s_f = (RmgType *)TradeImages::xymsms_s;
    xymsps_s_f = (RmgType *)TradeImages::xymsps_s;
    m0_s_f = (RmgType *)TradeImages::m0_s;
    m0_r_f = (RmgType *)TradeImages::m0_r;

    int ACTIVE_THREADS = 1, THREADS_PER_NODE;
    BaseThread T(0);

    THREADS_PER_NODE = T.get_threads_per_node();


    if(images > MAX_TRADE_IMAGES) {
       rmg_error_handler (__FILE__, __LINE__, "Images count too high in trade_imagesx_async. Modify and recompile may be required.\n");
    }

    tid = T.get_thread_tid();
    if(tid < 0) tid = 0;
    if(T.is_loop_over_states()) ACTIVE_THREADS = THREADS_PER_NODE;

    corner_node_stride = THREADS_PER_NODE * MAX_IMG3; 
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
      TradeImages::RMG_MPI_trade( &m0_r_f[0*corner_node_stride], ACTIVE_THREADS * img3, RMG_MPI_IRECV, -1, -1, -1, transition_get_grid_comm(), 18, &TradeImages::rreqs[18]);
      TradeImages::RMG_MPI_trade( &m0_r_f[1*corner_node_stride], ACTIVE_THREADS * img3, RMG_MPI_IRECV, -1, -1, 1, transition_get_grid_comm(), 19, &TradeImages::rreqs[19]);
      TradeImages::RMG_MPI_trade( &m0_r_f[2*corner_node_stride], ACTIVE_THREADS * img3, RMG_MPI_IRECV, -1, 1, -1, transition_get_grid_comm(), 20, &TradeImages::rreqs[20]);
      TradeImages::RMG_MPI_trade( &m0_r_f[3*corner_node_stride], ACTIVE_THREADS * img3, RMG_MPI_IRECV, -1, 1, 1, transition_get_grid_comm(), 21, &TradeImages::rreqs[21]);
      TradeImages::RMG_MPI_trade( &m0_r_f[4*corner_node_stride], ACTIVE_THREADS * img3, RMG_MPI_IRECV, 1, -1, -1, transition_get_grid_comm(), 22, &TradeImages::rreqs[22]);
      TradeImages::RMG_MPI_trade( &m0_r_f[5*corner_node_stride], ACTIVE_THREADS * img3, RMG_MPI_IRECV, 1, -1, 1, transition_get_grid_comm(), 23, &TradeImages::rreqs[23]);
      TradeImages::RMG_MPI_trade( &m0_r_f[6*corner_node_stride], ACTIVE_THREADS * img3, RMG_MPI_IRECV, 1, 1, -1, transition_get_grid_comm(), 24, &TradeImages::rreqs[24]);
      TradeImages::RMG_MPI_trade( &m0_r_f[7*corner_node_stride], ACTIVE_THREADS * img3, RMG_MPI_IRECV, 1, 1, 1, transition_get_grid_comm(), 25, &TradeImages::rreqs[25]);

        // The z planes
      TradeImages::RMG_MPI_trade(frdz2n_f, ACTIVE_THREADS * zlen, RMG_MPI_IRECV, 0, 0, 1, transition_get_grid_comm(), 0, &TradeImages::rreqs[0]);
      TradeImages::RMG_MPI_trade(frdz1n_f, ACTIVE_THREADS * zlen, RMG_MPI_IRECV, 0, 0, -1, transition_get_grid_comm(), 1, &TradeImages::rreqs[1]);

        // The y planes
      TradeImages::RMG_MPI_trade(frdy2n_f, ACTIVE_THREADS * ylen, RMG_MPI_IRECV, 0, 1, 0, transition_get_grid_comm(), 2, &TradeImages::rreqs[2]);
      TradeImages::RMG_MPI_trade(frdy1n_f, ACTIVE_THREADS * ylen, RMG_MPI_IRECV, 0, -1, 0, transition_get_grid_comm(), 3, &TradeImages::rreqs[3]);

        // The x planes
      TradeImages::RMG_MPI_trade(frdx2n_f, ACTIVE_THREADS * xlen, RMG_MPI_IRECV, 1, 0, 0, transition_get_grid_comm(), 4, &TradeImages::rreqs[4]);
      TradeImages::RMG_MPI_trade(frdx1n_f, ACTIVE_THREADS * xlen, RMG_MPI_IRECV, -1, 0, 0, transition_get_grid_comm(), 5, &TradeImages::rreqs[5]);

        // And the yz-plane edges
      TradeImages::RMG_MPI_trade(yzpsms_r_f, ACTIVE_THREADS * yzlen, RMG_MPI_IRECV, 0, 1, -1, transition_get_grid_comm(), 6, &TradeImages::rreqs[6]);
      TradeImages::RMG_MPI_trade(yzpsps_r_f, ACTIVE_THREADS * yzlen, RMG_MPI_IRECV, 0, 1, 1, transition_get_grid_comm(), 7, &TradeImages::rreqs[7]);
      TradeImages::RMG_MPI_trade(yzmsms_r_f, ACTIVE_THREADS * yzlen, RMG_MPI_IRECV, 0, -1, -1, transition_get_grid_comm(), 8, &TradeImages::rreqs[8]);
      TradeImages::RMG_MPI_trade(yzmsps_r_f, ACTIVE_THREADS * yzlen, RMG_MPI_IRECV, 0, -1, 1, transition_get_grid_comm(), 9, &TradeImages::rreqs[9]);

        // And the xy-plane edges
      TradeImages::RMG_MPI_trade(xypsms_r_f, ACTIVE_THREADS * xylen, RMG_MPI_IRECV, 1, -1, 0, transition_get_grid_comm(), 10, &TradeImages::rreqs[10]);
      TradeImages::RMG_MPI_trade(xypsps_r_f, ACTIVE_THREADS * xylen, RMG_MPI_IRECV, 1, 1, 0, transition_get_grid_comm(), 11, &TradeImages::rreqs[11]);
      TradeImages::RMG_MPI_trade(xymsms_r_f, ACTIVE_THREADS * xylen, RMG_MPI_IRECV, -1, -1, 0, transition_get_grid_comm(), 12, &TradeImages::rreqs[12]);
      TradeImages::RMG_MPI_trade(xymsps_r_f, ACTIVE_THREADS * xylen, RMG_MPI_IRECV, -1, 1, 0, transition_get_grid_comm(), 13, &TradeImages::rreqs[13]);

        // And the xz-plane edges
      TradeImages::RMG_MPI_trade(xzpsms_r_f, ACTIVE_THREADS * xzlen, RMG_MPI_IRECV, 1, 0, -1, transition_get_grid_comm(), 14, &TradeImages::rreqs[14]);
      TradeImages::RMG_MPI_trade(xzpsps_r_f, ACTIVE_THREADS * xzlen, RMG_MPI_IRECV, 1, 0, 1, transition_get_grid_comm(), 15, &TradeImages::rreqs[15]);
      TradeImages::RMG_MPI_trade(xzmsms_r_f, ACTIVE_THREADS * xzlen, RMG_MPI_IRECV, -1, 0, -1, transition_get_grid_comm(), 16, &TradeImages::rreqs[16]);
      TradeImages::RMG_MPI_trade(xzmsps_r_f, ACTIVE_THREADS * xzlen, RMG_MPI_IRECV, -1, 0, 1, transition_get_grid_comm(), 17, &TradeImages::rreqs[17]);

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
    T.thread_barrier_wait();
    if(tid == 0) {

        // Corners
      TradeImages::RMG_MPI_trade( &m0_s_f[0*corner_node_stride], ACTIVE_THREADS * img3, RMG_MPI_ISEND, 1, 1, 1, transition_get_grid_comm(), 18, &TradeImages::sreqs[18]);
      TradeImages::RMG_MPI_trade( &m0_s_f[1*corner_node_stride], ACTIVE_THREADS * img3, RMG_MPI_ISEND, 1, 1, -1, transition_get_grid_comm(), 19, &TradeImages::sreqs[19]);
      TradeImages::RMG_MPI_trade( &m0_s_f[2*corner_node_stride], ACTIVE_THREADS * img3, RMG_MPI_ISEND, 1, -1, 1, transition_get_grid_comm(), 20, &TradeImages::sreqs[20]);
      TradeImages::RMG_MPI_trade( &m0_s_f[3*corner_node_stride], ACTIVE_THREADS * img3, RMG_MPI_ISEND, 1, -1, -1, transition_get_grid_comm(), 21, &TradeImages::sreqs[21]);
      TradeImages::RMG_MPI_trade( &m0_s_f[4*corner_node_stride], ACTIVE_THREADS * img3, RMG_MPI_ISEND, -1, 1, 1, transition_get_grid_comm(), 22, &TradeImages::sreqs[22]);
      TradeImages::RMG_MPI_trade( &m0_s_f[5*corner_node_stride], ACTIVE_THREADS * img3, RMG_MPI_ISEND, -1, 1, -1, transition_get_grid_comm(), 23, &TradeImages::sreqs[23]);
      TradeImages::RMG_MPI_trade( &m0_s_f[6*corner_node_stride], ACTIVE_THREADS * img3, RMG_MPI_ISEND, -1, -1, 1, transition_get_grid_comm(), 24, &TradeImages::sreqs[24]);
      TradeImages::RMG_MPI_trade( &m0_s_f[7*corner_node_stride], ACTIVE_THREADS * img3, RMG_MPI_ISEND, -1, -1, -1, transition_get_grid_comm(), 25, &TradeImages::sreqs[25]);

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


    T.thread_barrier_wait();
    if(tid == 0) {

        // Send z planes
      TradeImages::RMG_MPI_trade(frdz1_f, ACTIVE_THREADS * zlen, RMG_MPI_ISEND, 0, 0, -1, transition_get_grid_comm(), 0, &TradeImages::sreqs[0]);
      TradeImages::RMG_MPI_trade(frdz2_f, ACTIVE_THREADS * zlen, RMG_MPI_ISEND, 0, 0, 1, transition_get_grid_comm(), 1, &TradeImages::sreqs[1]);

        // Send the north and south planes
      TradeImages::RMG_MPI_trade(frdy1_f, ACTIVE_THREADS * ylen, RMG_MPI_ISEND, 0, -1, 0, transition_get_grid_comm(), 2, &TradeImages::sreqs[2]);
      TradeImages::RMG_MPI_trade(frdy2_f, ACTIVE_THREADS * ylen, RMG_MPI_ISEND, 0, 1, 0, transition_get_grid_comm(), 3, &TradeImages::sreqs[3]);

        // Send the east and west planes
      TradeImages::RMG_MPI_trade(frdx1_f, ACTIVE_THREADS * xlen, RMG_MPI_ISEND, -1, 0, 0, transition_get_grid_comm(), 4, &TradeImages::sreqs[4]);
      TradeImages::RMG_MPI_trade(frdx2_f, ACTIVE_THREADS * xlen, RMG_MPI_ISEND, 1, 0, 0, transition_get_grid_comm(), 5, &TradeImages::sreqs[5]);

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
    T.thread_barrier_wait();
    if(tid == 0) {

        // Send yz-plane edges
      TradeImages::RMG_MPI_trade(yzpsms_s_f, ACTIVE_THREADS * yzlen, RMG_MPI_ISEND, 0, -1, 1, transition_get_grid_comm(), 6, &TradeImages::sreqs[6]);
      TradeImages::RMG_MPI_trade(yzpsps_s_f, ACTIVE_THREADS * yzlen, RMG_MPI_ISEND, 0, -1, -1, transition_get_grid_comm(), 7, &TradeImages::sreqs[7]);
      TradeImages::RMG_MPI_trade(yzmsms_s_f, ACTIVE_THREADS * yzlen, RMG_MPI_ISEND, 0, 1, 1, transition_get_grid_comm(), 8, &TradeImages::sreqs[8]);
      TradeImages::RMG_MPI_trade(yzmsps_s_f, ACTIVE_THREADS * yzlen, RMG_MPI_ISEND, 0, 1, -1, transition_get_grid_comm(), 9, &TradeImages::sreqs[9]);

        // Send xy-plane edges
      TradeImages::RMG_MPI_trade(xypsms_s_f, ACTIVE_THREADS * xylen, RMG_MPI_ISEND, -1, 1, 0, transition_get_grid_comm(), 10, &TradeImages::sreqs[10]);
      TradeImages::RMG_MPI_trade(xypsps_s_f, ACTIVE_THREADS * xylen, RMG_MPI_ISEND, -1, -1, 0, transition_get_grid_comm(), 11, &TradeImages::sreqs[11]);
      TradeImages::RMG_MPI_trade(xymsms_s_f, ACTIVE_THREADS * xylen, RMG_MPI_ISEND, 1, 1, 0, transition_get_grid_comm(), 12, &TradeImages::sreqs[12]);
      TradeImages::RMG_MPI_trade(xymsps_s_f, ACTIVE_THREADS * xylen, RMG_MPI_ISEND, 1, -1, 0, transition_get_grid_comm(), 13, &TradeImages::sreqs[13]);

        // Send xz-plane edges
      TradeImages::RMG_MPI_trade(xzpsms_s_f, ACTIVE_THREADS * xzlen, RMG_MPI_ISEND, -1, 0, 1, transition_get_grid_comm(), 14, &TradeImages::sreqs[14]);
      TradeImages::RMG_MPI_trade(xzpsps_s_f, ACTIVE_THREADS * xzlen, RMG_MPI_ISEND, -1, 0, -1, transition_get_grid_comm(), 15, &TradeImages::sreqs[15]);
      TradeImages::RMG_MPI_trade(xzmsms_s_f, ACTIVE_THREADS * xzlen, RMG_MPI_ISEND, 1, 0, 1, transition_get_grid_comm(), 16, &TradeImages::sreqs[16]);
      TradeImages::RMG_MPI_trade(xzmsps_s_f, ACTIVE_THREADS * xzlen, RMG_MPI_ISEND, 1, 0, -1, transition_get_grid_comm(), 17, &TradeImages::sreqs[17]);

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

            QMD_copy (dimz, &f[iys], ione, &w[iys2 + images], ione);

        }                       /* end for */

    }                           /* end for */


    // Wait for all the recvs to finish
    if(tid == 0) {
        retval = MPI_Waitall(26, TradeImages::rreqs, MPI_STATUSES_IGNORE);
        if(retval != MPI_SUCCESS) rmg_error_handler (__FILE__, __LINE__, "Error in MPI_Waitall.\n");
    }
    T.thread_barrier_wait();


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
    T.thread_barrier_wait();

} // end trade_imagesx_async



// Asynchronous image trades for central finite difference operators
template <typename RmgType>
void TradeImages::trade_imagesx_central_async (RmgType * f, RmgType * w, int dimx, int dimy, int dimz, int images)
{
    int ix, iy, iz, incx, incy, incx0, incy0, index, tim, ione = 1;
    int ixs, iys, ixs2, iys2, c1, idx;
    int xlen, ylen, zlen;
    int tid=0, retval;

    RmgType *frdx1_f, *frdx2_f, *frdy1_f, *frdy2_f, *frdz1_f, *frdz2_f;
    RmgType *frdx1n_f, *frdx2n_f, *frdy1n_f, *frdy2n_f, *frdz1n_f, *frdz2n_f;

    frdx1_f = (RmgType *)TradeImages::frdx1;
    frdx2_f = (RmgType *)TradeImages::frdx2;
    frdy1_f = (RmgType *)TradeImages::frdy1;
    frdy2_f = (RmgType *)TradeImages::frdy2;
    frdz1_f = (RmgType *)TradeImages::frdz1;
    frdz2_f = (RmgType *)TradeImages::frdz2;
    frdx1n_f = (RmgType *)TradeImages::frdx1n;
    frdx2n_f = (RmgType *)TradeImages::frdx2n;
    frdy1n_f = (RmgType *)TradeImages::frdy1n;
    frdy2n_f = (RmgType *)TradeImages::frdy2n;
    frdz1n_f = (RmgType *)TradeImages::frdz1n;
    frdz2n_f = (RmgType *)TradeImages::frdz2n;

    int ACTIVE_THREADS = 1, THREADS_PER_NODE;

    BaseThread T(0);
    THREADS_PER_NODE = T.get_threads_per_node();


    if(images > MAX_TRADE_IMAGES) {
       rmg_error_handler (__FILE__, __LINE__, "Images count too high in trade_imagesx_async. Modify and recompile may be required.\n");
    }

    tid = T.get_thread_tid();
    if(tid < 0) tid = 0;
    if(T.is_loop_over_states()) ACTIVE_THREADS = THREADS_PER_NODE;

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
      TradeImages::RMG_MPI_trade(frdz2n_f, ACTIVE_THREADS * zlen, RMG_MPI_IRECV, 0, 0, 1, transition_get_grid_comm(), 0, &TradeImages::rreqs[0]);
      TradeImages::RMG_MPI_trade(frdz1n_f, ACTIVE_THREADS * zlen, RMG_MPI_IRECV, 0, 0, -1, transition_get_grid_comm(), 1, &TradeImages::rreqs[1]);

        // The y planes
      TradeImages::RMG_MPI_trade(frdy2n_f, ACTIVE_THREADS * ylen, RMG_MPI_IRECV, 0, 1, 0, transition_get_grid_comm(), 2, &TradeImages::rreqs[2]);
      TradeImages::RMG_MPI_trade(frdy1n_f, ACTIVE_THREADS * ylen, RMG_MPI_IRECV, 0, -1, 0, transition_get_grid_comm(), 3, &TradeImages::rreqs[3]);

        // The x planes
      TradeImages::RMG_MPI_trade(frdx2n_f, ACTIVE_THREADS * xlen, RMG_MPI_IRECV, 1, 0, 0, transition_get_grid_comm(), 4, &TradeImages::rreqs[4]);
      TradeImages::RMG_MPI_trade(frdx1n_f, ACTIVE_THREADS * xlen, RMG_MPI_IRECV, -1, 0, 0, transition_get_grid_comm(), 5, &TradeImages::rreqs[5]);

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


    T.thread_barrier_wait();
    if(tid == 0) {

        // Send z planes
      TradeImages::RMG_MPI_trade(frdz1_f, ACTIVE_THREADS * zlen, RMG_MPI_ISEND, 0, 0, -1, transition_get_grid_comm(), 0, &TradeImages::sreqs[0]);
      TradeImages::RMG_MPI_trade(frdz2_f, ACTIVE_THREADS * zlen, RMG_MPI_ISEND, 0, 0, 1, transition_get_grid_comm(), 1, &TradeImages::sreqs[1]);

        // Send the north and south planes
      TradeImages::RMG_MPI_trade(frdy1_f, ACTIVE_THREADS * ylen, RMG_MPI_ISEND, 0, -1, 0, transition_get_grid_comm(), 2, &TradeImages::sreqs[2]);
      TradeImages::RMG_MPI_trade(frdy2_f, ACTIVE_THREADS * ylen, RMG_MPI_ISEND, 0, 1, 0, transition_get_grid_comm(), 3, &TradeImages::sreqs[3]);

        // Send the east and west planes
      TradeImages::RMG_MPI_trade(frdx1_f, ACTIVE_THREADS * xlen, RMG_MPI_ISEND, -1, 0, 0, transition_get_grid_comm(), 4, &TradeImages::sreqs[4]);
      TradeImages::RMG_MPI_trade(frdx2_f, ACTIVE_THREADS * xlen, RMG_MPI_ISEND, 1, 0, 0, transition_get_grid_comm(), 5, &TradeImages::sreqs[5]);

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

            QMD_copy (dimz, &f[iys], ione, &w[iys2 + images], ione);

        }                       /* end for */

    }                           /* end for */


    // Wait for all the recvs to finish
    if(tid == 0) {
        retval = MPI_Waitall(6, TradeImages::rreqs, MPI_STATUSES_IGNORE);
        if(retval != MPI_SUCCESS) rmg_error_handler (__FILE__, __LINE__, "Error in MPI_Waitall.\n");
    }
    T.thread_barrier_wait();


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
    T.thread_barrier_wait();

} // end trade_imagesx_central_async



// Used by coarse level multigrid routines which assume that the central data is already present
// and that the f array is sized [dimx+2][dimy+2][dimz+2]
//
template <typename RmgType>
void TradeImages::trade_images1_async (RmgType * f, int dimx, int dimy, int dimz)
{
    int ix, iy, iz, incx, incy, index;
    int ixs2, iys2, c1, c2, c3, idx, idx1;
    int xlen, ylen, zlen, yzlen, xylen, xzlen;
    int tid=0, corner_node_stride, node_idx, retval;

    RmgType *frdx1_f, *frdx2_f, *frdy1_f, *frdy2_f, *frdz1_f, *frdz2_f;
    RmgType *frdx1n_f, *frdx2n_f, *frdy1n_f, *frdy2n_f, *frdz1n_f, *frdz2n_f;
    RmgType *yzpsms_r_f, *yzpsps_r_f, *yzmsms_r_f, *yzmsps_r_f;
    RmgType *yzpsms_s_f, *yzpsps_s_f, *yzmsms_s_f, *yzmsps_s_f;
    RmgType *xzpsms_r_f, *xzpsps_r_f, *xzmsms_r_f, *xzmsps_r_f;
    RmgType *xzpsms_s_f, *xzpsps_s_f, *xzmsms_s_f, *xzmsps_s_f;
    RmgType *xypsms_r_f, *xypsps_r_f, *xymsms_r_f, *xymsps_r_f;
    RmgType *xypsms_s_f, *xypsps_s_f, *xymsms_s_f, *xymsps_s_f;
    RmgType *m0_s_f, *m0_r_f;

    frdx1_f = (RmgType *)TradeImages::frdx1;
    frdx2_f = (RmgType *)TradeImages::frdx2;
    frdy1_f = (RmgType *)TradeImages::frdy1;
    frdy2_f = (RmgType *)TradeImages::frdy2;
    frdz1_f = (RmgType *)TradeImages::frdz1;
    frdz2_f = (RmgType *)TradeImages::frdz2;
    frdx1n_f = (RmgType *)TradeImages::frdx1n;
    frdx2n_f = (RmgType *)TradeImages::frdx2n;
    frdy1n_f = (RmgType *)TradeImages::frdy1n;
    frdy2n_f = (RmgType *)TradeImages::frdy2n;
    frdz1n_f = (RmgType *)TradeImages::frdz1n;
    frdz2n_f = (RmgType *)TradeImages::frdz2n;
    yzpsms_r_f = (RmgType *)TradeImages::yzpsms_r;
    yzpsps_r_f = (RmgType *)TradeImages::yzpsps_r;
    yzmsms_r_f = (RmgType *)TradeImages::yzmsms_r;
    yzmsps_r_f = (RmgType *)TradeImages::yzmsps_r;
    yzpsms_s_f = (RmgType *)TradeImages::yzpsms_s;
    yzpsps_s_f = (RmgType *)TradeImages::yzpsps_s;
    yzmsms_s_f = (RmgType *)TradeImages::yzmsms_s;
    yzmsps_s_f = (RmgType *)TradeImages::yzmsps_s;
    xzpsms_r_f = (RmgType *)TradeImages::xzpsms_r;
    xzpsps_r_f = (RmgType *)TradeImages::xzpsps_r;
    xzmsms_r_f = (RmgType *)TradeImages::xzmsms_r;
    xzmsps_r_f = (RmgType *)TradeImages::xzmsps_r;
    xzpsms_s_f = (RmgType *)TradeImages::xzpsms_s;
    xzpsps_s_f = (RmgType *)TradeImages::xzpsps_s;
    xzmsms_s_f = (RmgType *)TradeImages::xzmsms_s;
    xzmsps_s_f = (RmgType *)TradeImages::xzmsps_s;
    xypsms_r_f = (RmgType *)TradeImages::xypsms_r;
    xypsps_r_f = (RmgType *)TradeImages::xypsps_r;
    xymsms_r_f = (RmgType *)TradeImages::xymsms_r;
    xymsps_r_f = (RmgType *)TradeImages::xymsps_r;
    xypsms_s_f = (RmgType *)TradeImages::xypsms_s;
    xypsps_s_f = (RmgType *)TradeImages::xypsps_s;
    xymsms_s_f = (RmgType *)TradeImages::xymsms_s;
    xymsps_s_f = (RmgType *)TradeImages::xymsps_s;
    m0_s_f = (RmgType *)TradeImages::m0_s;
    m0_r_f = (RmgType *)TradeImages::m0_r;

    int ACTIVE_THREADS = 1, THREADS_PER_NODE;
    BaseThread T(0);
    THREADS_PER_NODE = T.get_threads_per_node();

    tid = T.get_thread_tid();
    if(tid < 0) tid = 0;
    if(T.is_loop_over_states()) ACTIVE_THREADS = THREADS_PER_NODE;

    corner_node_stride = THREADS_PER_NODE * MAX_IMG3;

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
      TradeImages::RMG_MPI_trade( &m0_r_f[0*corner_node_stride], ACTIVE_THREADS, RMG_MPI_IRECV, -1, -1, -1, transition_get_grid_comm(), 18, &TradeImages::rreqs[18]);
      TradeImages::RMG_MPI_trade( &m0_r_f[1*corner_node_stride], ACTIVE_THREADS, RMG_MPI_IRECV, -1, -1, 1, transition_get_grid_comm(), 19, &TradeImages::rreqs[19]);
      TradeImages::RMG_MPI_trade( &m0_r_f[2*corner_node_stride], ACTIVE_THREADS, RMG_MPI_IRECV, -1, 1, -1, transition_get_grid_comm(), 20, &TradeImages::rreqs[20]);
      TradeImages::RMG_MPI_trade( &m0_r_f[3*corner_node_stride], ACTIVE_THREADS, RMG_MPI_IRECV, -1, 1, 1, transition_get_grid_comm(), 21, &TradeImages::rreqs[21]);
      TradeImages::RMG_MPI_trade( &m0_r_f[4*corner_node_stride], ACTIVE_THREADS, RMG_MPI_IRECV, 1, -1, -1, transition_get_grid_comm(), 22, &TradeImages::rreqs[22]);
      TradeImages::RMG_MPI_trade( &m0_r_f[5*corner_node_stride], ACTIVE_THREADS, RMG_MPI_IRECV, 1, -1, 1, transition_get_grid_comm(), 23, &TradeImages::rreqs[23]);
      TradeImages::RMG_MPI_trade( &m0_r_f[6*corner_node_stride], ACTIVE_THREADS, RMG_MPI_IRECV, 1, 1, -1, transition_get_grid_comm(), 24, &TradeImages::rreqs[24]);
      TradeImages::RMG_MPI_trade( &m0_r_f[7*corner_node_stride], ACTIVE_THREADS, RMG_MPI_IRECV, 1, 1, 1, transition_get_grid_comm(), 25, &TradeImages::rreqs[25]);

        // The z planes
      TradeImages::RMG_MPI_trade(frdz2n_f, ACTIVE_THREADS * zlen, RMG_MPI_IRECV, 0, 0, 1, transition_get_grid_comm(), 0, &TradeImages::rreqs[0]);
      TradeImages::RMG_MPI_trade(frdz1n_f, ACTIVE_THREADS * zlen, RMG_MPI_IRECV, 0, 0, -1, transition_get_grid_comm(), 1, &TradeImages::rreqs[1]);

        // The y planes
      TradeImages::RMG_MPI_trade(frdy2n_f, ACTIVE_THREADS * ylen, RMG_MPI_IRECV, 0, 1, 0, transition_get_grid_comm(), 2, &TradeImages::rreqs[2]);
      TradeImages::RMG_MPI_trade(frdy1n_f, ACTIVE_THREADS * ylen, RMG_MPI_IRECV, 0, -1, 0, transition_get_grid_comm(), 3, &TradeImages::rreqs[3]);

        // The x planes
      TradeImages::RMG_MPI_trade(frdx2n_f, ACTIVE_THREADS * xlen, RMG_MPI_IRECV, 1, 0, 0, transition_get_grid_comm(), 4, &TradeImages::rreqs[4]);
      TradeImages::RMG_MPI_trade(frdx1n_f, ACTIVE_THREADS * xlen, RMG_MPI_IRECV, -1, 0, 0, transition_get_grid_comm(), 5, &TradeImages::rreqs[5]);

        // And the yz-plane edges
      TradeImages::RMG_MPI_trade(yzpsms_r_f, ACTIVE_THREADS * yzlen, RMG_MPI_IRECV, 0, 1, -1, transition_get_grid_comm(), 6, &TradeImages::rreqs[6]);
      TradeImages::RMG_MPI_trade(yzpsps_r_f, ACTIVE_THREADS * yzlen, RMG_MPI_IRECV, 0, 1, 1, transition_get_grid_comm(), 7, &TradeImages::rreqs[7]);
      TradeImages::RMG_MPI_trade(yzmsms_r_f, ACTIVE_THREADS * yzlen, RMG_MPI_IRECV, 0, -1, -1, transition_get_grid_comm(), 8, &TradeImages::rreqs[8]);
      TradeImages::RMG_MPI_trade(yzmsps_r_f, ACTIVE_THREADS * yzlen, RMG_MPI_IRECV, 0, -1, 1, transition_get_grid_comm(), 9, &TradeImages::rreqs[9]);

        // And the xy-plane edges
      TradeImages::RMG_MPI_trade(xypsms_r_f, ACTIVE_THREADS * xylen, RMG_MPI_IRECV, 1, -1, 0, transition_get_grid_comm(), 10, &TradeImages::rreqs[10]);
      TradeImages::RMG_MPI_trade(xypsps_r_f, ACTIVE_THREADS * xylen, RMG_MPI_IRECV, 1, 1, 0, transition_get_grid_comm(), 11, &TradeImages::rreqs[11]);
      TradeImages::RMG_MPI_trade(xymsms_r_f, ACTIVE_THREADS * xylen, RMG_MPI_IRECV, -1, -1, 0, transition_get_grid_comm(), 12, &TradeImages::rreqs[12]);
      TradeImages::RMG_MPI_trade(xymsps_r_f, ACTIVE_THREADS * xylen, RMG_MPI_IRECV, -1, 1, 0, transition_get_grid_comm(), 13, &TradeImages::rreqs[13]);

        // And the xz-plane edges
      TradeImages::RMG_MPI_trade(xzpsms_r_f, ACTIVE_THREADS * xzlen, RMG_MPI_IRECV, 1, 0, -1, transition_get_grid_comm(), 14, &TradeImages::rreqs[14]);
      TradeImages::RMG_MPI_trade(xzpsps_r_f, ACTIVE_THREADS * xzlen, RMG_MPI_IRECV, 1, 0, 1, transition_get_grid_comm(), 15, &TradeImages::rreqs[15]);
      TradeImages::RMG_MPI_trade(xzmsms_r_f, ACTIVE_THREADS * xzlen, RMG_MPI_IRECV, -1, 0, -1, transition_get_grid_comm(), 16, &TradeImages::rreqs[16]);
      TradeImages::RMG_MPI_trade(xzmsps_r_f, ACTIVE_THREADS * xzlen, RMG_MPI_IRECV, -1, 0, 1, transition_get_grid_comm(), 17, &TradeImages::rreqs[17]);

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
    T.thread_barrier_wait();
    if(tid == 0) {

        // Corners
      TradeImages::RMG_MPI_trade( &m0_s_f[0*corner_node_stride], ACTIVE_THREADS, RMG_MPI_ISEND, 1, 1, 1, transition_get_grid_comm(), 18, &TradeImages::sreqs[18]);
      TradeImages::RMG_MPI_trade( &m0_s_f[1*corner_node_stride], ACTIVE_THREADS, RMG_MPI_ISEND, 1, 1, -1, transition_get_grid_comm(), 19, &TradeImages::sreqs[19]);
      TradeImages::RMG_MPI_trade( &m0_s_f[2*corner_node_stride], ACTIVE_THREADS, RMG_MPI_ISEND, 1, -1, 1, transition_get_grid_comm(), 20, &TradeImages::sreqs[20]);
      TradeImages::RMG_MPI_trade( &m0_s_f[3*corner_node_stride], ACTIVE_THREADS, RMG_MPI_ISEND, 1, -1, -1, transition_get_grid_comm(), 21, &TradeImages::sreqs[21]);
      TradeImages::RMG_MPI_trade( &m0_s_f[4*corner_node_stride], ACTIVE_THREADS, RMG_MPI_ISEND, -1, 1, 1, transition_get_grid_comm(), 22, &TradeImages::sreqs[22]);
      TradeImages::RMG_MPI_trade( &m0_s_f[5*corner_node_stride], ACTIVE_THREADS, RMG_MPI_ISEND, -1, 1, -1, transition_get_grid_comm(), 23, &TradeImages::sreqs[23]);
      TradeImages::RMG_MPI_trade( &m0_s_f[6*corner_node_stride], ACTIVE_THREADS, RMG_MPI_ISEND, -1, -1, 1, transition_get_grid_comm(), 24, &TradeImages::sreqs[24]);
      TradeImages::RMG_MPI_trade( &m0_s_f[7*corner_node_stride], ACTIVE_THREADS, RMG_MPI_ISEND, -1, -1, -1, transition_get_grid_comm(), 25, &TradeImages::sreqs[25]);

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


    T.thread_barrier_wait();
    if(tid == 0) {

        // Send z planes
      TradeImages::RMG_MPI_trade(frdz1_f, ACTIVE_THREADS * zlen, RMG_MPI_ISEND, 0, 0, -1, transition_get_grid_comm(), 0, &TradeImages::sreqs[0]);
      TradeImages::RMG_MPI_trade(frdz2_f, ACTIVE_THREADS * zlen, RMG_MPI_ISEND, 0, 0, 1, transition_get_grid_comm(), 1, &TradeImages::sreqs[1]);

        // Send the north and south planes
      TradeImages::RMG_MPI_trade(frdy1_f, ACTIVE_THREADS * ylen, RMG_MPI_ISEND, 0, -1, 0, transition_get_grid_comm(), 2, &TradeImages::sreqs[2]);
      TradeImages::RMG_MPI_trade(frdy2_f, ACTIVE_THREADS * ylen, RMG_MPI_ISEND, 0, 1, 0, transition_get_grid_comm(), 3, &TradeImages::sreqs[3]);

        // Send the east and west planes
      TradeImages::RMG_MPI_trade(frdx1_f, ACTIVE_THREADS * xlen, RMG_MPI_ISEND, -1, 0, 0, transition_get_grid_comm(), 4, &TradeImages::sreqs[4]);
      TradeImages::RMG_MPI_trade(frdx2_f, ACTIVE_THREADS * xlen, RMG_MPI_ISEND, 1, 0, 0, transition_get_grid_comm(), 5, &TradeImages::sreqs[5]);

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
    T.thread_barrier_wait();
    if(tid == 0) {

        // Send yz-plane edges
      TradeImages::RMG_MPI_trade(yzpsms_s_f, ACTIVE_THREADS * yzlen, RMG_MPI_ISEND, 0, -1, 1, transition_get_grid_comm(), 6, &TradeImages::sreqs[6]);
      TradeImages::RMG_MPI_trade(yzpsps_s_f, ACTIVE_THREADS * yzlen, RMG_MPI_ISEND, 0, -1, -1, transition_get_grid_comm(), 7, &TradeImages::sreqs[7]);
      TradeImages::RMG_MPI_trade(yzmsms_s_f, ACTIVE_THREADS * yzlen, RMG_MPI_ISEND, 0, 1, 1, transition_get_grid_comm(), 8, &TradeImages::sreqs[8]);
      TradeImages::RMG_MPI_trade(yzmsps_s_f, ACTIVE_THREADS * yzlen, RMG_MPI_ISEND, 0, 1, -1, transition_get_grid_comm(), 9, &TradeImages::sreqs[9]);

        // Send xy-plane edges
      TradeImages::RMG_MPI_trade(xypsms_s_f, ACTIVE_THREADS * xylen, RMG_MPI_ISEND, -1, 1, 0, transition_get_grid_comm(), 10, &TradeImages::sreqs[10]);
      TradeImages::RMG_MPI_trade(xypsps_s_f, ACTIVE_THREADS * xylen, RMG_MPI_ISEND, -1, -1, 0, transition_get_grid_comm(), 11, &TradeImages::sreqs[11]);
      TradeImages::RMG_MPI_trade(xymsms_s_f, ACTIVE_THREADS * xylen, RMG_MPI_ISEND, 1, 1, 0, transition_get_grid_comm(), 12, &TradeImages::sreqs[12]);
      TradeImages::RMG_MPI_trade(xymsps_s_f, ACTIVE_THREADS * xylen, RMG_MPI_ISEND, 1, -1, 0, transition_get_grid_comm(), 13, &TradeImages::sreqs[13]);

        // Send xz-plane edges
      TradeImages::RMG_MPI_trade(xzpsms_s_f, ACTIVE_THREADS * xzlen, RMG_MPI_ISEND, -1, 0, 1, transition_get_grid_comm(), 14, &TradeImages::sreqs[14]);
      TradeImages::RMG_MPI_trade(xzpsps_s_f, ACTIVE_THREADS * xzlen, RMG_MPI_ISEND, -1, 0, -1, transition_get_grid_comm(), 15, &TradeImages::sreqs[15]);
      TradeImages::RMG_MPI_trade(xzmsms_s_f, ACTIVE_THREADS * xzlen, RMG_MPI_ISEND, 1, 0, 1, transition_get_grid_comm(), 16, &TradeImages::sreqs[16]);
      TradeImages::RMG_MPI_trade(xzmsps_s_f, ACTIVE_THREADS * xzlen, RMG_MPI_ISEND, 1, 0, -1, transition_get_grid_comm(), 17, &TradeImages::sreqs[17]);

    }

    // Wait for all the recvs to finish
    if(tid == 0) {
        retval = MPI_Waitall(26, TradeImages::rreqs, MPI_STATUSES_IGNORE);
        if(retval != MPI_SUCCESS) rmg_error_handler (__FILE__, __LINE__, "Error in MPI_Waitall.\n");
    }
    T.thread_barrier_wait();


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
    T.thread_barrier_wait();

} // end trade_images1_async


// Used by coarse level multigrid routines which assume that the central data is already present
// and that the f array is sized [dimx+2][dimy+2][dimz+2]
//
template <typename RmgType>
void TradeImages::trade_images1_central_async (RmgType * f, int dimx, int dimy, int dimz)
{
    int ix, iy, iz, incx, incy, index;
    int ixs2, iys2, c1, idx;
    int xlen, ylen, zlen;
    int tid=0, retval;
    int ACTIVE_THREADS = 1, THREADS_PER_NODE;

    RmgType *frdx1_f, *frdx2_f, *frdy1_f, *frdy2_f, *frdz1_f, *frdz2_f;
    RmgType *frdx1n_f, *frdx2n_f, *frdy1n_f, *frdy2n_f, *frdz1n_f, *frdz2n_f;

    frdx1_f = (RmgType *)TradeImages::frdx1;
    frdx2_f = (RmgType *)TradeImages::frdx2;
    frdy1_f = (RmgType *)TradeImages::frdy1;
    frdy2_f = (RmgType *)TradeImages::frdy2;
    frdz1_f = (RmgType *)TradeImages::frdz1;
    frdz2_f = (RmgType *)TradeImages::frdz2;
    frdx1n_f = (RmgType *)TradeImages::frdx1n;
    frdx2n_f = (RmgType *)TradeImages::frdx2n;
    frdy1n_f = (RmgType *)TradeImages::frdy1n;
    frdy2n_f = (RmgType *)TradeImages::frdy2n;
    frdz1n_f = (RmgType *)TradeImages::frdz1n;
    frdz2n_f = (RmgType *)TradeImages::frdz2n;

    BaseThread T(0);
    THREADS_PER_NODE = T.get_threads_per_node();

    tid = T.get_thread_tid();
    if(tid < 0) tid = 0;
    if(T.is_loop_over_states()) ACTIVE_THREADS = THREADS_PER_NODE;


    incx = (dimy + 2) * (dimz + 2);
    incy = dimz + 2;

    zlen = dimx * dimy;
    ylen = dimx * dimz;
    xlen = dimy * dimz;

    // Thread 0 posts all of the receives
    if(tid == 0) {

        // The z planes
      TradeImages::RMG_MPI_trade(frdz2n_f, ACTIVE_THREADS * zlen, RMG_MPI_IRECV, 0, 0, 1, transition_get_grid_comm(), 0, &TradeImages::rreqs[0]);
      TradeImages::RMG_MPI_trade(frdz1n_f, ACTIVE_THREADS * zlen, RMG_MPI_IRECV, 0, 0, -1, transition_get_grid_comm(), 1, &TradeImages::rreqs[1]);

        // The y planes
      TradeImages::RMG_MPI_trade(frdy2n_f, ACTIVE_THREADS * ylen, RMG_MPI_IRECV, 0, 1, 0, transition_get_grid_comm(), 2, &TradeImages::rreqs[2]);
      TradeImages::RMG_MPI_trade(frdy1n_f, ACTIVE_THREADS * ylen, RMG_MPI_IRECV, 0, -1, 0, transition_get_grid_comm(), 3, &TradeImages::rreqs[3]);

        // The x planes
      TradeImages::RMG_MPI_trade(frdx2n_f, ACTIVE_THREADS * xlen, RMG_MPI_IRECV, 1, 0, 0, transition_get_grid_comm(), 4, &TradeImages::rreqs[4]);
      TradeImages::RMG_MPI_trade(frdx1n_f, ACTIVE_THREADS * xlen, RMG_MPI_IRECV, -1, 0, 0, transition_get_grid_comm(), 5, &TradeImages::rreqs[5]);

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


    T.thread_barrier_wait();
    if(tid == 0) {

        // Send z planes
      TradeImages::RMG_MPI_trade(frdz1_f, ACTIVE_THREADS * zlen, RMG_MPI_ISEND, 0, 0, -1, transition_get_grid_comm(), 0, &TradeImages::sreqs[0]);
      TradeImages::RMG_MPI_trade(frdz2_f, ACTIVE_THREADS * zlen, RMG_MPI_ISEND, 0, 0, 1, transition_get_grid_comm(), 1, &TradeImages::sreqs[1]);

        // Send the north and south planes
      TradeImages::RMG_MPI_trade(frdy1_f, ACTIVE_THREADS * ylen, RMG_MPI_ISEND, 0, -1, 0, transition_get_grid_comm(), 2, &TradeImages::sreqs[2]);
      TradeImages::RMG_MPI_trade(frdy2_f, ACTIVE_THREADS * ylen, RMG_MPI_ISEND, 0, 1, 0, transition_get_grid_comm(), 3, &TradeImages::sreqs[3]);

        // Send the east and west planes
      TradeImages::RMG_MPI_trade(frdx1_f, ACTIVE_THREADS * xlen, RMG_MPI_ISEND, -1, 0, 0, transition_get_grid_comm(), 4, &TradeImages::sreqs[4]);
      TradeImages::RMG_MPI_trade(frdx2_f, ACTIVE_THREADS * xlen, RMG_MPI_ISEND, 1, 0, 0, transition_get_grid_comm(), 5, &TradeImages::sreqs[5]);

    }

    // Wait for all the recvs to finish
    if(tid == 0) {
        retval = MPI_Waitall(6, TradeImages::rreqs, MPI_STATUSES_IGNORE);
        if(retval != MPI_SUCCESS) rmg_error_handler (__FILE__, __LINE__, "Error in MPI_Waitall.\n");
    }
    T.thread_barrier_wait();


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
    T.thread_barrier_wait();

} // end trade_images1_central_async


// C interfaces for transitional usage
extern "C" void init_TradeImages(void)
{
    TradeImages *T;
    T = new TradeImages();
}

extern "C"  void trade_imagesx (rmg_double_t * f, rmg_double_t * w, int dimx, int dimy, int dimz, int images, int type)
{
    TradeImages T;
    T.trade_imagesx<double>(f, w, dimx, dimy, dimz, images, type);
}
extern "C"  void trade_imagesx_f (rmg_float_t * f, rmg_float_t * w, int dimx, int dimy, int dimz, int images, int type)
{
    TradeImages T;
    T.trade_imagesx<float>(f, w, dimx, dimy, dimz, images, type);
}
extern "C" void trade_images (rmg_double_t *mat, int dimx, int dimy, int dimz, int type)
{
    TradeImages T;
    T.trade_images<double>(mat, dimx, dimy, dimz, type);
}

extern "C" void trade_images_f (rmg_float_t *mat, int dimx, int dimy, int dimz, int type)
{
    TradeImages T;
    T.trade_images<float>(mat, dimx, dimy, dimz, type);
}

