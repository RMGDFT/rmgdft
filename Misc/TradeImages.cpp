/************************** SVN Revision Information **************************
 **    $Id: trade_imagesx.c 2091 2014-02-05 15:48:31Z ebriggs $    **
******************************************************************************/

#include "TradeImages.h"
#include "transition.h"
#include "common_prototypes.h"
#include "BlasWrappers.h"

using namespace std;

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
         retval = MPI_Alloc_mem(6 * sizeof(rmg_double_t) * MAX_TRADE_IMAGES * T.get_threads_per_node() * grid_max1*grid_max2 , MPI_INFO_NULL, &swbuf1x);
         if(retval != MPI_SUCCESS) {
             rmg_error_handler("Error in MPI_Alloc_mem.\n");
         }
         retval = MPI_Alloc_mem(6 * sizeof(rmg_double_t) * MAX_TRADE_IMAGES * T.get_threads_per_node() * grid_max1*grid_max2 , MPI_INFO_NULL, &swbuf2x);
         if(retval != MPI_SUCCESS) {
             rmg_error_handler("Error in MPI_Alloc_mem.\n");
         }


         //my_malloc(swbuf1x, 6 * MAX_TRADE_IMAGES * grid_max1 * grid_max2 * ct.THREADS_PER_NODE, rmg_double_t);
         //my_malloc(swbuf2x, 6 * MAX_TRADE_IMAGES * grid_max1 * grid_max2 * ct.THREADS_PER_NODE, rmg_double_t);
         max_alloc = 6 * MAX_TRADE_IMAGES * grid_max1 * grid_max2 * T.get_threads_per_node();
         return;
    }

}

template <typename RmgType>
void TradeImages::CPP_trade_imagesx (RmgType *f, RmgType *w, int dimx, int dimy, int dimz, int images, int type)
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

#if ASYNC_TRADES
//    if(type == CENTRAL_FD) {
//        trade_imagesx_central_async_f (f, w, dimx, dimy, dimz, images);
//        return;
//    }
//    else {
//        trade_imagesx_async_f (f, w, dimx, dimy, dimz, images);
//        return;
//    }
#endif

#if MD_TIMERS
    rmg_double_t time1, time2, time3;
    time1 = my_crtc ();
#endif
    tid = 0;
#if HYBRID_MODEL
//    basetag = T.get_thread_basetag();
    tid = T.get_thread_tid();
    if(tid < 0) tid = 0;
    if(T.is_loop_over_states()) ACTIVE_THREADS = T.get_threads_per_node();
#endif

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
        rmg_error_handler("Not enough memory. This should never happen.");


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


#if MD_TIMERS
    time3 = my_crtc ();
#endif

    stop = zlen * ACTIVE_THREADS;
    T.scf_barrier_wait();
    if(tid == 0) {
        MPI_Sendrecv (&swbuf1x_f[0], factor*stop, MPI_BYTE, nb_ids[NB_D], (1>>16),
                  &swbuf2x_f[0], factor*stop, MPI_BYTE, nb_ids[NB_U], (1>>16), transition_get_grid_comm(), &mrstatus);

        MPI_Sendrecv (&swbuf1x_f[stop], factor*stop, MPI_BYTE, nb_ids[NB_U], (2>>16),
                  &swbuf2x_f[stop], factor*stop, MPI_BYTE, nb_ids[NB_D], (2>>16), transition_get_grid_comm(), &mrstatus);
    }
    T.scf_barrier_wait();

#if MD_TIMERS
    T.rmg_timings (TRADE_MPI_TIME, (my_crtc () - time3));
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


#if MD_TIMERS
    time3 = my_crtc ();
#endif
    stop = ylen * ACTIVE_THREADS;
    T.scf_barrier_wait();
    if(tid == 0) {
        MPI_Sendrecv (&swbuf1x_f[0], factor*stop, MPI_BYTE, nb_ids[NB_S], (3>>16),
                  &swbuf2x_f[0], factor*stop, MPI_BYTE, nb_ids[NB_N], (3>>16), transition_get_grid_comm(), &mrstatus);

        MPI_Sendrecv (&swbuf1x_f[stop], factor*stop, MPI_BYTE, nb_ids[NB_N], (4>>16),
                  &swbuf2x_f[stop], factor*stop, MPI_BYTE, nb_ids[NB_S], (4>>16), transition_get_grid_comm(), &mrstatus);
    }
    T.scf_barrier_wait();

#if MD_TIMERS
    T.rmg_timings (TRADE_MPI_TIME, (my_crtc () - time3));
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


#if MD_TIMERS
    time3 = my_crtc ();
#endif
    stop = xlen * ACTIVE_THREADS;
    T.scf_barrier_wait();
    if(tid == 0) {

        MPI_Sendrecv (&swbuf1x_f[0], factor*stop, MPI_BYTE, nb_ids[NB_W], (5>>16),
                   &swbuf2x_f[0], factor*stop, MPI_BYTE, nb_ids[NB_E], (5>>16), transition_get_grid_comm(), &mrstatus);

        MPI_Sendrecv (&swbuf1x_f[stop], factor*stop, MPI_BYTE, nb_ids[NB_E], (6>>16),
                  &swbuf2x_f[stop], factor*stop, MPI_BYTE, nb_ids[NB_W], (6>>16), transition_get_grid_comm(), &mrstatus);

    }
    T.scf_barrier_wait();

#if MD_TIMERS
    T.rmg_timings (TRADE_MPI_TIME, (my_crtc () - time3));
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
    T.rmg_timings (IMAGE_TIME, (time2 - time1));
#endif

} // end CPP_trade_imagesx



template <typename RmgType>
void TradeImages::CPP_trade_images (RmgType * mat, int dimx, int dimy, int dimz, int *nb_ids, int type)
{
    int i, j, ione=1;
    int incx, incy, incz;
    int xmax, ymax, zmax;
    int ix, iz, iys2;
    int alloc, alloc1, toffset;
    int ACTIVE_THREADS = 1, factor;
    BaseThread T(0);

    factor = sizeof(RmgType);

    MPI_Status mstatus;
    int idx, stop, tid=0;
    RmgType *nmat1, *nmat2;
    RmgType *swbuf1x_f, *swbuf2x_f;


    swbuf1x_f = (RmgType *)TradeImages::swbuf1x;
    swbuf2x_f = (RmgType *)TradeImages::swbuf2x;


#if ASYNC_TRADES
//    if(type == CENTRAL_FD) {
//        trade_images1_central_async_f (mat, dimx, dimy, dimz);
//        return;
//    }
#endif

#if MD_TIMERS
    rmg_double_t time1, time2, time3;
    time1 = my_crtc ();
#endif

    alloc = 4 * (dimx + 2) * (dimy + 2);
    alloc1 = 4 * (dimy + 2) * (dimz + 2);
    if(alloc1 > alloc)
        alloc = alloc1;
    alloc = alloc * T.get_threads_per_node();
    if(alloc > TradeImages::max_alloc)
        rmg_error_handler("Not enough memory. This should never happen.");


#if HYBRID_MODEL
    tid = T.get_thread_tid();
    if(tid < 0) tid = 0;
    if(T.is_loop_over_states()) ACTIVE_THREADS = T.get_threads_per_node();
#endif


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

#if MD_TIMERS
    time3 = my_crtc ();
#endif
    // We only want one thread to do the MPI call here.
    stop = dimx * dimy * ACTIVE_THREADS;
    T.scf_barrier_wait();
    if(tid == 0) {
        MPI_Sendrecv (swbuf1x_f, factor*stop, MPI_BYTE, nb_ids[NB_U], (1>>16), swbuf2x_f, factor*stop,
		  MPI_BYTE, nb_ids[NB_D], (1>>16), transition_get_grid_comm(), &mstatus);
    } 
    T.scf_barrier_wait();

#if MD_TIMERS
    T.rmg_timings (TRADE_MPI_TIME, (my_crtc () - time3));
#endif

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


#if MD_TIMERS
    time3 = my_crtc ();
#endif
    stop = dimx * dimy * ACTIVE_THREADS;
    T.scf_barrier_wait();
    if(tid == 0) {
        MPI_Sendrecv (swbuf1x_f, factor*stop, MPI_BYTE, nb_ids[NB_D], (2>>16), swbuf2x_f, factor*stop,
                      MPI_BYTE, nb_ids[NB_U], (2>>16), transition_get_grid_comm(), &mstatus);
    }
    T.scf_barrier_wait();
#if MD_TIMERS
    T.rmg_timings (TRADE_MPI_TIME, (my_crtc () - time3));
#endif


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
#if MD_TIMERS
    time3 = my_crtc ();
#endif
    T.scf_barrier_wait();
    if(tid == 0) {
        MPI_Sendrecv (swbuf1x_f, factor*stop * ACTIVE_THREADS, MPI_BYTE, nb_ids[NB_N], (3>>16), swbuf2x_f, factor*stop * ACTIVE_THREADS,
                  MPI_BYTE, nb_ids[NB_S], (3>>16), transition_get_grid_comm(), &mstatus);
    }
    T.scf_barrier_wait();
#if MD_TIMERS
    T.rmg_timings (TRADE_MPI_TIME, (my_crtc () - time3));
#endif

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
#if MD_TIMERS
    time3 = my_crtc ();
#endif
    T.scf_barrier_wait();
    if(tid == 0) {
        MPI_Sendrecv (swbuf1x_f, factor*stop * ACTIVE_THREADS, MPI_BYTE, nb_ids[NB_S], (4>>16), swbuf2x_f, factor*stop * ACTIVE_THREADS,
                  MPI_BYTE, nb_ids[NB_N], (4>>16), transition_get_grid_comm(), &mstatus);
    }
    T.scf_barrier_wait();
#if MD_TIMERS
    T.rmg_timings (TRADE_MPI_TIME, (my_crtc () - time3));
#endif

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
#if MD_TIMERS
    time3 = my_crtc ();
#endif
    T.scf_barrier_wait();
    if(tid == 0) {
        MPI_Sendrecv (swbuf1x_f, factor*stop * ACTIVE_THREADS, MPI_BYTE, nb_ids[NB_E], (5>>16), swbuf2x_f, factor*stop * ACTIVE_THREADS,
                  MPI_BYTE, nb_ids[NB_W], (5>>16), transition_get_grid_comm(), &mstatus);
    }
    T.scf_barrier_wait();
#if MD_TIMERS
    T.rmg_timings (TRADE_MPI_TIME, (my_crtc () - time3));
#endif
    QMD_copy (stop, &swbuf2x_f[tid * stop], ione, mat, ione);

    QMD_copy (stop, &mat[incx], ione, &swbuf1x_f[tid * stop], ione);
#if MD_TIMERS
    time3 = my_crtc ();
#endif
    T.scf_barrier_wait();
    if(tid == 0) {
        MPI_Sendrecv (swbuf1x_f, factor*stop * ACTIVE_THREADS, MPI_BYTE, nb_ids[NB_W], (6>>16), swbuf2x_f, factor*stop * ACTIVE_THREADS,
                  MPI_BYTE, nb_ids[NB_E], (6>>16), transition_get_grid_comm(), &mstatus);
    }
    T.scf_barrier_wait();
#if MD_TIMERS
    T.rmg_timings (TRADE_MPI_TIME, (my_crtc () - time3));
#endif
    QMD_copy (stop, &swbuf2x_f[tid * stop], ione, &mat[xmax + incx], ione);


    /* For clusters set the boundaries to zero -- this is wrong for the hartree
     * potential but we'll fix it up later. */
//    if ((ct.boundaryflag == CLUSTER) || (ct.boundaryflag == SURFACE))
//        set_bc (mat, dimx, dimy, dimz, 1, 0.0);


#  if MD_TIMERS
    time2 = my_crtc ();
    T.rmg_timings (IMAGE_TIME, (time2 - time1));
#  endif

}


// C interfaces for transitional usage
extern "C" void init_TradeImages(void)
{
    TradeImages *T;
    T = new TradeImages();
}

extern "C"  void trade_imagesx (rmg_double_t * f, rmg_double_t * w, int dimx, int dimy, int dimz, int images, int type)
{
    TradeImages T;
    T.CPP_trade_imagesx<double>(f, w, dimx, dimy, dimz, images, type);
}

extern "C"  void trade_imagesx_f (rmg_float_t * f, rmg_float_t * w, int dimx, int dimy, int dimz, int images, int type)
{
    TradeImages T;
    T.CPP_trade_imagesx<float>(f, w, dimx, dimy, dimz, images, type);
}

extern "C" void trade_images (rmg_double_t *mat, int dimx, int dimy, int dimz, int *nb_ids, int type)
{
    TradeImages T;
    T.CPP_trade_images<double>(mat, dimx, dimy, dimz, nb_ids, type);
}

extern "C" void trade_images_f (rmg_float_t *mat, int dimx, int dimy, int dimz, int *nb_ids, int type)
{
    TradeImages T;
    T.CPP_trade_images<float>(mat, dimx, dimy, dimz, nb_ids, type);
}

