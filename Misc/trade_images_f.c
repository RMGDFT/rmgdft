/************************** SVN Revision Information **************************
 **    $Id: trade_images.c 1870 2013-02-05 13:58:02Z ebriggs $    **
******************************************************************************/

/****f* QMD-MGDFT/trade_images.c *****
 * NAME
 *   Ab initio real space code with multigrid acceleration
 *   Quantum molecular dynamics package.
 *   Version: 2.1.5
 * COPYRIGHT
 *   Copyright (C) 1995  Emil Briggs
 *   Copyright (C) 1998  Emil Briggs, Charles Brabec, Mark Wensell, 
 *                       Dan Sullivan, Chris Rapcewicz, Jerzy Bernholc
 *   Copyright (C) 2001  Emil Briggs, Wenchang Lu,
 *                       Marco Buongiorno Nardelli,Charles Brabec, 
 *                       Mark Wensell,Dan Sullivan, Chris Rapcewicz,
 *                       Jerzy Bernholc
 * FUNCTION
 *   void trade_images( rmg_float_t *mat, int dimx, int dimy, int dimz, int *nb_ids )
 *   trades boundary information with neighboring PEs. Includes corner pieces, 
 *   as needed by restrict and expand in multi-grid
 * INPUTS
 *   mat[(dimx+2)*(dimy+2)*(dimz+2)]: +2 for images of neighboring cells
 *      there is no value on images as input
 *   dimx, dimy, dimz:  size of array except for images
 *   nb_ids:  array of neighboring PE numbers
 * OUTPUT
 *   image data in mat will be filled 
 * PARENTS
 *   app_cil.c app_cilr_bcc.c app_cilr_fcc.c app_cilr_hex.c app_cilr_ortho.c
 *   app_cir.c app_cir_bcc.c app_cir_fcc.c app_cir_hex.c app_cir_ortho.c
 *   mg_eig_state.c mgrid_solv.c subdiag_mpi.c xcgga.c
 * CHILDREN
 *   nothing
 * SOURCE
 */


#include "main.h"
#if HYBRID_MODEL
#include "hybrid.h"
#endif

#include <pthread.h>

static rmg_float_t *swbuf1=NULL;
static rmg_float_t *swbuf2=NULL;
static int max_alloc;


void trade_images_f (rmg_float_t * mat, int dimx, int dimy, int dimz, int *nb_ids)
{
    int i, j, ione=1;
    int incx, incy, incz;
    int xmax, ymax, zmax;
    int ix, iy, iz, iys, iys2;
    int alen, alloc, alloc1, toffset;
    int grid_max1, grid_max2;
    int ACTIVE_THREADS = 1;

    MPI_Status mstatus;
    MPI_Datatype newtype;
    int idx, stop, basetag=0, tid=0;
    rmg_float_t *nmat1, *nmat2;

    // To initialize call from init.c with NULL args
    if(swbuf1 == NULL) {
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
         my_malloc(swbuf1, 6 * grid_max1 * grid_max2 * ct.THREADS_PER_NODE, rmg_float_t);
         my_malloc(swbuf2, 6 * grid_max1 * grid_max2 * ct.THREADS_PER_NODE, rmg_float_t);
         max_alloc = 6 * grid_max1 * grid_max2 * ct.THREADS_PER_NODE;
         return;
    }

#if ASYNC_TRADES
//    trade_images1_async (mat, dimx, dimy, dimz);
//    return;
#endif

#if MD_TIMERS
    rmg_double_t time1, time2, time3;
    time1 = my_crtc ();
#endif

    alloc = 4 * (dimx + 2) * (dimy + 2);
    alloc1 = 4 * (dimy + 2) * (dimz + 2);
    if(alloc1 > alloc)
        alloc = alloc1;
    alloc = alloc * ct.THREADS_PER_NODE;
    if(alloc > max_alloc)
        error_handler("Not enough memory. This should never happen.");


#if HYBRID_MODEL
    basetag = get_thread_basetag();
    tid = get_thread_tid();
    if(tid < 0) tid = 0;
    if(is_loop_over_states()) ACTIVE_THREADS = ct.THREADS_PER_NODE;
#endif


/* precalc some boundaries */
    incx = (dimy + 2) * (dimz + 2);
    incy = (dimz + 2);
    incz = 1;
    xmax = dimx * incx;
    ymax = dimy * incy;
    zmax = dimz * incz;
    alen = dimz + 2;

// For they hybrid case we use a single array for all threads which lets us combine multiple
// Mpi requests into one
    nmat1 = &swbuf1[dimx * dimy * tid];
    nmat2 = &swbuf2[dimx * dimy * tid];

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
    scf_barrier_wait();
    if(tid == 0) {
        MPI_Sendrecv (swbuf1, stop, MPI_FLOAT, nb_ids[NB_U], (1>>16), swbuf2, stop,
		  MPI_FLOAT, nb_ids[NB_D], (1>>16), pct.grid_comm, &mstatus);
    } 
    scf_barrier_wait();

#if MD_TIMERS
    rmg_timings (TRADE_MPI_TIME, (my_crtc () - time3));
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
    scf_barrier_wait();
    if(tid == 0) {
        MPI_Sendrecv (swbuf1, stop, MPI_FLOAT, nb_ids[NB_D], (2>>16), swbuf2, stop,
                      MPI_FLOAT, nb_ids[NB_U], (2>>16), pct.grid_comm, &mstatus);
    }
    scf_barrier_wait();
#if MD_TIMERS
    rmg_timings (TRADE_MPI_TIME, (my_crtc () - time3));
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
            swbuf1[idx + toffset] = mat[iys2 + ymax + iz];
            idx++;
        }
    }
#if MD_TIMERS
    time3 = my_crtc ();
#endif
    scf_barrier_wait();
    if(tid == 0) {
        MPI_Sendrecv (swbuf1, stop * ACTIVE_THREADS, MPI_FLOAT, nb_ids[NB_N], (3>>16), swbuf2, stop * ACTIVE_THREADS,
                  MPI_FLOAT, nb_ids[NB_S], (3>>16), pct.grid_comm, &mstatus);
    }
    scf_barrier_wait();
#if MD_TIMERS
    rmg_timings (TRADE_MPI_TIME, (my_crtc () - time3));
#endif

    idx = 0;
    for (ix = 0; ix < dimx; ix++) {
        iys2 = (ix + 1) * incx;
        for (iz = 0; iz < dimz + 2; iz++) {
            mat[iys2 + iz] = swbuf2[idx + toffset];
            idx++;
        }
    }

    idx = 0;
    for (ix = 0; ix < dimx; ix++) {
        iys2 = (ix + 1) * incx;
        for (iz = 0; iz < dimz + 2; iz++) {
            swbuf1[idx + toffset] = mat[iys2 + incy + iz];
            idx++;
        }
    }
#if MD_TIMERS
    time3 = my_crtc ();
#endif
    scf_barrier_wait();
    if(tid == 0) {
        MPI_Sendrecv (swbuf1, stop * ACTIVE_THREADS, MPI_FLOAT, nb_ids[NB_S], (4>>16), swbuf2, stop * ACTIVE_THREADS,
                  MPI_FLOAT, nb_ids[NB_N], (4>>16), pct.grid_comm, &mstatus);
    }
    scf_barrier_wait();
#if MD_TIMERS
    rmg_timings (TRADE_MPI_TIME, (my_crtc () - time3));
#endif

    idx = 0;
    for (ix = 0; ix < dimx; ix++) {
        iys2 = (ix + 1) * incx;
        for (iz = 0; iz < dimz + 2; iz++) {
            mat[iys2 + ymax + incy + iz] = swbuf2[idx + toffset];
            idx++;
        }
    }


/*
 * Finally, do East - West Trade (+/- x)
 *  (trading (dimy+2) * (dimz+2) array of data, twice)
 */

    stop = (dimy + 2) * (dimz + 2);

    QMD_scopy (stop, &mat[xmax], ione, &swbuf1[tid * stop], ione);
#if MD_TIMERS
    time3 = my_crtc ();
#endif
    scf_barrier_wait();
    if(tid == 0) {
        MPI_Sendrecv (swbuf1, stop * ACTIVE_THREADS, MPI_FLOAT, nb_ids[NB_E], (5>>16), swbuf2, stop * ACTIVE_THREADS,
                  MPI_FLOAT, nb_ids[NB_W], (5>>16), pct.grid_comm, &mstatus);
    }
    scf_barrier_wait();
#if MD_TIMERS
    rmg_timings (TRADE_MPI_TIME, (my_crtc () - time3));
#endif
    QMD_scopy (stop, &swbuf2[tid * stop], ione, mat, ione);

    QMD_scopy (stop, &mat[incx], ione, &swbuf1[tid * stop], ione);
#if MD_TIMERS
    time3 = my_crtc ();
#endif
    scf_barrier_wait();
    if(tid == 0) {
        MPI_Sendrecv (swbuf1, stop * ACTIVE_THREADS, MPI_FLOAT, nb_ids[NB_W], (6>>16), swbuf2, stop * ACTIVE_THREADS,
                  MPI_FLOAT, nb_ids[NB_E], (6>>16), pct.grid_comm, &mstatus);
    }
    scf_barrier_wait();
#if MD_TIMERS
    rmg_timings (TRADE_MPI_TIME, (my_crtc () - time3));
#endif
    QMD_scopy (stop, &swbuf2[tid * stop], ione, &mat[xmax + incx], ione);


    /* For clusters set the boundaries to zero -- this is wrong for the hartree
     * potential but we'll fix it up later. */
//    if ((ct.boundaryflag == CLUSTER) || (ct.boundaryflag == SURFACE))
//        set_bc (mat, dimx, dimy, dimz, 1, 0.0);


#  if MD_TIMERS
    time2 = my_crtc ();
    rmg_timings (IMAGE_TIME, (time2 - time1));
#  endif

}


