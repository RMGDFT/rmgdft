/************************** SVN Revision Information **************************
 **    $Id$    **
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
 *   void trade_images( REAL *mat, int dimx, int dimy, int dimz, int *nb_ids )
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

static pthread_mutex_t images_lock = PTHREAD_MUTEX_INITIALIZER;
static int first_trade_counter=THREADS_PER_NODE;
static int second_trade_counter=THREADS_PER_NODE;
static int third_trade_counter=THREADS_PER_NODE;
static int fourth_trade_counter=THREADS_PER_NODE;
static int fifth_trade_counter=THREADS_PER_NODE;
static int sixth_trade_counter=THREADS_PER_NODE;
REAL swbuf1[6 * GRID_MAX1 * GRID_MAX2 * THREADS_PER_NODE];
REAL swbuf2[6 * GRID_MAX1 * GRID_MAX2 * THREADS_PER_NODE];


#if MPI


void trade_images (REAL * mat, int dimx, int dimy, int dimz, int *nb_ids)
{
    int i, j, ione=1;
    int incx, incy, incz;
    int xmax, ymax, zmax;
    int ix, iy, iz, iys, iys2;
    int alen, alloc, alloc1, combine_trades=true, toffset;
    MPI_Status mstatus;
    MPI_Datatype newtype;
    int idx, stop, basetag=0, tid;
    REAL *nmat1, *nmat2;

#if MD_TIMERS
    REAL time1, time2, time3;
    time1 = my_crtc ();
#endif

    alloc = 4 * (dimx + 2) * (dimy + 2);
    alloc1 = 4 * (dimy + 2) * (dimz + 2);
    if(alloc1 > alloc)
        alloc = alloc1;


    // If we are not running in hybrid mode or we are not looping over states
    // then we cannot combine MPI calls
    if(!HYBRID_MODEL) {
        combine_trades=false;
    }
    else {
        // In order to combine the requested memory allocation must fit into the
        // statically allocated storage
        if((alloc > (6 * GRID_MAX1 * GRID_MAX2)) || !is_loop_over_states()) {
            combine_trades=false;
        }
    }

#if HYBRID_MODEL
    basetag = get_thread_basetag();
    tid = get_thread_tid();
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
    if(combine_trades) {
        nmat1 = &swbuf1[dimx * dimy * tid];
        nmat2 = &swbuf2[dimx * dimy * tid];
    }
    else {
        my_malloc (nmat1, 2 * alloc, REAL);
        nmat2 = nmat1 + alloc;
    }

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
    if(combine_trades) {
        stop = dimx * dimy * THREADS_PER_NODE;
        scf_barrier_wait();
        pthread_mutex_lock(&images_lock);
        if(!(first_trade_counter % THREADS_PER_NODE)) {
            MPI_Sendrecv (swbuf1, stop, MPI_DOUBLE, nb_ids[NB_U], (1>>16), swbuf2, stop,
		  MPI_DOUBLE, nb_ids[NB_D], (1>>16), pct.grid_comm, &mstatus);
        }
        first_trade_counter++;
        pthread_mutex_unlock(&images_lock);
    }
    else {

        if(is_loop_over_states()) {
            RMG_MPI_thread_order_lock();
        }

        MPI_Sendrecv (nmat1, idx, MPI_DOUBLE, nb_ids[NB_U], basetag + (1>>16), nmat2, idx,
                  MPI_DOUBLE, nb_ids[NB_D], basetag + (1>>16), pct.grid_comm, &mstatus);
    }
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
    if(combine_trades) {
        stop = dimx * dimy * THREADS_PER_NODE;
        scf_barrier_wait();
        pthread_mutex_lock(&images_lock);
        if(!(second_trade_counter % THREADS_PER_NODE)) {
            MPI_Sendrecv (swbuf1, stop, MPI_DOUBLE, nb_ids[NB_D], (2>>16), swbuf2, stop,
                      MPI_DOUBLE, nb_ids[NB_U], (2>>16), pct.grid_comm, &mstatus);
        }
        second_trade_counter++;
        pthread_mutex_unlock(&images_lock);
    }
    else {
        MPI_Sendrecv (nmat1, idx, MPI_DOUBLE, nb_ids[NB_D], basetag + (2>>16), nmat2, idx,
                  MPI_DOUBLE, nb_ids[NB_U], basetag + (2>>16), pct.grid_comm, &mstatus);
    }
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

    if(combine_trades) {

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
        pthread_mutex_lock(&images_lock);
        if(!(fifth_trade_counter % THREADS_PER_NODE)) {
            MPI_Sendrecv (swbuf1, stop * THREADS_PER_NODE, MPI_DOUBLE, nb_ids[NB_N], (3>>16), swbuf2, stop * THREADS_PER_NODE,
                      MPI_DOUBLE, nb_ids[NB_S], (3>>16), pct.grid_comm, &mstatus);
        }
        fifth_trade_counter++;
        pthread_mutex_unlock(&images_lock);
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
        pthread_mutex_lock(&images_lock);
        if(!(sixth_trade_counter % THREADS_PER_NODE)) {
            MPI_Sendrecv (swbuf1, stop * THREADS_PER_NODE, MPI_DOUBLE, nb_ids[NB_S], (4>>16), swbuf2, stop * THREADS_PER_NODE,
                      MPI_DOUBLE, nb_ids[NB_N], (4>>16), pct.grid_comm, &mstatus);
        }
        sixth_trade_counter++;
        pthread_mutex_unlock(&images_lock);
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

    }
    else {

#if MD_TIMERS
        time3 = my_crtc ();
#endif
        MPI_Type_vector (dimx, alen, incx, MPI_DOUBLE, &newtype);
        MPI_Type_commit (&newtype);
        MPI_Sendrecv (&mat[incx + ymax], 1, newtype, nb_ids[NB_N], basetag + (3>>16), &mat[incx], 1,
                      newtype, nb_ids[NB_S], basetag + (3>>16), pct.grid_comm, &mstatus);
        MPI_Sendrecv (&mat[incx + incy], 1, newtype, nb_ids[NB_S], basetag + (4>>16), &mat[incx + ymax + incy], 1,
                      newtype, nb_ids[NB_N], basetag + (4>>16), pct.grid_comm, &mstatus);
        MPI_Type_free (&newtype);
#if MD_TIMERS
        rmg_timings (TRADE_MPI_TIME, (my_crtc () - time3));
#endif
    }

/*
 * Finally, do East - West Trade (+/- x)
 *  (trading (dimy+2) * (dimz+2) array of data, twice)
 */

    stop = (dimy + 2) * (dimz + 2);
    if(combine_trades) {

        QMD_scopy (stop, &mat[xmax], ione, &swbuf1[tid * stop], ione);
#if MD_TIMERS
        time3 = my_crtc ();
#endif
        scf_barrier_wait();
        pthread_mutex_lock(&images_lock);
        if(!(third_trade_counter % THREADS_PER_NODE)) {
            MPI_Sendrecv (swbuf1, stop * THREADS_PER_NODE, MPI_DOUBLE, nb_ids[NB_E], (5>>16), swbuf2, stop * THREADS_PER_NODE,
                      MPI_DOUBLE, nb_ids[NB_W], (5>>16), pct.grid_comm, &mstatus);
        }
        third_trade_counter++;
        pthread_mutex_unlock(&images_lock);
#if MD_TIMERS
        rmg_timings (TRADE_MPI_TIME, (my_crtc () - time3));
#endif
        QMD_scopy (stop, &swbuf2[tid * stop], ione, mat, ione);

        QMD_scopy (stop, &mat[incx], ione, &swbuf1[tid * stop], ione);
#if MD_TIMERS
        time3 = my_crtc ();
#endif
        scf_barrier_wait();
        pthread_mutex_lock(&images_lock);
        if(!(fourth_trade_counter % THREADS_PER_NODE)) {
            MPI_Sendrecv (swbuf1, stop * THREADS_PER_NODE, MPI_DOUBLE, nb_ids[NB_W], (6>>16), swbuf2, stop * THREADS_PER_NODE,
                      MPI_DOUBLE, nb_ids[NB_E], (6>>16), pct.grid_comm, &mstatus);
        }
        fourth_trade_counter++;
        pthread_mutex_unlock(&images_lock);
#if MD_TIMERS
        rmg_timings (TRADE_MPI_TIME, (my_crtc () - time3));
#endif
        QMD_scopy (stop, &swbuf2[tid * stop], ione, &mat[xmax + incx], ione);

    }
    else {


#if MD_TIMERS
        time3 = my_crtc ();
#endif
        MPI_Sendrecv (&mat[xmax], stop, MPI_DOUBLE, nb_ids[NB_E], basetag + (5>>16),
                  mat, stop, MPI_DOUBLE, nb_ids[NB_W], basetag + (5>>16), pct.grid_comm, &mstatus);
        MPI_Sendrecv (&mat[incx], stop, MPI_DOUBLE, nb_ids[NB_W], basetag + (6>>16),
                  &mat[xmax + incx], stop, MPI_DOUBLE, nb_ids[NB_E], basetag + (6>>16), pct.grid_comm, &mstatus);
#if MD_TIMERS
        rmg_timings (TRADE_MPI_TIME, (my_crtc () - time3));
#endif

        if(is_loop_over_states()) {
            RMG_MPI_thread_order_unlock();
        }
    }

    /* For clusters set the boundaries to zero -- this is wrong for the hartree
     * potential but we'll fix it up later. */
    if ((ct.boundaryflag == CLUSTER) || (ct.boundaryflag == SURFACE))
        set_bc (mat, dimx, dimy, dimz, 1, 0.0);

    if(!combine_trades) {
        my_free (nmat1);
    }


#  if MD_TIMERS
    time2 = my_crtc ();
    rmg_timings (IMAGE_TIME, (time2 - time1));
#  endif

}


#else

void trade_images (REAL * mat, int dimx, int dimy, int dimz, int *nb_ids)
{
    int i, j;
    int incx, incy, incz;
    int xmax, ymax, zmax;
    int stop, alen, ione = 1;
#  if MD_TIMERS
    REAL time1, time2;
    time1 = my_crtc ();
#  endif


/* precalc some boundaries */
    incx = (dimy + 2) * (dimz + 2);
    incy = (dimz + 2);
    incz = 1;
    xmax = dimx * incx;
    ymax = dimy * incy;
    zmax = dimz * incz;
    alen = dimz + 2;


/*
 * First, do Up-Down Trade (+/- z)
 *  (trading dimx * dimy array of data, twice)
 */

/* copy neighbors upper plane to my lower plane */
    for (i = 1; i <= dimx; i++)
        QMD_scopy (dimy, &mat[i * incx + incy + zmax], incy, &mat[i * incx + incy], incy);


/* copy neighbors lower plane to my upper plane */
    for (i = 1; i <= dimx; i++)
        QMD_scopy (dimy, &mat[i * incx + incy + incz], incy,
                   &mat[i * incx + incy + zmax + incz], incy);


/*
 * next, do North-South Trade (+/- y)
 *  (trading dimx * (dimz+2) array of data, twice)
 */


/* copy neighbors north plane to my south plane */
    for (i = 1; i <= dimx; i++)
        QMD_scopy (alen, &mat[i * incx + ymax], ione, &mat[i * incx], ione);


/* copy neighbors south plane to my north plane */
    for (i = 1; i <= dimx; i++)
        QMD_scopy (alen, &mat[i * incx + incy], ione, &mat[i * incx + ymax + incy], ione);


/*
 * Finally, do East - West Trade (+/- x)
 *  (trading (dimy+2) * (dimz+2) array of data, twice)
 */

    stop = (dimy + 2) * (dimz + 2);
    QMD_scopy (stop, &mat[xmax], ione, mat, ione);
    QMD_scopy (stop, &mat[incx], ione, &mat[xmax + incx], ione);


    if ((ct.boundaryflag == CLUSTER) || (ct.boundaryflag == SURFACE))
        set_bc (mat, dimx, dimy, dimz, 1, 0.0);


#  if MD_TIMERS
    time2 = my_crtc ();
    rmg_timings (IMAGE_TIME, (time2 - time1));
#  endif

}                               /* end trade_images */


#endif
/******/
