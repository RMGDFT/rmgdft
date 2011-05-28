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
#include "hybrid.h"


#if MPI


void trade_images (REAL * mat, int dimx, int dimy, int dimz, int *nb_ids)
{
    int i, j;
    int incx, incy, incz;
    int xmax, ymax, zmax;
    int alen, alloc;
    MPI_Status mstatus;
    MPI_Datatype newtype;
    int idx, stop, basetag=0;
    REAL *nmat1, *nmat2;

#if MD_TIMERS
    REAL time1, time2;
    time1 = my_crtc ();
#endif

#if HYBRID_MODEL
    basetag = get_thread_basetag();
#endif

    alloc = 4 * (dimx + 2) * (dimy + 2);
    my_malloc (nmat1, 2 * alloc, REAL);
    nmat2 = nmat1 + alloc;

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


    idx = 0;
    for (i = incx; i <= incx * dimx; i += incx)
    {
        for (j = 1; j <= dimy; j++)
        {
            nmat1[idx] = mat[i + j * incy + zmax];
            idx++;
        }
    }

    MPI_Sendrecv (nmat1, idx, MPI_DOUBLE, nb_ids[NB_U], basetag + (1>>16), nmat2, idx,
                  MPI_DOUBLE, nb_ids[NB_D], basetag + (1>>16), pct.grid_comm, &mstatus);

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


    MPI_Sendrecv (nmat1, idx, MPI_DOUBLE, nb_ids[NB_D], basetag + (2>>16), nmat2, idx,
                  MPI_DOUBLE, nb_ids[NB_U], basetag + (2>>16), pct.grid_comm, &mstatus);

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


    MPI_Type_vector (dimx, alen, incx, MPI_DOUBLE, &newtype);
    MPI_Type_commit (&newtype);
    MPI_Sendrecv (&mat[incx + ymax], 1, newtype, nb_ids[NB_N], basetag + (3>>16), &mat[incx], 1,
                  newtype, nb_ids[NB_S], basetag + (3>>16), pct.grid_comm, &mstatus);
    MPI_Sendrecv (&mat[incx + incy], 1, newtype, nb_ids[NB_S], basetag + (4>>16), &mat[incx + ymax + incy], 1,
                  newtype, nb_ids[NB_N], basetag + (4>>16), pct.grid_comm, &mstatus);

    MPI_Type_free (&newtype);

/*
 * Finally, do East - West Trade (+/- x)
 *  (trading (dimy+2) * (dimz+2) array of data, twice)
 */


    stop = (dimy + 2) * (dimz + 2);
    MPI_Sendrecv (&mat[xmax], stop, MPI_DOUBLE, nb_ids[NB_E], basetag + (5>>16),
                  mat, stop, MPI_DOUBLE, nb_ids[NB_W], basetag + (5>>16), pct.grid_comm, &mstatus);
    MPI_Sendrecv (&mat[incx], stop, MPI_DOUBLE, nb_ids[NB_W], basetag + (6>>16),
                  &mat[xmax + incx], stop, MPI_DOUBLE, nb_ids[NB_E], basetag + (6>>16), pct.grid_comm, &mstatus);


    /* For clusters set the boundaries to zero -- this is wrong for the hartree
     * potential but we'll fix it up later. */
    if ((ct.boundaryflag == CLUSTER) || (ct.boundaryflag == SURFACE))
        set_bc (mat, dimx, dimy, dimz, 1, 0.0);


    my_free (nmat1);

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
