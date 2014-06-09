/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
/****f* QMD-MGDFT/trade_images2.c *****
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
 *   void trade_images2(rmg_double_t *f, rmg_double_t *w, int dimx, int dimy, int dimz)
 *   trades boundary information with neighboring PEs (to second neighbor)   
 * INPUTS
 *   f[dimx*dimy*dimz]: its image data are not filled
 * OUTPUT
 *   w[(dimx+4)*(dimy+4)*(dimz+4)]: image data are filled after call
 * PARENTS
 *   app6_del2.c app_grad.c
 * CHILDREN
 *   nothing
 * SOURCE
 */


#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "main.h"
#include "prototypes_on.h"



void trade_images2(rmg_double_t * f, rmg_double_t * w, int dimx, int dimy, int dimz)
{
    int ix, iy, iz;
    int ixs, iys;
    int xlen, ylen, zlen;
    int *nb_ids;
    MPI_Status mrstatus;
    rmg_double_t *frdx1, *frdx2, *frdy1, *frdy2, *frdz1, *frdz2;
    rmg_double_t *frdx1n, *frdx2n, *frdy1n, *frdy2n, *frdz1n, *frdz2n;

    int incx = dimy * dimz;
    int incy = dimz;
    int incx2 = (dimy + 4) * (dimz + 4);
    int incy2 = dimz + 4;


    zlen = dimx * dimy * 2;
    ylen = dimx * 2 * (dimz + 4);
    xlen = 2 * (dimy + 4) * (dimz + 4);


    my_malloc_init( frdx1, 2 * (xlen + ylen + zlen), rmg_double_t );
    frdx2 = frdx1 + xlen;
    frdy1 = frdx2 + xlen;
    frdy2 = frdy1 + ylen;
    frdz1 = frdy2 + ylen;
    frdz2 = frdz1 + zlen;

    my_malloc_init( frdx1n, 2 * (xlen + ylen + zlen), rmg_double_t );
    frdx2n = frdx1n + xlen;
    frdy1n = frdx2n + xlen;
    frdy2n = frdy1n + ylen;
    frdz1n = frdy2n + ylen;
    frdz2n = frdz1n + zlen;


    nb_ids = &pct.neighbors[0];


    /* Load up w with the basic stuff */
    for (ix = 0; ix < dimx; ix++)
    {

        for (iy = 0; iy < dimy; iy++)
        {

            for (iz = 0; iz < dimz; iz++)
            {

                w[(ix + 2) * incx2 + (iy + 2) * incy2 + (iz + 2)] = f[ix * incx + iy * incy + iz];

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */


    /* Collect the positive z-plane and negative z-planes */
    for (ix = 0; ix < dimx; ix++)
    {
        ixs = ix * dimy * 2;
        for (iy = 0; iy < dimy; iy++)
        {
            iys = iy * 2;
            for (iz = 0; iz < 2; iz++)
            {

                frdz1[ixs + iys + iz] = w[(ix + 2) * incx2 + (iy + 2) * incy2 + iz + 2];
                frdz2[ixs + iys + iz] = w[(ix + 2) * incx2 + (iy + 2) * incy2 + dimz + iz];

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */


    MPI_Sendrecv(frdz1, zlen, MPI_DOUBLE, nb_ids[NB_D], 1,
                 frdz2n, zlen, MPI_DOUBLE, nb_ids[NB_U], 1, pct.grid_comm, &mrstatus);

    MPI_Sendrecv(frdz2, zlen, MPI_DOUBLE, nb_ids[NB_U], 2,
                 frdz1n, zlen, MPI_DOUBLE, nb_ids[NB_D], 2, pct.grid_comm, &mrstatus);


    /* Unpack them */
    for (ix = 0; ix < dimx; ix++)
    {

        ixs = ix * dimy * 2;
        for (iy = 0; iy < dimy; iy++)
        {
            iys = iy * 2;
            for (iz = 0; iz < 2; iz++)
            {

                w[(ix + 2) * incx2 + (iy + 2) * incy2 + iz] = frdz1n[ixs + iys + iz];
                w[(ix + 2) * incx2 + (iy + 2) * incy2 + iz + dimz + 2] = frdz2n[ixs + iys + iz];

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */



    /* Collect the north and south planes */
    for (ix = 0; ix < dimx; ix++)
    {
        ixs = ix * 2 * (dimz + 4);
        for (iy = 0; iy < 2; iy++)
        {

            iys = iy * (dimz + 4);
            for (iz = 0; iz < dimz + 4; iz++)
            {

                frdy1[ixs + iys + iz] = w[(ix + 2) * incx2 + (iy + 2) * incy2 + iz];
                frdy2[ixs + iys + iz] = w[(ix + 2) * incx2 + (dimy + iy) * incy2 + iz];

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */


    MPI_Sendrecv(frdy1, ylen, MPI_DOUBLE, nb_ids[NB_S], 1,
                 frdy2n, ylen, MPI_DOUBLE, nb_ids[NB_N], 1, pct.grid_comm, &mrstatus);

    MPI_Sendrecv(frdy2, ylen, MPI_DOUBLE, nb_ids[NB_N], 2,
                 frdy1n, ylen, MPI_DOUBLE, nb_ids[NB_S], 2, pct.grid_comm, &mrstatus);


    /* Unpack them */
    for (ix = 0; ix < dimx; ix++)
    {
        ixs = ix * 2 * (dimz + 4);
        for (iy = 0; iy < 2; iy++)
        {
            iys = iy * (dimz + 4);
            for (iz = 0; iz < dimz + 4; iz++)
            {

                w[(ix + 2) * incx2 + iy * incy2 + iz] = frdy1n[ixs + iys + iz];
                w[(ix + 2) * incx2 + (iy + dimy + 2) * incy2 + iz] = frdy2n[ixs + iys + iz];

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */


    /* Collect the east and west planes */
    for (ix = 0; ix < 2; ix++)
    {
        ixs = ix * (dimy + 4) * (dimz + 4);
        for (iy = 0; iy < dimy + 4; iy++)
        {
            iys = iy * (dimz + 4);
            for (iz = 0; iz < dimz + 4; iz++)
            {

                frdx1[ixs + iys + iz] = w[(ix + 2) * incx2 + iy * incy2 + iz];
                frdx2[ixs + iys + iz] = w[(dimx + ix) * incx2 + iy * incy2 + iz];

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */


    MPI_Sendrecv(frdx1, xlen, MPI_DOUBLE, nb_ids[NB_W], 1,
                 frdx2n, xlen, MPI_DOUBLE, nb_ids[NB_E], 1, pct.grid_comm, &mrstatus);

    MPI_Sendrecv(frdx2, xlen, MPI_DOUBLE, nb_ids[NB_E], 2,
                 frdx1n, xlen, MPI_DOUBLE, nb_ids[NB_W], 2, pct.grid_comm, &mrstatus);


    /* Unpack them */
    for (ix = 0; ix < 2; ix++)
    {
        ixs = ix * (dimy + 4) * (dimz + 4);
        for (iy = 0; iy < dimy + 4; iy++)
        {
            iys = iy * (dimz + 4);
            for (iz = 0; iz < dimz + 4; iz++)
            {

                w[ix * incx2 + iy * incy2 + iz] = frdx1n[ixs + iys + iz];
                w[(ix + dimx + 2) * incx2 + iy * incy2 + iz] = frdx2n[ixs + iys + iz];

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */


    my_free(frdx1n);
    my_free(frdx1);


}                               /* end trade_images2 */


/******/
