/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

/****f* QMD-MGDFT/set_bc.c *****
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
 *   void set_bc( REAL *mat, int dimx, int dimy, int dimz, REAL val, int *nb_ids )
 *   Sets the boundaries of a smoothing grid to a specified value
 *   only usefule for SURFACE or CLUSTER conditions
 * INPUTS
 *   mat[dimx*dimy*dimz]: 
 *   dimx, dimy, dimz: size of array (excluding images) 
 *   images : number of images on each side (i.e. for trade_image2 images should be 2) 
 *   val: boundary data
 * OUTPUT
 *   mat: its boundary data are set to val 
 * PARENTS
 *   get_vh.c trade_images.c trade_image_mpi_sp.c getpoi_bc.c
 * CHILDREN
 *   pe2xyz.c
 * SOURCE
 */



#include "main.h"
#include <stdlib.h>

/* This sets values at boundary points either to b_array array, or to single value (usually 0) if
 * b_array is a NULL pointer
 * This version of set_bc works for multiple boundary points, i.e. when doing trade_images2 or higher*/


void set_bcx (REAL * mat, int dimx, int dimy, int dimz, int images, REAL val)
{
    int ix, iy, iz, i, tim;
    int dimx_max, dimy_max, dimz_max;
    int pex, pey, pez;
    int incx, incy;
    int xstop, ystop, zstop;
    int ibnd;


    /* Figure out which processor we are */
    /*pe2xyz(pct.thispe, &pex, &pey, &pez); */
    pex = pct.pe_x;
    pey = pct.pe_y;
    pez = pct.pe_z;


    ibnd = FALSE;
    if ((pex == 0) || (pey == 0) || (pez == 0))
        ibnd = TRUE;
    if ((pex == (PE_X - 1)) || (pey == (PE_Y - 1)) || (pez == (PE_Z - 1)))
        ibnd = TRUE;

    if (!ibnd)
        return;

    tim = 2 * images;

    /* precalc some boundaries */
    dimx_max = dimx + tim;
    dimy_max = dimy + tim;
    dimz_max = dimz + tim;

    incx = dimy_max * dimz_max;
    incy = dimz_max;

    xstop = dimx_max - images;
    ystop = dimy_max - images;
    zstop = dimz_max - images;




    /* If the boundary condition flag is SURFACE then we only set the z-boundary */
    if (ct.boundaryflag != SURFACE)
    {

        if (pex == 0)
        {

            for (iy = 0; iy < dimy_max; iy++)
                for (iz = 0; iz < dimz_max; iz++)

                    for (i = 0; i < images; i++)
                        mat[i * incx + iy * incy + iz] = val;


        }                       /* end if */


        if (pex == (PE_X - 1))
        {

            for (iy = 0; iy < dimy_max; iy++)
                for (iz = 0; iz < dimz_max; iz++)

                    for (i = xstop; i < dimx_max; i++)
                        mat[i * incx + iy * incy + iz] = val;


        }                       /* end if */


        if (pey == 0)
        {

            for (ix = 0; ix < dimx_max; ix++)
                for (iz = 0; iz < dimz_max; iz++)

                    for (i = 0; i < images; i++)
                        mat[ix * incx + i * incy + iz] = val;


        }                       /* end if */


        if (pey == (PE_Y - 1))
        {

            for (ix = 0; ix < dimx_max; ix++)
                for (iz = 0; iz < dimz_max; iz++)

                    for (i = ystop; i < dimy_max; i++)
                        mat[ix * incx + i * incy + iz] = val;


        }                       /* end if */

    }                           /* end if */


    if (pez == 0)
    {

        for (ix = 0; ix < dimx_max; ix++)
            for (iy = 0; iy < dimy_max; iy++)

                for (i = 0; i < images; i++)
                    mat[ix * incx + iy * incy + i] = val;


    }                           /* end if */


    if (pez == (PE_Z - 1))
    {

        for (ix = 0; ix < dimx_max; ix++)
            for (iy = 0; iy < dimy_max; iy++)

                for (i = zstop; i < dimz_max; i++)
                    mat[ix * incx + iy * incy + i] = val;


    }                           /* end if */







}                               /* end set_bc */

/******/
