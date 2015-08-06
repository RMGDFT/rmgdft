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
 *   void set_bc( double *mat, int dimx, int dimy, int dimz, double val, int *nb_ids )
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



#include "grid.h"
#include "main.h"
#include "common_prototypes.h"



void set_bc (double * mat, int dimx, int dimy, int dimz, int images, double val)
{
    int ix, iy, iz, tim;
    int dimx_max, dimy_max, dimz_max;
    int pex, pey, pez;
    int incx, incy;
    int xmax, ymax, zmax;
    int ibnd;


    /* Figure out what processor we are */
    pe2xyz (pct.gridpe, &pex, &pey, &pez);


    ibnd = FALSE;
    if ((pex == 0) || (pey == 0) || (pez == 0))
        ibnd = TRUE;
    if ((pex == (get_PE_X() - 1)) || (pey == (get_PE_Y() - 1)) || (pez == (get_PE_Z() - 1)))
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

    xmax = (dimx + tim - 1) * incx;
    ymax = (dimy + tim - 1) * incy;
    zmax = (dimz + tim - 1);



    /* If the boundary condition flag is SURFACE then we only set the z-boundary */
    if (ct.boundaryflag != SURFACE)
    {

        if (pex == 0)
        {

            for (iy = 0; iy < dimy_max; iy++)
                for (iz = 0; iz < dimz_max; iz++)
                    mat[iy * incy + iz] = val;

        }                       /* end if */


        if (pex == (get_PE_X() - 1))
        {

            for (iy = 0; iy < dimy_max; iy++)
                for (iz = 0; iz < dimz_max; iz++)
                    mat[xmax + iy * incy + iz] = val;


        }                       /* end if */


        if (pey == 0)
        {

            for (ix = 0; ix < dimx_max; ix++)
                for (iz = 0; iz < dimz_max; iz++)
                    mat[ix * incx + iz] = val;

        }                       /* end if */


        if (pey == (get_PE_Y() - 1))
        {

            for (ix = 0; ix < dimx_max; ix++)
                for (iz = 0; iz < dimz_max; iz++)
                    mat[ix * incx + ymax + iz] = val;

        }                       /* end if */

    }                           /* end if */


    if (pez == 0)
    {

        for (ix = 0; ix < dimx_max; ix++)
            for (iy = 0; iy < dimy_max; iy++)
                mat[ix * incx + iy * incy] = val;

    }                           /* end if */


    if (pez == (get_PE_Z() - 1))
    {

        for (ix = 0; ix < dimx_max; ix++)
            for (iy = 0; iy < dimy_max; iy++)
                mat[ix * incx + iy * incy + zmax] = val;

    }                           /* end if */

}                               /* end set_bc */


/******/
