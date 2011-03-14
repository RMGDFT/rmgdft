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
 *   mat[(dimx+2)*(dimy+2)*(dimz+2)]: 
 *   dimx, dimy, dimzr: size of array except for images 
 *   nb_ids: array of neighboring PE numbers   
 *   val: boundart data
 * OUTPUT
 *   mat: its boundary data are set to val 
 * PARENTS
 *   get_vh.c trade_images.c trade_image_mpi_sp.c getpoi_bc.c
 * CHILDREN
 *   pe2xyz.c
 * SOURCE
 */



#include "md.h"



void set_bc(REAL * mat, int dimx, int dimy, int dimz, REAL val, int *nb_ids)
{
    int ix, iy, iz;
    int pex, pey, pez;
    int incx, incy, incz;
    int xmax, ymax;
    int ibnd, idx, stop;


    /* Figure out what processor we are */
    pe2xyz(pct.gridpe, &pex, &pey, &pez);


    ibnd = FALSE;
    if ((pex == 0) || (pey == 0) || (pez == 0))
        ibnd = TRUE;
    if ((pex == (pct.pe_x - 1)) || (pey == (pct.pe_y - 1)) || (pez == (pct.pe_z - 1)))
        ibnd = TRUE;

    if (!ibnd)
        return;


/* precalc some boundaries */
    incx = (dimy + 2) * (dimz + 2);
    incy = (dimz + 2);
    incz = 1;
    xmax = (dimx + 1) * incx;
    ymax = (dimy + 1) * incy;


    /* If the boundary condition flag is SURFACE then we only set the z-boundary */
    if (ct.boundaryflag != SURFACE)
    {

        if (pex == 0)
        {

            for (iy = 0; iy < dimy + 2; iy++)
            {
                for (iz = 0; iz < dimz + 2; iz++)
                {

                    mat[iy * incy + iz] = val;

                }               /* end for */
            }                   /* end for */

        }                       /* end if */


        if (pex == (pct.pe_x - 1))
        {

            for (iy = 0; iy < dimy + 2; iy++)
            {
                for (iz = 0; iz < dimz + 2; iz++)
                {

                    mat[xmax + iy * incy + iz] = val;

                }               /* end for */
            }                   /* end for */

        }                       /* end if */


        if (pey == 0)
        {

            for (ix = 0; ix < dimx + 2; ix++)
            {
                for (iz = 0; iz < dimz + 2; iz++)
                {

                    mat[ix * incx + iz] = val;

                }               /* end for */
            }                   /* end for */

        }                       /* end if */


        if (pey == (pct.pe_y - 1))
        {

            for (ix = 0; ix < dimx + 2; ix++)
            {
                for (iz = 0; iz < dimz + 2; iz++)
                {

                    mat[ix * incx + ymax + iz] = val;

                }               /* end for */
            }                   /* end for */

        }                       /* end if */

    }                           /* end if */


    if (pez == 0)
    {

        for (ix = 0; ix < dimx + 2; ix++)
        {
            for (iy = 0; iy < dimy + 2; iy++)
            {

                mat[ix * incx + iy * incy] = val;

            }                   /* end for */
        }                       /* end for */

    }                           /* end if */


    if (pez == (pct.pe_z - 1))
    {

        for (ix = 0; ix < dimx + 2; ix++)
        {
            for (iy = 0; iy < dimy + 2; iy++)
            {

                mat[ix * incx + iy * incy + dimz + 1] = val;

            }                   /* end for */
        }                       /* end for */

    }                           /* end if */

}                               /* end set_bc */


/******/
