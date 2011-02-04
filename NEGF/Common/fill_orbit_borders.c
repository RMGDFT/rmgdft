/************************** SVN Revision Information **************************
 **    $Id: fill_orbit_borders.c 1242 2011-02-02 18:55:23Z luw $    **
******************************************************************************/
 
/****f* QMD-MGDFT/pack_ptos.c *****
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
 *   void fill_orbit_borders(REAL *sg, int dimx, int dimy, int dimz)
 *   fill the border values for a orbit 
 * INPUTS
 *   sg[(dimx+2)*(dimy+2)*(dimz+2)]: its elements except for borders are copied from pg
 *               by calling pack_ptos(...)
 *   dimx, dimy, dimz: dimensions of the array
 * OUTPUT
 *   sg[(dimx+2)*(dimy+2)*(dimz+2)]: its border elements are filled
 * PARENTS
 *   
 * CHILDREN
 *   nothing
 * SOURCE
 */


#include "md.h"
#include <float.h>
#include <math.h>



void fill_orbit_borders (REAL * sg, int dimx, int dimy, int dimz)
{

#if 0
    int ix, iy, iz;
    int incx, incy;
    incy = dimz + 2;
    incx = (dimy + 2) * (dimz + 2);


    if (ct.ions[0].xfold == 1)
    {
        for (iy = 1; iy < dimy + 1; iy++)
            for (iz = 1; iz < dimz + 1; iz++)
            {
                ix = 0;
                sg[ix * incx + iy * incy + iz] = sg[dimx * incx + iy * incy + iz];
                sg[(dimx + 1) * incx + iy * incy + iz] = sg[(ix + 1) * incx + iy * incy + iz];
            }
    }

    if (ct.ions[0].yfold == 1)
    {
        for (ix = 0; ix < dimx + 2; ix++)
            for (iz = 1; iz < dimz + 1; iz++)
            {
                iy = 0;
                sg[ix * incx + iy * incy + iz] = sg[ix * incx + dimy * incy + iz];
                sg[ix * incx + (dimy + 1) * incy + iz] = sg[ix * incx + (iy + 1) * incy + iz];
            }
    }

    if (ct.ions[0].zfold == 1)
    {
        for (ix = 0; ix < dimx + 2; ix++)
            for (iy = 0; iy < dimy + 2; iy++)
            {
                iz = 0;
                sg[ix * incx + iy * incy + iz] = sg[ix * incx + iy * incy + dimz];
                sg[ix * incx + iy * incy + dimz + 1] = sg[ix * incx + iy * incy + iz + 1];
            }
    }

#endif

}                               /* end fill_orbit_borders.c */

/******/
