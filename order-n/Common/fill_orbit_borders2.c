/************************** SVN Revision Information **************************
 **    $Id: fill_orbit_borders2.c 587 2006-08-18 22:57:56Z miro $    **
******************************************************************************/
 
/****f* QMD-MGDFT/pack_ptos.c *****
 * NAME
 *   Ab initio real space code with multigrid acceleration
 *   Quantum molecular dynamics package.
 *   Version: 2.1.5
 * COPYRIGHT
 * FUNCTION
 *   void fill_orbit_borders2(REAL *sg, int dimx, int dimy, int dimz)
 *   fill the border values for a orbit 
 * INPUTS
 *   sg[(dimx+4)*(dimy+4)*(dimz+4)]: its elements except for borders are copied from pg
 *               by calling pack_ptos(...)
 *   dimx, dimy, dimz: dimensions of the array
 * OUTPUT
 *   sg[(dimx+4)*(dimy+4)*(dimz+4)]: its border elements are filled
 * PARENTS
 *   
 * CHILDREN
 *   nothing
 * SOURCE
 */


#include "md.h"
#include <float.h>
#include <math.h>



void fill_orbit_borders2(REAL * sg, REAL * pg, int dimx, int dimy, int dimz)
{

#if  0
    int ix, iy, iz, ixh, iyh;
    int incx, incy, incxs, incys;
    int ione = 1;
    incys = dimz;
    incxs = dimy * dimz;
    incx = (dimy + 4) * (dimz + 4);
    incy = dimz + 4;



/* Load up original values from pg  */

    for (ix = 0; ix < dimx; ix++)
        for (iy = 0; iy < dimy; iy++)
            for (iz = 0; iz < dimz; iz++)
                sg[(ix + 2) * incx + (iy + 2) * incy + iz + 2] = pg[ix * incxs + iy * incys + iz];


/* fill the border value */

    if (ct.ions[0].xfold == 1)
    {
        for (iy = 2; iy < dimy + 2; iy++)
            for (iz = 2; iz < dimz + 2; iz++)
            {
                sg[0 * incx + iy * incy + iz] = sg[dimx * incx + iy * incy + iz];
                sg[1 * incx + iy * incy + iz] = sg[(dimx + 1) * incx + iy * incy + iz];
                sg[(dimx + 2) * incx + iy * incy + iz] = sg[2 * incx + iy * incy + iz];
                sg[(dimx + 3) * incx + iy * incy + iz] = sg[3 * incx + iy * incy + iz];


            }
    }

    if (ct.ions[0].yfold == 1)
    {
        for (ix = 0; ix < dimx + 4; ix++)
            for (iz = 2; iz < dimz + 2; iz++)
            {
                sg[ix * incx + 0 * incy + iz] = sg[ix * incx + dimy * incy + iz];
                sg[ix * incx + 1 * incy + iz] = sg[ix * incx + (dimy + 1) * incy + iz];
                sg[ix * incx + (dimy + 2) * incy + iz] = sg[ix * incx + 2 * incy + iz];
                sg[ix * incx + (dimy + 3) * incy + iz] = sg[ix * incx + 3 * incy + iz];
            }
    }

    if (ct.ions[0].zfold == 1)
    {
        for (ix = 0; ix < dimx + 4; ix++)
            for (iy = 0; iy < dimy + 4; iy++)
            {
                sg[ix * incx + iy * incy + 0] = sg[ix * incx + iy * incy + dimz];
                sg[ix * incx + iy * incy + 1] = sg[ix * incx + iy * incy + dimz + 1];
                sg[ix * incx + iy * incy + dimz + 2] = sg[ix * incx + iy * incy + 2];
                sg[ix * incx + iy * incy + dimz + 3] = sg[ix * incx + iy * incy + 3];
            }
    }


#endif


}                               /* end fill_orbit_borders2.c */

/******/
