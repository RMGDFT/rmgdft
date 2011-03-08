/************************** SVN Revision Information **************************
 **    $Id$    **
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



void fill_orbit_borders3(REAL * sg, REAL * pg, int dimx, int dimy, int dimz)
{

    int ix, iy, iz, ixh, iyh;
    int incx, incy, incxs, incys;
    int ione = 1;
    incys = dimz;
    incxs = dimy * dimz;
    incx = (dimy + 6) * (dimz + 6);
    incy = dimz + 6;


    for (ix = 0; ix < dimx + 6; ix++)
        for (iy = 0; iy < dimy + 6; iy++)
            for (iz = 0; iz < dimz + 6; iz++)
                sg[ix * incx + iy * incy + iz] = 0.0;

/* Load up original values from pg  */

    for (ix = 0; ix < dimx; ix++)
        for (iy = 0; iy < dimy; iy++)
            for (iz = 0; iz < dimz; iz++)
                sg[(ix + 3) * incx + (iy + 3) * incy + iz + 3] = pg[ix * incxs + iy * incys + iz];


}                               /* end fill_orbit_borders3.c */

/******/
