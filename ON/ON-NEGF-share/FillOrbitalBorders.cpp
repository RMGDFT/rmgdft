/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
#include <float.h>
#include <math.h>
#include <assert.h>



void FillOrbitalBorders(double * sg, double * pg, int dimx, int dimy, int dimz, int order)
{

    int ix, iy, iz;
    int incx, incy, incxs, incys;
    incys = dimz;
    incxs = dimy * dimz;
    incx = (dimy + order) * (dimz + order);
    incy = dimz + order;

    assert(order/2 * 2 == order);

    for (ix = 0; ix < dimx + order; ix++)
        for (iy = 0; iy < dimy + order; iy++)
            for (iz = 0; iz < dimz + order; iz++)
                sg[ix * incx + iy * incy + iz] = 0.0;

/* Load up original values from pg  */

    for (ix = 0; ix < dimx; ix++)
        for (iy = 0; iy < dimy; iy++)
            for (iz = 0; iz < dimz; iz++)
                sg[(ix + order/2) * incx + (iy + order/2) * incy + iz + order/2] = pg[ix * incxs + iy * incys + iz];


}                               /* end fill_orbit_borders3.c */

/******/
