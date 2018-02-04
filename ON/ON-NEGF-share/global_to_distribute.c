/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
/*


distribute_to_global.c

from distributed array PX0_BASIS * PY0_BASIS * PZ0_BASIS
get  global array  get_NX_GRID() * get_NY_GRID() * get_NZ_GRID()

*/

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "main.h"
#include "prototypes_on.h"


void global_to_distribute(double * global_array, double * distr_array)
{

    int ix, iy, iz, ii, jj, kk;
    int idx2, idx1, incx, incx1, incy, incy1;
    int dimx = get_PX0_GRID();
    int dimy = get_PY0_GRID();
    int dimz = get_PZ0_GRID();

    incx = dimy * dimz;
    incy = dimz;
    incx1 = get_NY_GRID() * get_NZ_GRID();
    incy1 = get_NZ_GRID();

    pe2xyz(pct.gridpe, &ii, &jj, &kk);
    ii = get_PX_OFFSET();
    jj = get_PY_OFFSET();
    kk = get_PZ_OFFSET();

    idx1 = get_NX_GRID() * get_NY_GRID() * get_NZ_GRID();

    for (ix = 0; ix < dimx; ix++)
        for (iy = 0; iy < dimy; iy++)
            for (iz = 0; iz < dimz; iz++)
            {
                idx1 = ix * incx + iy * incy + iz;
                idx2 = (ix + ii) * incx1 + (iy + jj) * incy1 + iz + kk;
                distr_array[idx1] = global_array[idx2];
            }

}
