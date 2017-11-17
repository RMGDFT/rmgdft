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


void distribute_to_global(double * distr_array, double * global_array)
{

    int ix, iy, iz, ii, jj, kk;
    int idx2, idx1, incx, incx1, incy, incy1;
    int dimx =  get_PX0_GRID();
    int dimy =  get_PY0_GRID();
    int dimz =  get_PZ0_GRID();
    int global_basis = get_NX_GRID() * get_NY_GRID() * get_NZ_GRID();

    incx = get_PY0_GRID() * get_PZ0_GRID();
    incy = get_PZ0_GRID();
    incx1 = get_NY_GRID() * get_NZ_GRID();
    incy1 = get_NZ_GRID();

    ii =  get_PX_OFFSET();
    jj =  get_PY_OFFSET();
    kk =  get_PZ_OFFSET();

    for (idx1 = 0; idx1 < global_basis; idx1++)
        global_array[idx1] = 0.0;

    for (ix = 0; ix < dimx; ix++)
        for (iy = 0; iy < dimy; iy++)
            for (iz = 0; iz < dimz; iz++)
            {
                idx1 = ix * incx + iy * incy + iz;
                idx2 = (ix + ii) * incx1 + (iy + jj) * incy1 + iz + kk;
                global_array[idx2] = distr_array[idx1];
            }

    global_sums(global_array, &global_basis, pct.grid_comm);

}
