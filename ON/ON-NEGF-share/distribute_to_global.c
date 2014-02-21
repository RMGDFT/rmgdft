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
#include "main_on.h"


void distribute_to_global(rmg_double_t * distr_array, rmg_double_t * global_array)
{

    int ix, iy, iz, ii, jj, kk;
    int idx2, idx1, incx, incx1, incy, incy1;

    incx = get_PY0_GRID() * get_PZ0_GRID();
    incy = get_PZ0_GRID();
    incx1 = get_NY_GRID() * get_NZ_GRID();
    incy1 = get_NZ_GRID();

    ii =  pct.PX_OFFSET;
    jj =  pct.PY_OFFSET;
    kk =  pct.PZ_OFFSET;

    for (idx1 = 0; idx1 < get_NX_GRID() * get_NY_GRID() * get_NZ_GRID(); idx1++)
        global_array[idx1] = 0.0;

    for (ix = 0; ix < get_PX0_GRID(); ix++)
        for (iy = 0; iy < get_PY0_GRID(); iy++)
            for (iz = 0; iz < get_PZ0_GRID(); iz++)
            {
                idx1 = ix * incx + iy * incy + iz;
                idx2 = (ix + ii) * incx1 + (iy + jj) * incy1 + iz + kk;
                global_array[idx2] = distr_array[idx1];
            }

    idx1 = get_NX_GRID() * get_NY_GRID() * get_NZ_GRID();
    global_sums(global_array, &idx1, pct.grid_comm);

}
