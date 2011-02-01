/************************** SVN Revision Information **************************
 **    $Id: distribute_to_global.c 587 2006-08-18 22:57:56Z miro $    **
******************************************************************************/
 
/*


distribute_to_global.c

from distributed array PX0_BASIS * PY0_BASIS * PZ0_BASIS
get  global array  NX_GRID * NY_GRID * NZ_GRID

*/

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "md.h"


void distribute_to_global(REAL * distr_array, REAL * global_array)
{

    int ix, iy, iz, ii, jj, kk;
    int idx2, idx1, incx, incx1, incy, incy1;

    incx = PY0_GRID * PZ0_GRID;
    incy = PZ0_GRID;
    incx1 = NY_GRID * NZ_GRID;
    incy1 = NZ_GRID;

    pe2xyz(pct.thispe, &ii, &jj, &kk);
    ii *= PX0_GRID;
    jj *= PY0_GRID;
    kk *= PZ0_GRID;

    for (idx1 = 0; idx1 < NX_GRID * NY_GRID * NZ_GRID; idx1++)
        global_array[idx1] = 0.0;

    for (ix = 0; ix < PX0_GRID; ix++)
        for (iy = 0; iy < PY0_GRID; iy++)
            for (iz = 0; iz < PZ0_GRID; iz++)
            {
                idx1 = ix * incx + iy * incy + iz;
                idx2 = (ix + ii) * incx1 + (iy + jj) * incy1 + iz + kk;
                global_array[idx2] = distr_array[idx1];
            }

    idx1 = NX_GRID * NY_GRID * NZ_GRID;
    global_sums(global_array, &idx1);

}
