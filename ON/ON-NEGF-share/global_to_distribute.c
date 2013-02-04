/************************** SVN Revision Information **************************
 **    $Id$    **
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


void global_to_distribute(REAL * global_array, REAL * distr_array)
{

    int ix, iy, iz, ii, jj, kk;
    int idx2, idx1, incx, incx1, incy, incy1;

    incx = PY0_GRID * PZ0_GRID;
    incy = PZ0_GRID;
    incx1 = NY_GRID * NZ_GRID;
    incy1 = NZ_GRID;

    pe2xyz(pct.gridpe, &ii, &jj, &kk);
    ii = pct.PX_OFFSET;
    jj = pct.PY_OFFSET;
    kk = pct.PZ_OFFSET;

    idx1 = NX_GRID * NY_GRID * NZ_GRID;

    for (ix = 0; ix < PX0_GRID; ix++)
        for (iy = 0; iy < PY0_GRID; iy++)
            for (iz = 0; iz < PZ0_GRID; iz++)
            {
                idx1 = ix * incx + iy * incy + iz;
                idx2 = (ix + ii) * incx1 + (iy + jj) * incy1 + iz + kk;
                distr_array[idx1] = global_array[idx2];
            }

}
