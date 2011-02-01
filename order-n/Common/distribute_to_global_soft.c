/************************** SVN Revision Information **************************
 **    $Id: distribute_to_global_soft.c 587 2006-08-18 22:57:56Z miro $    **
******************************************************************************/
 
/*


distribute_to_global_soft.c

from distributed array FPX0_BASIS * FPY0_BASIS * FPZ0_BASIS
get  global array  FNX_GRID * FNY_GRID * FNZ_GRID

*/

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "md.h"


void distribute_to_global_soft(REAL * distr_array, REAL * global_array)
{

    int ix, iy, iz, ii, jj, kk;
    int idx2, idx1, incx, incx1, incy, incy1;

    incx = FPY0_GRID * FPZ0_GRID;
    incy = FPZ0_GRID;
    incx1 = FNY_GRID * FNZ_GRID;
    incy1 = FNZ_GRID;

    pe2xyz(pct.thispe, &ii, &jj, &kk);
    ii *= FPX0_GRID;
    jj *= FPY0_GRID;
    kk *= FPZ0_GRID;

    for (idx1 = 0; idx1 < FNX_GRID * FNY_GRID * FNZ_GRID; idx1++)
        global_array[idx1] = 0.0;

    for (ix = 0; ix < FPX0_GRID; ix++)
        for (iy = 0; iy < FPY0_GRID; iy++)
            for (iz = 0; iz < FPZ0_GRID; iz++)
            {
                idx1 = ix * incx + iy * incy + iz;
                idx2 = (ix + ii) * incx1 + (iy + jj) * incy1 + iz + kk;
                global_array[idx2] = distr_array[idx1];
            }

    idx1 = FNX_GRID * FNY_GRID * FNZ_GRID;
    global_sums(global_array, &idx1);

}
