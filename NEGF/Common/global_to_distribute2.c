/************************** SVN Revision Information **************************
 **    $Id: global_to_distribute2.c 834 2007-08-07 17:32:44Z luw $    **
******************************************************************************/
 
/*
In case of Fine grid, this routine helps to move data 
from the global array FNX_GRID * FNY_GRID * FNZ_GRID 
to the distributed array FPX0_GRID * FPY0_GRID * FPZ0_GRID
*/

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "md.h"


void global_to_distribute2 (double *global_array, double *distr_array)
{

    int ix, iy, iz, ii, jj, kk;
    int idx2, idx1, incx, incx1, incy, incy1;

    incx = FPY0_GRID * FPZ0_GRID;
    incy = FPZ0_GRID;
    incx1 = FNY_GRID * FNZ_GRID;
    incy1 = FNZ_GRID;

    pe2xyz (pct.thispe, &ii, &jj, &kk);
    ii *= FPX0_GRID;
    jj *= FPY0_GRID;
    kk *= FPZ0_GRID;

    for (idx1 = 0; idx1 < FP0_BASIS; idx1++)
        distr_array[idx1] = 0.0;

    for (ix = 0; ix < FPX0_GRID; ix++)
        for (iy = 0; iy < FPY0_GRID; iy++)
            for (iz = 0; iz < FPZ0_GRID; iz++)
            {
                idx1 = ix * incx + iy * incy + iz;
                idx2 = (ix + ii) * incx1 + (iy + jj) * incy1 + iz + kk;
                distr_array[idx1] = global_array[idx2];
            }

}
