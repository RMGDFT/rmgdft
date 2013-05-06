/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
/*
In case of Fine grid, this routine helps to move data 
from the global array FNX_GRID * FNY_GRID 
to the distributed array FPX0_GRID * FPY0_GRID 
*/

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "main.h"


void global_to_distribute3 (double *global_array, double *distr_array)
{

    int ix, iy, ii, jj, kk;
    int idx2, idx1, incy, incy1;

    incy = FPY0_GRID;
    incy1 = FNY_GRID;

    ii = pct.FPX_OFFSET;
    jj = pct.FPY_OFFSET;

    for (idx1 = 0; idx1 < FPX0_GRID * FPY0_GRID; idx1++)
        distr_array[idx1] = 0.0;

    for (ix = 0; ix < FPX0_GRID; ix++)
        for (iy = 0; iy < FPY0_GRID; iy++)
            {
                idx1 = ix * incy + iy;
                idx2 = (ix + ii) * incy1 + (iy + jj);
                distr_array[idx1] = global_array[idx2];
            }

}
