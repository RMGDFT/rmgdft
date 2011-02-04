/************************** SVN Revision Information **************************
 **    $Id: distribute_to_X_soft.c 1242 2011-02-02 18:55:23Z luw $    **
******************************************************************************/
 
/*



from distributed array FPX0_GRID * FPY0_GRID * FPZ0_GRID
get X  array  FNX_GRID * FPY0_GRID * FPZ0_GRID


*/

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "md.h"


void distribute_to_X_soft (REAL * distr_array, REAL * global_array)
{

    int ix, iy, iz, ii, jj, kk;
    int idx2, idx1, incx, incx1, incy, incy1;

    incx = FPY0_GRID * FPZ0_GRID;
    incy = FPZ0_GRID;

    pe2xyz (pct.thispe, &ii, &jj, &kk);
    ii *= FPX0_GRID;

    for (idx1 = 0; idx1 < FNX_GRID * FPY0_GRID * FPZ0_GRID; idx1++)
        global_array[idx1] = 0.0;

    for (ix = 0; ix < FPX0_GRID; ix++)
        for (iy = 0; iy < FPY0_GRID * FPZ0_GRID; iy++)
        {
            idx1 = ix * incx + iy;
            idx2 = (ix + ii) * incx + iy;
            global_array[idx2] = distr_array[idx1];
        }

    idx1 = FNX_GRID * FPY0_GRID * FPZ0_GRID;
    global_sums_X (global_array, &idx1);

}
