/************************** SVN Revision Information **************************
 **    $Id: X_to_distribute_soft.c 1242 2011-02-02 18:55:23Z luw $    **
******************************************************************************/
 
/*


from distributed array FPX0_BASIS * FPY0_BASIS * FPZ0_BASIS
get  global array  FNX_GRID * FNY_GRID * FNZ_GRID

*/

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "md.h"


void X_to_distribute_soft (REAL * global_array, REAL * distr_array)
{

    int ix, iy, iz, ii, jj, kk;
    int idx2, idx1, incx, incx1, incy, incy1;

    incx = FPY0_GRID * FPZ0_GRID;

    pe2xyz (pct.thispe, &ii, &jj, &kk);
    ii *= FPX0_GRID;


    for (ix = 0; ix < FPX0_GRID; ix++)
    {
        for (iy = 0; iy < FPY0_GRID * FPZ0_GRID; iy++)
        {
            idx1 = ix * incx + iy;
            idx2 = (ix + ii) * incx + iy;
            distr_array[idx1] = global_array[idx2];
        }
    }
}
