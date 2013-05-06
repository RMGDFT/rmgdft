/************************** SVN Revision Information **************************
 **    $Id$    **
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
#include "main.h"


void X_to_distribute_soft (REAL * global_array, REAL * distr_array)
{

    int ix, iy, iz, ii;
    int idx2, idx1, incx, incx1, incy, incy1;

    incx = FPY0_GRID * FPZ0_GRID;

    ii = pct.FPX_OFFSET;


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
