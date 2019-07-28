#include "negf_prototypes.h"
/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
/*
In case of Fine grid, this routine helps to move data 
from the global array get_FNX_GRID() * get_FNY_GRID() 
to the distributed array get_FPX0_GRID() * get_FPY0_GRID() 
*/

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "main.h"
#include "init_var.h"
#include "LCR.h"


void global_to_distribute3 (double *global_array, double *distr_array)
{

    int ix, iy, ii, jj, kk;
    int idx2, idx1, incy, incy1;

    incy = get_FPY0_GRID();
    incy1 = get_FNY_GRID();

    ii = get_FPX_OFFSET();
    jj = get_FPY_OFFSET();

    for (idx1 = 0; idx1 < get_FPX0_GRID() * get_FPY0_GRID(); idx1++)
        distr_array[idx1] = 0.0;

    for (ix = 0; ix < get_FPX0_GRID(); ix++)
        for (iy = 0; iy < get_FPY0_GRID(); iy++)
            {
                idx1 = ix * incy + iy;
                idx2 = (ix + ii) * incy1 + (iy + jj);
                distr_array[idx1] = global_array[idx2];
            }

}
