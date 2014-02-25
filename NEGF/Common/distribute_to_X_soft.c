/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
/*



from distributed array get_FPX0_GRID() * get_FPY0_GRID() * get_FPZ0_GRID()
get X  array  get_FNX_GRID() * get_FPY0_GRID() * get_FPZ0_GRID()


*/

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "main.h"
#include "init_var.h"
#include "LCR.h"


void distribute_to_X_soft (rmg_double_t * distr_array, rmg_double_t * global_array)
{

    int ix, iy, iz, ii, jj, kk;
    int idx2, idx1, incx, incx1, incy, incy1;

    incx = get_FPY0_GRID() * get_FPZ0_GRID();
    incy = get_FPZ0_GRID();

    ii = get_FPX_OFFSET();

    for (idx1 = 0; idx1 < get_FNX_GRID() * get_FPY0_GRID() * get_FPZ0_GRID(); idx1++)
        global_array[idx1] = 0.0;

    for (ix = 0; ix < get_FPX0_GRID(); ix++)
        for (iy = 0; iy < get_FPY0_GRID() * get_FPZ0_GRID(); iy++)
        {
            idx1 = ix * incx + iy;
            idx2 = (ix + ii) * incx + iy;
            global_array[idx2] = distr_array[idx1];
        }

    idx1 = get_FNX_GRID() * get_FPY0_GRID() * get_FPZ0_GRID();
    global_sums_X (global_array, &idx1);

}
