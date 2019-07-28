#include "negf_prototypes.h"
/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
/*


from distributed array FPX0_BASIS * FPY0_BASIS * FPZ0_BASIS
get  global array  get_FNX_GRID() * get_FNY_GRID() * get_FNZ_GRID()

*/

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "main.h"
#include "init_var.h"
#include "LCR.h"


void X_to_distribute_soft (double * global_array, double * distr_array)
{

    int ix, iy, iz, ii;
    int idx2, idx1, incx, incx1, incy, incy1;

    incx = get_FPY0_GRID() * get_FPZ0_GRID();

    ii = get_FPX_OFFSET();


    for (ix = 0; ix < get_FPX0_GRID(); ix++)
    {
        for (iy = 0; iy < get_FPY0_GRID() * get_FPZ0_GRID(); iy++)
        {
            idx1 = ix * incx + iy;
            idx2 = (ix + ii) * incx + iy;
            distr_array[idx1] = global_array[idx2];
        }
    }
}
