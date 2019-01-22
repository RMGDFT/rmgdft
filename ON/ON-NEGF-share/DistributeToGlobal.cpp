/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
/*


distribute_to_global.c

from distributed array PX0_BASIS * PY0_BASIS * PZ0_BASIS
get  global array  get_NX_GRID() * get_NY_GRID() * get_NZ_GRID()

*/

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "rmgtypedefs.h"
#include "params.h"
#include "typedefs.h"
#include "RmgTimer.h"
#include "transition.h"
#include "blas.h"
#include "RmgSumAll.h"

#include "prototypes_on.h"
#include "init_var.h"


void DistributeToGlobal(double * distr_array, double * global_array)
{

    int ix, iy, iz, ii, jj, kk;
    int idx2, idx1, incx, incx1, incy, incy1;
    int dimx =  Rmg_G->get_PX0_GRID(1);
    int dimy =  Rmg_G->get_PY0_GRID(1);
    int dimz =  Rmg_G->get_PZ0_GRID(1);
    int global_basis = Rmg_G->get_NX_GRID(1) * Rmg_G->get_NY_GRID(1) * Rmg_G->get_NZ_GRID(1);

    incx = Rmg_G->get_PY0_GRID(1) * Rmg_G->get_PZ0_GRID(1);
    incy = Rmg_G->get_PZ0_GRID(1);
    incx1 = Rmg_G->get_NY_GRID(1) * Rmg_G->get_NZ_GRID(1);
    incy1 = Rmg_G->get_NZ_GRID(1);

    ii =  Rmg_G->get_PX_OFFSET(1);
    jj =  Rmg_G->get_PY_OFFSET(1);
    kk =  Rmg_G->get_PZ_OFFSET(1);

    for (idx1 = 0; idx1 < global_basis; idx1++)
        global_array[idx1] = 0.0;

    for (ix = 0; ix < dimx; ix++)
        for (iy = 0; iy < dimy; iy++)
            for (iz = 0; iz < dimz; iz++)
            {
                idx1 = ix * incx + iy * incy + iz;
                idx2 = (ix + ii) * incx1 + (iy + jj) * incy1 + iz + kk;
                global_array[idx2] = distr_array[idx1];
            }

    // global_sums(global_array, &global_basis, pct.grid_comm);
    MPI_Allreduce(MPI_IN_PLACE, global_array, global_basis, MPI_DOUBLE, MPI_SUM, pct.grid_comm);


}
