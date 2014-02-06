/************************** SVN Revision Information **************************
 **    $Id: app_cir_sixth.c 1840 2013-01-22 13:52:34Z ebriggs $    **
******************************************************************************/

#include "const.h"
#include "common_prototypes.h"
#include "fixed_dims.h"
#include "rmg_alloc.h"


#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <pthread.h>
#include "hybrid.h"


// Compilers can generate much better code if they know the loop dimensions at compile
// as opposed to run time. Therefore since most of the finite difference stencils
// are applied at the global level we check at the top level to see if the grid
// dimensions correpsond to the global case. If so we call a routine with those
// dimensions set at compile time. If not we just fall through to the general case.

// For the global grid case we also implement an optimized version of the operator
// that takes advantage of certain symmetries to improve performance.

void app_cir_sixth_f (rmg_float_t * a, rmg_float_t * b, int dimx, int dimy, int dimz)
{

    int numgrid, P0_BASIS;
    rmg_float_t *rptr;

    int sbasis = (dimx + 4) * (dimy + 4) * (dimz + 4);

    P0_BASIS = get_P0_BASIS();

    my_malloc (rptr, sbasis + 64, rmg_float_t);

    trade_imagesx_f (a, rptr, dimx, dimy, dimz, 2, FULL_FD);

    // first check for fixed dim case
    numgrid = dimx * dimy * dimz;
    if(numgrid == P0_BASIS) {
          FD_app_cir_sixth_global_rmg_float(rptr, b);
    }
    else {
          FD_app_cir_sixth_standard_rmg_float(rptr, b, dimx, dimy, dimz);
    }

    my_free(rptr);

}

