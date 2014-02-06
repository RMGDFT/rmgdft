/************************** SVN Revision Information **************************
 **    $Id$    **
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

void app_cir_sixth (rmg_double_t * a, rmg_double_t * b, int dimx, int dimy, int dimz)
{

    int numgrid, P0_BASIS;
    rmg_double_t *rptr=NULL;

    int sbasis = (dimx + 4) * (dimy + 4) * (dimz + 4);

    P0_BASIS = get_P0_BASIS();

    my_malloc (rptr, sbasis + 64, rmg_double_t);

    trade_imagesx (a, rptr, dimx, dimy, dimz, 2, FULL_FD);

    // first check for fixed dim case
    numgrid = dimx * dimy * dimz;
    if(numgrid == P0_BASIS) {
//        app_cir_sixth_global (rptr, b);
        FD_app_cir_sixth_global_rmg_double(rptr, b);

    }
    else {
//        app_cir_sixth_standard (rptr, b, dimx, dimy, dimz);
        FD_app_cir_sixth_standard_rmg_double(rptr, b, dimx, dimy, dimz);
    }

    my_free(rptr);

}

