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

double FD_app_cil_sixth_standard_rmg_double(rmg_double_t *rptr, rmg_double_t *b, int dimx, int dimy, int dimz, rmg_double_t gridhx, rmg_double_t gridhy, rmg_double_t gridhz);
double FD_app_cil_sixth_standard_rmg_float(rmg_float_t *rptr, rmg_float_t *b, int dimx, int dimy, int dimz, rmg_double_t gridhx, rmg_double_t gridhy, rmg_double_t gridhz);
double FD_app_cil_sixth_global_rmg_double(rmg_double_t *rptr, rmg_double_t *b, rmg_double_t gridhx, rmg_double_t gridhy, rmg_double_t gridhz);
double FD_app_cil_sixth_global_rmg_float(rmg_float_t *rptr, rmg_float_t *b, rmg_double_t gridhx, rmg_double_t gridhy, rmg_double_t gridhz);


// Compilers can generate much better code if they know the loop dimensions at compile
// as opposed to run time. Therefore since most of the finite difference stencils
// are applied at the global level we check at the top level to see if the grid
// dimensions correpsond to the global case. If so we call a routine with those
// dimensions set at compile time. If not we just fall through to the general case.

rmg_double_t app_cil_sixth (rmg_double_t * psi, rmg_double_t * b, int dimx, int dimy, int dimz,
                    rmg_double_t gridhx, rmg_double_t gridhy, rmg_double_t gridhz)
{

    int numgrid, P0_BASIS;
    rmg_double_t cc;
    rmg_double_t *rptr;

    int sbasis = (dimx + 4) * (dimy + 4) * (dimz + 4);

    P0_BASIS = get_P0_BASIS();

    my_malloc (rptr, sbasis + 64, rmg_double_t);

    trade_imagesx (psi, rptr, dimx, dimy, dimz, 2, FULL_FD);

    // first check for fixed dim case  
    numgrid = dimx * dimy * dimz;
    if(numgrid == P0_BASIS) {
        cc = FD_app_cil_sixth_global_rmg_double (rptr, b, gridhx, gridhy, gridhz);
    }
    else {
        cc = FD_app_cil_sixth_standard_rmg_double (rptr, b, dimx, dimy, dimz, gridhx, gridhy, gridhz);
    }

    my_free(rptr);
    return cc;

}


