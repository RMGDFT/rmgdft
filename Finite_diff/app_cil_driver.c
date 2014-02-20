/************************** SVN Revision Information **************************
 **    $Id: app_cil_driver.c 1722 2012-05-07 19:38:37Z ebriggs $    **
******************************************************************************/

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include "const.h"
#include "common_prototypes.h"
#include "rmg_alloc.h"
#include "rmg_error.h"


double FD_app_cil_sixth_standard_rmg_double(rmg_double_t *rptr, rmg_double_t *b, int dimx, int dimy, int dimz, rmg_double_t gridhx, rmg_double_t gridhy, rmg_double_t gridhz);
double FD_app_cil_sixth_standard_rmg_float(rmg_float_t *rptr, rmg_float_t *b, int dimx, int dimy, int dimz, rmg_double_t gridhx, rmg_double_t gridhy, rmg_double_t gridhz);
double FD_app_cil_sixth_global_rmg_double(rmg_double_t *rptr, rmg_double_t *b, rmg_double_t gridhx, rmg_double_t gridhy, rmg_double_t gridhz);
double FD_app_cil_sixth_global_rmg_float(rmg_float_t *rptr, rmg_float_t *b, rmg_double_t gridhx, rmg_double_t gridhy, rmg_double_t gridhz);


rmg_double_t app_cil_driver (rmg_double_t * a, rmg_double_t * b, int dimx, int dimy, int dimz, rmg_double_t gridhx, rmg_double_t gridhy, rmg_double_t gridhz, int order)
{

    int P0_BASIS, numgrid, sbasis;
    rmg_double_t cc;
    rmg_double_t *rptr;

    P0_BASIS = get_P0_BASIS();
    numgrid = dimx * dimy * dimz;

    sbasis = (dimx + 4) * (dimy + 4) * (dimz + 4);
    my_malloc (rptr, sbasis + 64, rmg_double_t);

    if(order == APP_CI_FOURTH) {
        my_free(rptr);
        return app_cil_fourth(a, b, dimx, dimy, dimz, gridhx, gridhy, gridhz);
    }
    if(order == APP_CI_SIXTH) {
        trade_imagesx (a, rptr, dimx, dimy, dimz, 2, FULL_FD);
        if(numgrid == P0_BASIS) {
            cc = FD_app_cil_sixth_global_rmg_double (rptr, b, gridhx, gridhy, gridhz);
        }
        else {
            cc = FD_app_cil_sixth_standard_rmg_double (rptr, b, dimx, dimy, dimz, gridhx, gridhy, gridhz);
        }

        my_free(rptr);
        return cc;
    }
    
    rmg_error_handler("APP_CIL order not programmed yet in app_cil_driver.\n");

}

rmg_double_t app_cil_driver_f (rmg_float_t * a, rmg_float_t * b, int dimx, int dimy, int dimz, rmg_double_t gridhx, rmg_double_t gridhy, rmg_double_t gridhz, int order)
{

    int P0_BASIS, numgrid, sbasis;
    rmg_double_t cc;
    rmg_float_t *rptr;

    P0_BASIS = get_P0_BASIS();
    numgrid = dimx * dimy * dimz;

    sbasis = (dimx + 4) * (dimy + 4) * (dimz + 4);
    my_malloc (rptr, sbasis + 64, rmg_double_t);

    if(order == APP_CI_FOURTH) {
        my_free(rptr);
        return app_cil_fourth_f(a, b, dimx, dimy, dimz, gridhx, gridhy, gridhz);
    }
    if(order == APP_CI_SIXTH) {
        trade_imagesx_f (a, rptr, dimx, dimy, dimz, 2, FULL_FD);
        if(numgrid == P0_BASIS) {
            cc = FD_app_cil_sixth_global_rmg_float (rptr, b, gridhx, gridhy, gridhz);
        }
        else {
            cc = FD_app_cil_sixth_standard_rmg_float (rptr, b, dimx, dimy, dimz, gridhx, gridhy, gridhz);
        }
        my_free(rptr);
        return cc;
    }
    
    rmg_error_handler("APP_CIL order not programmed yet in app_cil_driver.\n");

}
