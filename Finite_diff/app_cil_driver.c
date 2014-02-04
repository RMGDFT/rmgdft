/************************** SVN Revision Information **************************
 **    $Id: app_cil_driver.c 1722 2012-05-07 19:38:37Z ebriggs $    **
******************************************************************************/

#include "common_prototypes.h"
#include <float.h>
#include <math.h>
#include <stdlib.h>



rmg_double_t app_cil_driver (rmg_double_t * a, rmg_double_t * b, int dimx, int dimy, int dimz, rmg_double_t gridhx, rmg_double_t gridhy, rmg_double_t gridhz, int order)
{


    if(order == APP_CI_FOURTH) {
        return app_cil_fourth(a, b, dimx, dimy, dimz, gridhx, gridhy, gridhz);
    }
    if(order == APP_CI_SIXTH) {
        return app_cil_sixth(a, b, dimx, dimy, dimz, gridhx, gridhy, gridhz);
    }
    
    error_handler("APP_CIL order not programmed yet in app_cil_driver.\n");

}

rmg_double_t app_cil_driver_f (rmg_float_t * a, rmg_float_t * b, int dimx, int dimy, int dimz, rmg_double_t gridhx, rmg_double_t gridhy, rmg_double_t gridhz, int order)
{


    if(order == APP_CI_FOURTH) {
        return app_cil_fourth_f(a, b, dimx, dimy, dimz, gridhx, gridhy, gridhz);
    }
    if(order == APP_CI_SIXTH) {
        return app_cil_sixth_f(a, b, dimx, dimy, dimz, gridhx, gridhy, gridhz);
    }
    
    error_handler("APP_CIL order not programmed yet in app_cil_driver.\n");

}
