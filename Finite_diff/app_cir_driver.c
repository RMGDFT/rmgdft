/************************** SVN Revision Information **************************
 **    $Id: app_cir_driver.c 1722 2012-05-07 19:38:37Z ebriggs $    **
******************************************************************************/

#include "common_prototypes.h"
#include <float.h>
#include <math.h>
#include <stdlib.h>



void app_cir_driver (rmg_double_t * a, rmg_double_t * b, int dimx, int dimy, int dimz, int order)
{


    if(order == APP_CI_FOURTH) {
        app_cir_fourth(a, b, dimx, dimy, dimz);
        return;
    }
    if(order == APP_CI_SIXTH) {
        app_cir_sixth(a, b, dimx, dimy, dimz);
        return;
    }
    
    error_handler("APP_CIR order not programmed yet in app_cir_driver.\n");

}

void app_cir_driver_f (rmg_float_t * a, rmg_float_t * b, int dimx, int dimy, int dimz, int order)
{


    if(order == APP_CI_FOURTH) {
        app_cir_fourth_f(a, b, dimx, dimy, dimz);
        return;
    }
    if(order == APP_CI_SIXTH) {
        app_cir_sixth_f(a, b, dimx, dimy, dimz);
        return;
    }

    error_handler("APP_CIR order not programmed yet in app_cir_driver.\n");

}

void app_cir_beta_driver (rmg_double_t * a, rmg_double_t * b, int dimx, int dimy, int dimz, int order)
{


    if(order == APP_CI_FOURTH) {
        app_cir_beta_fourth(a, b, dimx, dimy, dimz);
        return;
    }
    if(order == APP_CI_SIXTH) {
        app_cir_beta_sixth(a, b, dimx, dimy, dimz);
        return;
    }
    
    error_handler("APP_CIR order not programmed yet in app_cir_driver.\n");

}
