/************************** SVN Revision Information **************************
 **    $Id: app_cir_driver.c 1722 2012-05-07 19:38:37Z ebriggs $    **
******************************************************************************/

#include "const.h"
#include "common_prototypes.h"
#include "main.h"
#include "rmg_alloc.h"


#include <float.h>
#include <math.h>
#include <stdlib.h>
#include "FiniteDiff.h"



void app_cir_beta_driver (double * a, double * b, int dimx, int dimy, int dimz, int order)
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
