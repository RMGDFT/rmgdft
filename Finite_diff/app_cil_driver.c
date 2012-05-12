/************************** SVN Revision Information **************************
 **    $Id: app_cil_driver.c 1722 2012-05-07 19:38:37Z ebriggs $    **
******************************************************************************/

#include "main.h"
#include <float.h>
#include <math.h>
#include <stdlib.h>



REAL app_cil_driver (REAL * a, REAL * b, int dimx, int dimy, int dimz, REAL gridhx, REAL gridhy, REAL gridhz, int order)
{


    if(order == APP_CI_FOURTH) {
        return app_cil_fourth(a, b, dimx, dimy, dimz, gridhx, gridhy, gridhz);
    }
    if(order == APP_CI_SIXTH) {
        return app_cil_sixth(a, b, dimx, dimy, dimz, gridhx, gridhy, gridhz);
    }
    
    error_handler("APP_CIL order not programmed yet in app_cil_driver.\n");

}
