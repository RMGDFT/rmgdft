/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

#include "main.h"
#include <float.h>
#include <math.h>
#include <stdlib.h>



void app_cilr_driver (REAL * psi, REAL * a_psi, REAL *b_psi, REAL *vtot_eig_s, 
    int dimx, int dimy, int dimz, REAL hx, REAL hy, REAL hz, int order)
{


    if(order == APP_CI_FOURTH) {
        app_cilr_fourth(psi, a_psi, b_psi, vtot_eig_s, dimx, dimy, dimz, hx, hy, hz);
        return;
    }
    if(order == APP_CI_SIXTH) {
        app_cilr_sixth(psi, a_psi, b_psi, vtot_eig_s, dimx, dimy, dimz, hx, hy, hz);
        return;
    }
    
    error_handler("APP_CIR order not programmed yet in app_cir_driver.\n");

}

