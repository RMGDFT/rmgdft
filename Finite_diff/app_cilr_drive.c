/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

#include "common_prototypes.h"
#include <float.h>
#include <math.h>
#include <stdlib.h>



void app_cilr_driver (rmg_double_t * psi, rmg_double_t * a_psi, rmg_double_t *b_psi, rmg_double_t *vtot_eig_s, 
    int dimx, int dimy, int dimz, rmg_double_t hx, rmg_double_t hy, rmg_double_t hz, int order)
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

