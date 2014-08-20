/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

#include "const.h"
#include "common_prototypes.h"
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include "FiniteDiff.h"




void app_cilr_driver (double * psi, double * a_psi, double *b_psi, double *vtot_eig_s, 
    int dimx, int dimy, int dimz, double hx, double hy, double hz, int order)
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

