/*

  Driver routine that calls the correct version of mg_eig_state

*/

#include "main.h"

void mg_eig_state_driver (STATE * sp, int tid, REAL * vtot_psi, int precision)
{

#if !GAMMA_PT

    // Single precision code path not programmed for non gamma calculations
    mg_eig_state (sp, tid, vtot_psi);

#else

    if(precision == sizeof(rmg_double_t)) {

        mg_eig_state (sp, tid, vtot_psi);

    }
    else {

        mg_eig_state_f (sp, tid, vtot_psi);

    }

#endif

}
