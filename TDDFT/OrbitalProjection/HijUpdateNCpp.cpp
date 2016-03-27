/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
#include <float.h>
#include <stdio.h>
#include <assert.h>



#include "const.h"
#include "RmgTimer.h"
#include "RmgException.h"
#include "rmgtypedefs.h"
#include "params.h"
#include "typedefs.h"
#include "rmg_error.h"
#include "Kpoint.h"
#include "InputKey.h"
#include "blas.h"
#include "init_var.h"
#include "transition.h"
#include "prototypes_on.h"
#include "Kbpsi.h"


#include "../Headers/common_prototypes.h"
#include "../Headers/common_prototypes1.h"
#include "prototypes_tddft.h"


void HijUpdateNCpp (STATE * states_distribute, double *vtot_diff, double *mat_local, double *A_glob)
{
    int idx, st1, st2, idx1, idx2;
    int st11, st22;
    int maxst, n2;
    STATE *sp;
    int ione = 1;
    double tem, tem1;
    int ixx, iyy, izz;
    char msg[100];
    double *psi, one = 1.0, zero = 0.0;

    int ix, iy,iz;

    maxst = ct.num_states;
    int pbasis = get_P0_BASIS();



    if(pct.num_local_orbit >0)
    {

        psi = new double[pct.num_local_orbit * get_P0_BASIS()];
        //  calculate V|psi>
        for (st1 = 0; st1 < pct.num_local_orbit; st1++)
            for(idx1 = 0; idx1 <get_P0_BASIS(); idx1++)
                psi[st1 * get_P0_BASIS() + idx1] = states_distribute[st1].psiR[idx1] * vtot_diff[idx1];
        // <psi|V|psi>
        dgemm ("T", "N", &pct.num_local_orbit, &pct.num_local_orbit,
                &pbasis, &one, psi, &pbasis,
                states_distribute[0].psiR, &pbasis, &zero, mat_local, &pct.num_local_orbit);

        delete [] psi;

    }


    n2 = pct.num_local_orbit * pct.num_local_orbit;
    double vel = get_vel();
    dscal (&n2, &vel, mat_local, &ione);

    MatrixToGlobal(states_distribute, mat_local, A_glob);


}
