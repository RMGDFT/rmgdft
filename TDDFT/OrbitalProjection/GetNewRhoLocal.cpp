/************************** SVN Revision Information **************************
 **    $Id: get_new_rho_local.c 3140 2015-08-06 15:48:24Z luw $    **
 ******************************************************************************/

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include "const.h"
#include "RmgTimer.h"
#include "RmgException.h"
#include "rmgtypedefs.h"
#include "params.h"
#include "typedefs.h"
#include "rmg_error.h"
#include "blas.h"
#include "init_var.h"
#include "transition.h"
#include "prototypes_on.h"
#include "Kbpsi.h"


#include "../Headers/common_prototypes.h"
#include "../Headers/common_prototypes1.h"
#include "prototypes_tddft.h"



void GetNewRhoLocal (STATE * states_distribute, double *rho, double *mat_local, double *rho_matrix)
{
    int idx, ione = 1;
    double t2;
    register double tcharge;

    /* for parallel libraries */

    double *psi1, *psi2, scale;
    int i, st1, st2, proc1, proc2, st11;
//    int loop, state_per_proc, num_send, num_recv, num_sendrecv, size1, size2;
    double *rho_temp;
    int ix, iy; 
    double tem;
    char filename[MAX_PATH];

    double *psi, one = 1.0, zero = 0.0;
    int pbasis = get_P0_BASIS();

    rho_temp = new double[pbasis];

    for(idx = 0; idx < pbasis; idx++)rho_temp[idx] = 0.0;

    if(pct.num_local_orbit > 0)
    {
        
        psi = new double[pbasis * pct.num_local_orbit];

        dgemm ("N", "N", &pbasis, &pct.num_local_orbit, &pct.num_local_orbit, &one, 
                states_distribute[0].psiR, &pbasis, mat_local, &pct.num_local_orbit, 
                &zero, psi, &pbasis);

        for(st1 = 0; st1 < pct.num_local_orbit; st1++)
            for(idx = 0; idx < pbasis; idx++)
                rho_temp[idx] += states_distribute[st1].psiR[idx] * psi[st1 * pbasis + idx];

        delete [] psi;
    }


    mg_prolong_MAX10 (rho, rho_temp, get_FPX0_GRID(), get_FPY0_GRID(), get_FPZ0_GRID(),
            get_PX0_GRID(), get_PY0_GRID(), get_PZ0_GRID(), get_FG_RATIO(), 6);


    delete [] rho_temp;
    RhoAugmented(rho, rho_matrix);


    MPI_Barrier(pct.img_comm);

}

