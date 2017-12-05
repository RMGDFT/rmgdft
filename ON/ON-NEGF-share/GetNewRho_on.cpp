/************************** SVN Revision Information **************************
 **    $Id$    **
 ******************************************************************************/

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>


#include "make_conf.h"
#include "params.h"

#include "rmgtypedefs.h"
#include "typedefs.h"
#include "init_var.h"
#include "RmgTimer.h"

#include "transition.h"
#include "prototypes_on.h"

#include "my_scalapack.h"
#include "blas.h"
#include "RmgParallelFft.h"
#include "BaseGrid.h"



extern "C" void GetNewRho_on_c(STATE * states, double *rho, double *rho_matrix)
{
    GetNewRho_on(states, rho, rho_matrix);
}
void GetNewRho_on(STATE * states, double *rho, double *rho_matrix)
{
    int idx, ione = 1;
    double t2;
    register double tcharge;

    /* for parallel libraries */

    double *psi1, *psi2, scale;
    int i, st1, st2;
    int loop, state_per_proc, num_recv;
    double *rho_temp;
    double tem;
    int global_basis = Rmg_G->get_NX_GRID(1) * Rmg_G->get_NY_GRID(1) * Rmg_G->get_NZ_GRID(1);
    int global_fbasis = get_FPX0_GRID() *get_FPY0_GRID() *get_FPZ0_GRID();

    int st11;

    RmgTimer *RT0 = new RmgTimer("3-get_new_rho");
    state_per_proc = ct.state_per_proc + 2;

    rho_temp = new double[Rmg_G->get_P0_BASIS(1)];


    for (idx = 0; idx < global_basis; idx++)
        rho_global[idx] = 0.;

    RmgTimer *RT1 = new RmgTimer("3-get_new_rho: states in this proc");
    for (st1 = ct.state_begin; st1 < ct.state_end; st1++)
        for (st2 = st1; st2 < ct.state_end; st2++)
        {
            st11 = st1 - ct.state_begin;
            if (st1 == st2)
                scale =  rho_matrix[st11 * ct.num_states + st2];
            if (st1 != st2) scale = 2.0 * rho_matrix[st11 * ct.num_states + st2];
            psi1 = states[st1].psiR;
            psi2 = states[st2].psiR;

            if (state_overlap_or_not[st11 * ct.num_states + st2 ] == 1)
                density_orbit_X_orbit(st1, st2, scale, psi1, psi2,
                        rho_global, 0, states, orbit_overlap_region);

        }

    delete(RT1);

    RmgTimer *RT2 = new RmgTimer("3-get_new_rho: states other proc");


    for (loop = 0; loop < num_sendrecv_loop; loop++)
    {
        num_recv = recv_from[loop * state_per_proc + 1];

        for (i = 0; i < num_recv; i++)
        {
            st2 = recv_from[loop * state_per_proc + i + 2];

            psi2= states[st2].psiR;
            for (st1 = ct.state_begin; st1 < ct.state_end; st1++)
            {
                psi1 = states[st1].psiR;

                st11 = st1 - ct.state_begin;

                if (state_overlap_or_not[st11 * ct.num_states + st2] == 1)
                {
                    psi1 = states[st1].psiR;
                    scale = 2.0 * rho_matrix[st11 * ct.num_states + st2];
                    density_orbit_X_orbit(st1, st2, scale, psi1, psi2,
                            rho_global, 0, states, orbit_overlap_region);
                }
            }
        }

    }                           /* end of loop  */


    delete(RT2);


    RmgTimer *RT3 = new RmgTimer("3-get_new_rho: distribution");
    MPI_Allreduce(MPI_IN_PLACE, rho_global, global_basis, MPI_DOUBLE, MPI_SUM, pct.grid_comm);
    global_to_distribute(rho_global, rho_temp);
    delete(RT3);

    RmgTimer *RT4 = new RmgTimer("3-get_new_rho: interpolation");
    /* Interpolate onto fine grid, result will be stored in rho*/
    switch (ct.interp_flag)
    {
        case PROLONG_INTERPOLATION:
            mg_prolong_MAX10 (rho, rho_temp, get_FPX0_GRID(), get_FPY0_GRID(), get_FPZ0_GRID(), get_PX0_GRID(), get_PY0_GRID(), get_PZ0_GRID(), get_FG_RATIO(), 6);
            break;
        case FFT_INTERPOLATION:
            //FftFilter(work, *coarse_pwaves, ct.cparm, LOW_PASS);  // limit to G-vectors within the inscribed sphere
            FftInterpolation (*Rmg_G, rho_temp, rho, Rmg_G->default_FG_RATIO, ct.sqrt_interpolation);
            break;
        default:
            //Dprintf ("charge interpolation is set to %d", ct.interp_flag);
            printf("\n ct.interp_flag = %d", ct.interp_flag);
            rmg_error_handler (__FILE__, __LINE__, "ct.interp_flag is set to an invalid value.");


    }

    delete [] rho_temp;

    delete(RT4);

    RmgTimer *RT5 = new RmgTimer("3-get_new_rho: augmented");

    RhoAugmented(rho, rho_matrix);

    tem = 0.0;
    for(st1 = 0; st1 < global_fbasis; st1++)
        tem += rho[st1];


    delete(RT5);

    delete(RT0);


}


