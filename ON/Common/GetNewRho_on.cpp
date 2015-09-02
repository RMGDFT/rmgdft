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

    int st11;

    RmgTimer *RT0 = new RmgTimer("3-get_new_rho");
    state_per_proc = ct.state_per_proc + 2;

    rho_temp = new double[Rmg_G->get_P0_BASIS(1)];


    for (idx = 0; idx < Rmg_G->get_NX_GRID(1) * Rmg_G->get_NY_GRID(1) * Rmg_G->get_NZ_GRID(1); idx++)
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
    idx = get_NX_GRID() * get_NY_GRID() * get_NZ_GRID();
    global_sums(rho_global, &idx, pct.grid_comm);

    global_to_distribute(rho_global, rho_temp);

    delete(RT3);

    RmgTimer *RT4 = new RmgTimer("3-get_new_rho: interpolation");
    mg_prolong_MAX10 (rho, rho_temp, get_FPX0_GRID(), get_FPY0_GRID(), get_FPZ0_GRID(), get_PX0_GRID(), get_PY0_GRID(), get_PZ0_GRID(), get_FG_RATIO(), 6);

    delete [] rho_temp;

    delete(RT4);

    RmgTimer *RT5 = new RmgTimer("3-get_new_rho: augmented");

    RhoAugmented(rho, rho_matrix);

    delete(RT5);

    int iii = get_FP0_BASIS();

    tcharge = 0.0;
    for (idx = 0; idx < get_FP0_BASIS(); idx++)
        tcharge += rho[idx];
    ct.tcharge = real_sum_all(tcharge, pct.grid_comm);
    ct.tcharge = real_sum_all(ct.tcharge, pct.spin_comm);


    ct.tcharge *= get_vel_f();

    t2 = ct.nel / ct.tcharge;
    dscal(&iii, &t2, &rho[0], &ione);


    if(fabs(t2 -1.0) > 1.0e-6 && pct.gridpe == 0)
        printf("\n Warning: total charge Normalization constant = %e  \n", t2);

    delete(RT0);


}
