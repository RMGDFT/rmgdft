/************************** SVN Revision Information **************************
 **    $Id$    **
 ******************************************************************************/

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>



#include "params.h"

#include "rmgtypedefs.h"
#include "typedefs.h"
#include "init_var.h"
#include "RmgTimer.h"

#include "transition.h"
#include "prototypes_on.h"

#include "blas.h"
#include "RmgParallelFft.h"
#include "BaseGrid.h"

extern std::vector<ORBITAL_PAIR> OrbitalPairs;

extern "C" void GetNewRho_on_c(STATE * states, double *rho, double *rho_matrix)
{
    GetNewRho_on(states, rho, rho_matrix);
}
void GetNewRho_on(STATE * states, double *rho, double *rho_matrix)
{
    int idx;

    /* for parallel libraries */

    double *rho_temp;
    int global_basis = Rmg_G->get_NX_GRID(1) * Rmg_G->get_NY_GRID(1) * Rmg_G->get_NZ_GRID(1);

    int pair;

    RmgTimer *RT0 = new RmgTimer("3-get_new_rho");

    rho_temp = new double[Rmg_G->get_P0_BASIS(1)];


    for (idx = 0; idx < global_basis; idx++)
        rho_global[idx] = 0.;

    RmgTimer *RT1 = new RmgTimer("3-get_new_rho: states in this proc");
#pragma omp parallel private(pair)
    {
        double *rho_global_private = new double[global_basis]();
#pragma omp barrier
#pragma omp for schedule(dynamic) nowait
        for(pair = 0; pair < (int)OrbitalPairs.size(); pair++)
        {
            ORBITAL_PAIR onepair;
            onepair = OrbitalPairs[pair];
            int st1 = onepair.orbital1;
            int st2 = onepair.orbital2;

            double scale;
            int st11 = st1 - ct.state_begin;
            scale =  rho_matrix[st11 * ct.num_states + st2];
            double *psi1 = states[st1].psiR;
            double *psi2 = states[st2].psiR;

            if (state_overlap_or_not[st11 * ct.num_states + st2 ] == 1)
                density_orbit_X_orbit(st1, st2, scale, psi1, psi2,
                        rho_global_private, 0, states, onepair);

        }
#pragma omp critical
        for(int idx = 0;idx < global_basis;idx++)rho_global[idx] += rho_global_private[idx];

        delete [] rho_global_private;
    }

    delete(RT1);


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


    delete(RT5);

    delete(RT0);


}


