/************************** SVN Revision Information **************************
 **    $Id: 
******************************************************************************/

/* calculating Hamiltonian matrix Hij and 
 * overlap matrix matB togather
 */

 
#include <float.h>
#include <stdio.h>
#include <assert.h>

#include "params.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "RmgTimer.h"
#include "transition.h"
#include "prototypes_on.h"
#include "init_var.h"
#include "blas.h"
#include "Kbpsi.h"
#include "FiniteDiff.h"





void GetHS(STATE * states, STATE * states1, double *vtot_c, double *Hij_00, double *Bij_00)
{
    int ione = 1;
    int maxst = ct.num_states;

    double hxgrid = Rmg_G->get_hxgrid(1);
    double hygrid = Rmg_G->get_hygrid(1);
    double hzgrid = Rmg_G->get_hzgrid(1);

    int order = ct.kohn_sham_fd_order;
    int item = (ct.max_orbit_nx+order) *(ct.max_orbit_ny+order) *(ct.max_orbit_nz+order);


    RmgTimer *RT = new RmgTimer("4-get_HS");
    distribute_to_global(vtot_c, vtot_global);

    for (int st1 = 0; st1 < (ct.state_end-ct.state_begin) * ct.num_states; st1++) Hij_00[st1] = 0.0;
    for (int st1 = 0; st1 < (ct.state_end-ct.state_begin) * ct.num_states; st1++) Bij_00[st1] = 0.0;


    /* Loop over states */
    /* calculate the H |phi> on this processor and stored in states1[].psiR[] */

    FiniteDiff FD(&Rmg_L);
    int st1;
    double *orbital_border, *orbit_tem;
    RmgTimer *RT0 = new RmgTimer("4-get_HS: apply fd");
#pragma omp parallel private(st1,orbital_border,orbit_tem)
{
    orbital_border = new double[2*item];
    orbit_tem = new double[2*item];
#pragma omp for schedule(static,1) nowait
    for (st1 = ct.state_begin; st1 < ct.state_end; st1++)
    {
        STATE *sp = &states[st1];
        int ixx = states[st1].orbit_nx;
        int iyy = states[st1].orbit_ny;
        int izz = states[st1].orbit_nz;

        /* Generate 2*V*psi and store it  in orbit_tem */
        genvlocpsi(states[st1].psiR, st1, states1[st1].psiR, vtot_global, states);

        /* A operating on psi stored in orbit_tem */
        /* Eighth-order finite-differencing for Laplacian operating on psi stored in orbit_tem */
        if(sp->radius > 0) 
        {
            FillOrbitalBorders(orbital_border, sp->psiR, ixx, iyy, izz, order);
        }
        else
        {
            Rmg_T->trade_imagesx_central_local(sp->psiR, orbital_border, ixx, iyy, izz, order/2);
        }

        if(ct.laplacian_offdiag || ct.laplacian_autocoeff)
            FiniteDiffLap (orbital_border, orbit_tem, ixx, iyy, izz, LC);
        else
            FD.app8_del2 (orbital_border, orbit_tem, ixx, iyy, izz, hxgrid, hygrid, hzgrid);
        /* A |psi > + 0.5 (B V|psi> + V B |psi>) */

        for (int idx = 0; idx < ixx * iyy * izz; idx++)
        {
            states1[st1].psiR[idx] = 0.5 * states1[st1].psiR[idx] - 0.5 * orbit_tem[idx];
        }                       

    }                           /* end for st1 = .. */

    delete [] orbit_tem;
    delete [] orbital_border;
}
    delete RT0;

    /* print_sum(pct.psi_size, states1[ct.state_begin].psiR, "states1 sum get_Hij");
     * print_state_sum(states1); 
     */

    /* calculate the < states.psiR | states1.psiR>  */


    RmgTimer *RT1 = new RmgTimer("4-get_HS: orbit_dot_orbit");
    OrbitDotOrbit(states, states1, Hij_00, Bij_00);
    delete(RT1);


    RmgTimer *RT3 = new RmgTimer("4-get_HS: Hvnlij");
    GetHvnlij(Hij_00, Bij_00);
    delete(RT3);


    int n2 = (ct.state_end-ct.state_begin) * ct.num_states;

    double t1 = (double) (Rmg_G->get_NX_GRID(1) * Rmg_G->get_NY_GRID(1) * Rmg_G->get_NZ_GRID(1));
    double vel = Rmg_L.get_omega() / t1;

    dscal (&n2, &vel, Hij_00, &ione);
    dscal (&n2, &vel, Bij_00, &ione);

    if (pct.gridpe == 0)
    {
        print_matrix(Hij_00, 6, maxst);
        print_matrix(Bij_00, 6, maxst);
    }

    delete(RT);

}


