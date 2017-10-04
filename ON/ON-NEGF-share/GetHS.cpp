/************************** SVN Revision Information **************************
 **    $Id: 
******************************************************************************/

/* calculating Hamiltonian matrix Hij and 
 * overlap matrix matB togather
 */

 
#include <float.h>
#include <stdio.h>
#include <assert.h>

#include "make_conf.h"
#include "params.h"

#include "rmgtypedefs.h"
#include "typedefs.h"
#include "RmgTimer.h"
#include "transition.h"
#include "prototypes_on.h"
#include "init_var.h"



#include "my_scalapack.h"
#include "blas.h"
#include "Kbpsi.h"
#include "FiniteDiff.h"


extern BaseGrid *OG;
extern FiniteDiff *OFD;


void GetHS(STATE * states, STATE * states1, double *vtot_c, double *Hij_00, double *Bij_00)
{
    int idx, st1;
    int maxst, n2;
    STATE *sp;
    int ione = 1;
    unsigned int ion, num_orbital_thision, num_proj;
    int ip, iip1;

    double hxgrid, hygrid, hzgrid;

    int order = ct.kohn_sham_fd_order;

    int ixx = states[0].ixmax - states[0].ixmin + 1;
    int iyy = states[0].iymax - states[0].iymin + 1;
    int izz = states[0].izmax - states[0].izmin + 1;
    if(OG == NULL) OG = new BaseGrid(ixx, iyy, izz, 1, 1, 1, 0, 1);
    if(OFD == NULL) OFD = new FiniteDiff(&Rmg_L, OG, CLUSTER, CLUSTER, CLUSTER, 1, order);


    int item = (ct.max_orbit_nx+order) *(ct.max_orbit_ny+order) *(ct.max_orbit_nz+order);
    hxgrid = Rmg_G->get_hxgrid(1);
    hygrid = Rmg_G->get_hygrid(1);
    hzgrid = Rmg_G->get_hzgrid(1);



    RmgTimer *RT = new RmgTimer("4-get_HS");


    maxst = ct.num_states;

    distribute_to_global(vtot_c, vtot_global);

    for (st1 = 0; st1 < (ct.state_end-ct.state_begin) * ct.num_states; st1++)
    {
        Hij_00[st1] = 0.;
        Bij_00[st1] = 0.;
    }

    /* Loop over states */
    /* calculate the H |phi> on this processor and stored in states1[].psiR[] */

    for (st1 = ct.state_begin; st1 < ct.state_end; st1++)
    {
        sp = &states[st1];

        /* Generate 2*V*psi and store it  in orbit_tem */
        genvlocpsi(states[st1].psiR, st1, states1[st1].psiR, vtot_global, states);


        /* A operating on psi stored in orbit_tem */

        /* Eighth-order finite-differenital method for Laplatian operating on psi stored in orbit_tem */
        OFD->app_del2_np(sp->psiR, orbit_tem, hxgrid, hygrid, hzgrid);


        /* A |psi > + 0.5 (B V|psi> + V B |psi>) */

        //ZeroBoundary(orbit_tem, ixx, iyy, izz);
        for (idx = 0; idx < ixx * iyy * izz; idx++)
        {
            states1[st1].psiR[idx] = 0.5 * states1[st1].psiR[idx] - 0.5 * orbit_tem[idx];

        }                       
        ZeroBoundary(states1[st1].psiR, ixx, iyy, izz);

    }                           /* end for st1 = .. */

    /* print_sum(pct.psi_size, states1[ct.state_begin].psiR, "states1 sum get_Hij");
     * print_state_sum(states1); 
     */

    /* calculate the < states.psiR | states1.psiR>  */


    RmgTimer *RT1 = new RmgTimer("4-get_HS: orbit_dot_orbit");
    orbit_dot_orbit(states, states1, Hij_00, Bij_00);
    delete(RT1);


    RmgTimer *RT3 = new RmgTimer("4-get_HS: Hvnlij");
    GetHvnlij(Hij_00, Bij_00);
    delete(RT3);


    n2 = (ct.state_end-ct.state_begin) * ct.num_states;

    double t1 = (double) (Rmg_G->get_NX_GRID(1) * Rmg_G->get_NY_GRID(1) * Rmg_G->get_NZ_GRID(1));
    double vel = Rmg_L.get_omega() / t1;

    dscal (&n2, &vel, Hij_00, &ione);
    dscal (&n2, &vel, Bij_00, &ione);

    if (pct.gridpe == 0)
    {
        print_matrix(Hij_00, 5, maxst);
        print_matrix(Bij_00, 4, maxst);
    }

    delete(RT);

}


