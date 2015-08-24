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



void GetHS(STATE * states, STATE * states1, double *vtot_c, double *Hij_00, double *Bij_00)
{
    int idx, st1;
    int maxst, n2;
    STATE *sp;
    int ione = 1;
    int ixx, iyy, izz;
    unsigned int ion, num_orbital_thision, num_proj;
    int ip, iip1;

    double hxgrid, hygrid, hzgrid;

    hxgrid = Rmg_G->get_hxgrid(1);
    hygrid = Rmg_G->get_hygrid(1);
    hzgrid = Rmg_G->get_hzgrid(1);


    RmgTimer *RT0 = new RmgTimer("orbital_comm");
    orbital_comm(states);
    delete(RT0);



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
        ixx = states[st1].ixmax - states[st1].ixmin + 1;
        iyy = states[st1].iymax - states[st1].iymin + 1;
        izz = states[st1].izmax - states[st1].izmin + 1;

        /* Generate 2*V*psi and store it  in orbit_tem */
        genvlocpsi(states[st1].psiR, st1, states1[st1].psiR, vtot_global, states);


        /* A operating on psi stored in orbit_tem */

        /* Eighth-order finite-differenital method for Laplatian operating on psi stored in orbit_tem */
        app10_del2(sp->psiR, orbit_tem, ixx, iyy, izz, hxgrid, hygrid, hzgrid);

        /* A |psi > + 0.5 (B V|psi> + V B |psi>) */

        for (idx = 0; idx < ixx * iyy * izz; idx++)
        {
            states1[st1].psiR[idx] = 0.5 * states1[st1].psiR[idx] - 0.5 * orbit_tem[idx];

        }                       
    }                           /* end for st1 = .. */

    /* print_sum(pct.psi_size, states1[ct.state_begin].psiR, "states1 sum get_Hij");
     * print_state_sum(states1); 
     */

    /* calculate the < states.psiR | states1.psiR>  */

    my_barrier();

    RmgTimer *RT1 = new RmgTimer("4-get_HS: orbit_dot_orbit");
    orbit_dot_orbit(states, states1, Hij_00, Bij_00);
    delete(RT1);


    my_barrier();

    RmgTimer *RT2 = new RmgTimer("4-get_HS: kbpsi");
    get_all_kbpsi(states, states, ion_orbit_overlap_region_nl, projectors, kbpsi);
    delete(RT2);

    for(ion = 0; ion < pct.n_ion_center; ion++)
    {
        Kbpsi_str.kbpsi_ion[ion].clear();
        num_orbital_thision = Kbpsi_str.num_orbital_thision[ion];
        num_proj = pct.prj_per_ion[pct.ionidx[ion]];
        Kbpsi_str.orbital_index[ion].resize(num_orbital_thision);


        // uopdate values from this process
        for(idx = 0; idx < num_orbital_thision; idx++)
        {
            st1 =Kbpsi_str.orbital_index[ion][idx];
            iip1 = (st1-ct.state_begin) * pct.n_ion_center * ct.max_nl + ion * ct.max_nl;
            for(ip = 0; ip < num_proj; ip++)
                Kbpsi_str.kbpsi_ion[ion].emplace_back(kbpsi[iip1 + ip]);
        }


    }

    KbpsiComm();

    RmgTimer *RT3 = new RmgTimer("4-get_HS: Hvnlij");
//    for (st1 = 0; st1 < (ct.state_end-ct.state_begin) * ct.num_states; st1++)
//    {
 //       Hij_00[st1] = 0.;
  //      Bij_00[st1] = 0.;
   // }
    GetHvnlij(Hij_00, Bij_00);

/*
    for (st1 = 0; st1 < (ct.state_end-ct.state_begin) * ct.num_states; st1++)
    {
        Hij_00[st1] = 0.;
        Bij_00[st1] = 0.;
    }
    get_Hvnlij(Hij_00, Bij_00);
 */
    delete(RT3);

    fflush(NULL);
  //  exit(0);

    n2 = (ct.state_end-ct.state_begin) * ct.num_states;

    double t1 = (double) (Rmg_G->get_NX_GRID(1) * Rmg_G->get_NY_GRID(1) * Rmg_G->get_NZ_GRID(1));
    double vel = Rmg_L.get_omega() / t1;

    dscal (&n2, &vel, Hij_00, &ione);
    dscal (&n2, &vel, Bij_00, &ione);

    if (pct.gridpe == 0)
    {
        print_matrix(Hij_00, 5, maxst);
        print_matrix(Bij_00, 5, maxst);
    }

    delete(RT);

}
