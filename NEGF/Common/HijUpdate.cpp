/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
#include <float.h>
#include <stdio.h>
#include <assert.h>


#include "params.h"

#include "rmgtypedefs.h"
#include "typedefs.h"
#include "RmgTimer.h"
#include "transition.h"
#include "LCR.h"
#include "prototypes_on.h"
#include "prototypes_negf.h"
#include "init_var.h"



#include "Scalapack.h"
#include "blas.h"
#include "Kbpsi.h"
#include "FiniteDiff.h"

#include "pmo.h"



void HijUpdate (STATE * states, double *vtot_c, double *Aij)
{
    int idx, st1, st2, idx1, idx2;
    int st11, st22;
    int maxst, n2;
    STATE *sp;
    int ione = 1;
    double t1;
    double tem, tem1;
    int ixx, iyy, izz;
    char msg[100];
    double *psi, one = 1.0, zero = 0.0;

    int ix, iy,iz;

    maxst = ct.num_states;
    int pbasis = get_P0_BASIS();

    RmgTimer *RT = new RmgTimer("2-SCF: HijUpdate");

    for (st1 = 0; st1 < ct.num_states * (ct.state_end-ct.state_begin); st1++)
    {
        Hij_00[st1] = 0.;
        Bij_00[st1] = 0.;
    }
    distribute_to_global(vtot_c, vtot_global);

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
        n2 = ixx * iyy * izz;
        t1 = 0.5;
        dscal(&n2, &t1, states1[st1].psiR, &ione);

    }                           /* end for st1 = .. */


    /* calculate the < states.psiR | states1.psiR>  */


    RmgTimer *RT1 = new RmgTimer("4-get_HS: orbit_dot_orbit");
    orbit_dot_orbit(states, states1, Hij_00, Bij_00);
    delete(RT1);


    RmgTimer *RT3 = new RmgTimer("4-get_HS: Hvnlij");
    GetHvnlij(Hij_00, Bij_00);

    delete(RT3);

    fflush(NULL);

    n2 = (ct.state_end-ct.state_begin) * ct.num_states;

    t1 = (double) (Rmg_G->get_NX_GRID(1) * Rmg_G->get_NY_GRID(1) * Rmg_G->get_NZ_GRID(1));
    double vel = Rmg_L.get_omega() / t1;

    dscal (&n2, &vel, Hij_00, &ione);
    dscal (&n2, &vel, Bij_00, &ione);
    row_to_tri_p (lcr[0].Htri, Hij_00, ct.num_blocks, ct.block_dim);

    if (pct.gridpe == 0)
    {
        print_matrix(Hij_00, 5, maxst);
        print_matrix(Bij_00, 5, maxst);
    }

    delete(RT);

}


