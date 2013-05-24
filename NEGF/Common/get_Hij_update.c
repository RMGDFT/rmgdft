/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
#include <float.h>
#include <stdio.h>
#include <assert.h>
#include "main.h"


#if USE_DIS_MAT

#include "my_scalapack.h"

#endif


void get_Hij_update (STATE * states, STATE * states1, double *vtot_c, double *Aij)
{
    int idx, st1, st2, idx1, idx2;
    int maxst, n2;
    STATE *sp;
    int ione = 1;
    REAL tem;
    int ixx, iyy, izz;
    char msg[100];
    double time1, time2, time3, time4;

    time1 = my_crtc ();

    n2 = ct.num_states * ct.num_states;
    maxst = ct.num_states;

    my_malloc_init( vtot_global, NX_GRID * NY_GRID * NZ_GRID, REAL );

    distribute_to_global (vtot_c, vtot_global);


    for (st1 = 0; st1 < ct.num_states * ct.num_states; st1++)
        Aij[st1] = 0.;

    /* Loop over states */
    /* calculate the H |phi> on this processor and stored in states1[].psiR[] */

    time3 = my_crtc ();

    for (st1 = ct.state_begin; st1 < ct.state_end; st1++)
    {
        sp = &states[st1];
        ixx = states[st1].ixmax - states[st1].ixmin + 1;
        iyy = states[st1].iymax - states[st1].iymin + 1;
        izz = states[st1].izmax - states[st1].izmin + 1;

        /* Generate 2*V*psi and store it  in orbit_tem */
        genvlocpsi (states[st1].psiR, st1, states1[st1].psiR, vtot_global, states);

/*  sprintf(msg, "STATE1: %d sum", st1);
 *  print_sum(states[st1].size, states1[st1].psiR, msg);
 */



        for (idx = 0; idx < ixx * iyy * izz; idx++)
        {
            states1[st1].psiR[idx] = 0.5 * (states1[st1].psiR[idx]);

        }
    }                           /* end for st1 = .. */

    time4 = my_crtc ();
    rmg_timings (H_psi_TIME, (time4 - time3));

    /* calculate the < states.psiR | states1.psiR>  */

    my_barrier ();
    if (pct.gridpe == 0)
        printf ("\n AAAAAAAAAAA %f", Aij[0]);
    time3 = my_crtc ();
    orbit_dot_orbit (states, states1, work_matrix);
    my_barrier ();
    time4 = my_crtc ();
    rmg_timings (ORBIT_DOT_ORBIT_H, (time4 - time3));

    time3 = my_crtc ();


    get_Hvnlij (Aij);

    time4 = my_crtc ();
    rmg_timings (get_Hnl_TIME, (time4 - time3));


    /* symmetrize the Aij */
    for (st1 = 0; st1 < ct.num_states - 1; st1++)
        for (st2 = st1 + 1; st2 < ct.num_states; st2++)
        {
            idx1 = st1 + st2 * ct.num_states;
            idx2 = st2 + st1 * ct.num_states;
            Aij[idx1] = 0.5 * (Aij[idx1] + Aij[idx2]);
            Aij[idx2] = Aij[idx1];
        }

    global_sums (Aij, &n2, pct.grid_comm);     /* sum up Aij contributions */

    dscal (&n2, &ct.vel, Aij, &ione);

    if (pct.gridpe == 0)
    {
        printf (" matrix Hij\n");
        print_matrix (work_matrix, 5, maxst);
        print_sum (n2, work_matrix, "work_matrix before leaving get_Hij ddd ");
    }

    time2 = my_crtc ();

    rmg_timings (GET_Hij_TIME, (time2 - time1));

    my_free(vtot_global);

}
