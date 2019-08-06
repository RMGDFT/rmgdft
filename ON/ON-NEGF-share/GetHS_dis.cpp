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
#include "LocalObject.h"
#include "blacs.h"
#include "LdaU_on.h"


void mat_global_to_dist(double *a_dist, int *desca, double *a_glob);
void GetHS_dis(LocalObject<double> &Phi, LocalObject<double> &H_Phi, 
        double *vtot_c, double *Hij_dist, double *Bij_dist, double *Kbpsi_mat)
{
    int ione = 1;

    double hxgrid = Rmg_G->get_hxgrid(1);
    double hygrid = Rmg_G->get_hygrid(1);
    double hzgrid = Rmg_G->get_hzgrid(1);
    int pbasis = Rmg_G->get_P0_BASIS(1);

    int order = ct.kohn_sham_fd_order;

    RmgTimer *RT = new RmgTimer("4-get_HS");

    double *Hij_glob = new double[Phi.num_tot * Phi.num_tot];
    double *Bij_glob = new double[Phi.num_tot * Phi.num_tot];

    for (int st1 = 0; st1 < Phi.num_tot * Phi.num_tot; st1++) Hij_glob[st1] = 0.0;
    for (int st1 = 0; st1 < Phi.num_tot * Phi.num_tot; st1++) Bij_glob[st1] = 0.0;


    /* Loop over states */
    /* calculate the H |phi> on this processor and stored in states1[].psiR[] */

    FiniteDiff FD(&Rmg_L);
    RmgTimer *RT0 = new RmgTimer("4-get_HS_dis: apply fd");

    for(int st1 = 0; st1 < Phi.num_thispe; st1++)
    {
        double *a_phi = &Phi.storage_proj[st1 * pbasis];
        double *h_phi = &H_Phi.storage_proj[st1 * pbasis];
        ApplyAOperator (a_phi, h_phi, "Coarse");

        for (int idx = 0; idx < pbasis; idx++)
        {
            h_phi[idx] = a_phi[idx] * vtot_c[idx] - 0.5 * h_phi[idx];
        }
    }

    if(ct.num_ldaU_ions > 0 )
        ldaU_on->app_vhubbard(H_Phi, *Rmg_G);

    delete RT0;

    /* print_sum(pct.psi_size, states1[ct.state_begin].psiR, "states1 sum get_Hij");
     * print_state_sum(states1); 
     */

    RmgTimer *RT1 = new RmgTimer("4-get_HS: orbit_dot_orbit");
    LO_x_LO(Phi, H_Phi, Hij_glob, *Rmg_G);
    LO_x_LO(Phi, Phi, Bij_glob, *Rmg_G);

    delete(RT1);

    RmgTimer *RT3 = new RmgTimer("4-get_HS: Hvnlij");
    GetHvnlij_dis(Hij_glob, Bij_glob, Kbpsi_mat, Phi.num_tot, LocalProj->num_tot);
    delete(RT3);

    int n2 = Phi.num_tot * Phi.num_tot;

    double t1 = (double) (Rmg_G->get_NX_GRID(1) * Rmg_G->get_NY_GRID(1) * Rmg_G->get_NZ_GRID(1));
    double vel = Rmg_L.get_omega() / t1;

    //  dscal (&n2, &vel, Hij_glob, &ione);
    //  dscal (&n2, &vel, Bij_glob, &ione);

    if (pct.gridpe == 0)
    {
        print_matrix(Hij_glob, 6, Phi.num_tot);
        print_matrix(Bij_glob, 6, Phi.num_tot);
    }


    mat_global_to_dist(Hij_dist, pct.desca, Hij_glob);
    mat_global_to_dist(Bij_dist, pct.desca, Bij_glob);
    delete [] Hij_glob;
    delete [] Bij_glob;

    delete(RT);

}




