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
#include "RmgException.h"


void mat_global_to_dist(double *a_dist, int *desca, double *a_glob);
void GetHS_proj(LocalObject<double> &Phi1, LocalObject<double> &Phi2, LocalObject<double> &H_Phi, 
        double *vtot_c, double *Hij_glob, double *Bij_glob, double *Kbpsi_mat1, double *Kbpsi_mat2)
{

    int pbasis = Rmg_G->get_P0_BASIS(1);

    // Phi2 and H_Phi need to have the same data structure
    if(Phi2.num_tot != H_Phi.num_tot || Phi2.num_thispe != H_Phi.num_thispe)
    {
        throw RmgFatalException() << "H_Phi and Phi2 have different data structure"  << __FILE__ << " at line " << __LINE__ << "\n";
    }

    RmgTimer *RT = new RmgTimer("4-get_HS");
    for (int st1 = 0; st1 < Phi1.num_tot * Phi2.num_tot; st1++) Hij_glob[st1] = 0.0;
    for (int st1 = 0; st1 < Phi1.num_tot * Phi2.num_tot; st1++) Bij_glob[st1] = 0.0;


    /* Loop over states */
    /* calculate the H |phi> on this processor and stored in states1[].psiR[] */

    FiniteDiff FD(&Rmg_L);
    RmgTimer *RT0 = new RmgTimer("4-get_HS_proj: apply fd");

    for(int st1 = 0; st1 < Phi2.num_thispe; st1++)
    {
        double *a_phi = &Phi2.storage_proj[st1 * pbasis];
        double *h_phi = &H_Phi.storage_proj[st1 * pbasis];
        ApplyAOperator (a_phi, h_phi, "Coarse");

        for (int idx = 0; idx < pbasis; idx++)
        {
            h_phi[idx] = a_phi[idx] * vtot_c[idx] - 0.5 * h_phi[idx];
        }
    }


    //print_sum_square(Phi.num_tot * pbasis, H_Phi.storage_proj, "HPhi before");
    if(ct.num_ldaU_ions > 0 )
        ldaU_on->app_vhubbard(H_Phi, *Rmg_G);
    //print_sum_square(Phi.num_tot * pbasis, H_Phi.storage_proj, "HPhi aftere");

    delete RT0;

    /* print_sum(pct.psi_size, states1[ct.state_begin].psiR, "states1 sum get_Hij");
     * print_state_sum(states1); 
     */

    RmgTimer *RT1 = new RmgTimer("4-get_HS: orbit_dot_orbit");
    LO_x_LO(Phi1, H_Phi, Hij_glob, *Rmg_G);
    LO_x_LO(Phi1, Phi2, Bij_glob, *Rmg_G);

    delete(RT1);


    RmgTimer *RT3 = new RmgTimer("4-get_HS: Hvnlij");
    GetHvnlij_proj(Hij_glob, Bij_glob, Kbpsi_mat1, Kbpsi_mat2, Phi1.num_tot, Phi2.num_tot, LocalProj->num_tot);
    delete(RT3);


    if (pct.gridpe == 0)
    {
        print_matrix(Hij_glob, 4, Phi1.num_tot);
        print_matrix(Bij_glob, 4, Phi1.num_tot);
    }

    delete(RT);

}




