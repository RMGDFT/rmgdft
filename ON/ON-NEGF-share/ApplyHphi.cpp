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


void ApplyHphi(LocalObject<double> &Phi, LocalObject<double> &H_Phi, double *vtot_c)
{

    int pbasis = Rmg_G->get_P0_BASIS(1);

    /* Loop over states */
    /* calculate the H |phi> on this processor and stored in states1[].psiR[] */

    FiniteDiff FD(&Rmg_L);

    for(int st1 = 0; st1 < Phi.num_thispe; st1++)
    {
        double *a_phi = &Phi.storage_cpu[st1 * pbasis];
        double *h_phi = &H_Phi.storage_cpu[st1 * pbasis];
        ApplyAOperator (a_phi, h_phi);

        for (int idx = 0; idx < pbasis; idx++)
        {
            h_phi[idx] = a_phi[idx] * vtot_c[idx] - 0.5 * h_phi[idx];
        }
    }


    //print_sum_square(Phi.num_tot * pbasis, H_Phi.storage_proj, "HPhi before");
    if(ct.num_ldaU_ions > 0 )
        ldaU_on->app_vhubbard(H_Phi, *Rmg_G);
    //print_sum_square(Phi.num_tot * pbasis, H_Phi.storage_proj, "HPhi aftere");


}




