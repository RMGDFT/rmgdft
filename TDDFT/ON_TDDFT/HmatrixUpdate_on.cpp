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


void HmatrixUpdate_on(LocalObject<double> &Phi, LocalObject<double> &H_Phi, 
        double *vtot_c, double *Hij_glob)
{
    int ione = 1;

    double hxgrid = Rmg_G->get_hxgrid(1);
    double hygrid = Rmg_G->get_hygrid(1);
    double hzgrid = Rmg_G->get_hzgrid(1);
    int pbasis = Rmg_G->get_P0_BASIS(1);


    RmgTimer *RT = new RmgTimer("4-Hupdate");


    for (int st1 = 0; st1 < Phi.num_tot * Phi.num_tot; st1++) Hij_glob[st1] = 0.0;

    for(int st1 = 0; st1 < Phi.num_thispe; st1++)
    {
        double *a_phi = &Phi.storage_proj[st1 * pbasis];
        double *h_phi = &H_Phi.storage_proj[st1 * pbasis];

        for (int idx = 0; idx < pbasis; idx++)
        {
            h_phi[idx] = a_phi[idx] * vtot_c[idx];
        }
    }

    
    LO_x_LO(Phi, H_Phi, Hij_glob, *Rmg_G);

    delete(RT);

}




