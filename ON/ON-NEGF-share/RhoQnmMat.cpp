/************************** SVN Revision Information **************************
 **    $Id$    **
 ******************************************************************************/

/*
rho_Qnm_mat:

Get the elements of the charge density due to the augmented function
and add them into Aij.


 */
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>


#include "params.h"

#include "rmgtypedefs.h"
#include "typedefs.h"
#include "RmgTimer.h"

#include "prototypes_on.h"
#include "init_var.h"
#include "Kbpsi.h"

// Aij = <beta_i| phi_k1> * X_k1,k2 * <beta_j|phi_k2>


void RhoQnmMat(double *Aij, double * global_mat_X)
{
    int ion, ip1, ip2, st1, st2, ist;
    int ion1;
    int nh;
    int st11;

    int idx1, idx2, num_orb, tot_orb;

    double *mat_kbpsi = new double[ct.max_nl];

    for (ion = 0; ion < pct.n_ion_center; ion++)
    {
        ion1 = pct.ionidx[ion];
        nh = pct.prj_per_ion[ion1];
        num_orb = Kbpsi_str.num_orbital_thision[ion];
        tot_orb = Kbpsi_str.orbital_index[ion].size();

        for(idx1 = 0; idx1 < num_orb; idx1++)
        {
            st1 = Kbpsi_str.orbital_index[ion][idx1];
            st11 = st1 - ct.state_begin;

            for(ip1 = 0; ip1 < nh; ip1++) mat_kbpsi[ip1] = 0.0;

            for(idx2 = 0; idx2 < tot_orb; idx2++)
            {
                st2 = Kbpsi_str.orbital_index[ion][idx2];

                for (ip2 = 0; ip2 < nh; ip2++)
                {
                    mat_kbpsi[ip2] += global_mat_X[st11 * ct.num_states + st2] *
                        Kbpsi_str.kbpsi_ion[ion][idx2 * nh + ip2];
                }
            }


            for (ip1 = 0; ip1 < nh; ip1++)
            {
                for (ip2 = 0; ip2 < nh; ip2++)
                {
                    ist = ion1 * ct.max_nl * ct.max_nl + ip1 * ct.max_nl + ip2;
                    Aij[ist] += mat_kbpsi[ip2] *
                        Kbpsi_str.kbpsi_ion[ion][idx1 * nh + ip1];
                }
            }
        }             
    }                

    delete [] mat_kbpsi;
}
