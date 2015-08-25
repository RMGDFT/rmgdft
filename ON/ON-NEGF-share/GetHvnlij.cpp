/************************** SVN Revision Information **************************
 **    $Id$    **
 ******************************************************************************/

/*
get_Hvnlij:

Get the elements of the Hamiltonian matrix due to the non-local
potential, and add them into Aij.


 */
#include <float.h>
#include <math.h>
#include <stdlib.h>
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



void GetHvnlij(double *Aij, double *Bij)
{
    int nh, ion, st1, st2, ist;
    double alpha, zero = 0.0, one = 1.0;
    int ion1;
    double *dnmI, *qnmI;
    double *temA, *temB;
    double *double_ptr;
    int *int_ptr;
    int num_orb, tot_orb, idx1, idx2, max_state, max_nl;




    double t1 = (double) (Rmg_G->get_NX_GRID(1) * Rmg_G->get_NY_GRID(1) * Rmg_G->get_NZ_GRID(1));
    alpha = t1/Rmg_L.get_omega();

    max_state = 0;
    max_nl = 0;
    for (ion = 0; ion < pct.n_ion_center; ion++)
    {
        nh = pct.prj_per_ion[pct.ionidx[ion]];
        tot_orb = Kbpsi_str.orbital_index[ion].size();
        max_state = std::max(max_state, tot_orb);
        max_nl = std::max(max_nl, nh);
    }

    temA = new double[(ct.state_end-ct.state_begin) * max_nl];
    temB = new double[(ct.state_end-ct.state_begin) * max_state];

    for (ion = 0; ion < pct.n_ion_center; ion++)
    {
        /* begin shuchun wang */
        ion1 = pct.ionidx[ion];
        nh = pct.prj_per_ion[ion1];
        dnmI = pct.dnmI[ion1];
        qnmI = pct.qqq[ion1];
        double_ptr = Kbpsi_str.kbpsi_ion[ion].data();
        int_ptr = Kbpsi_str.orbital_index[ion].data();
        num_orb = Kbpsi_str.num_orbital_thision[ion]; 
        tot_orb = Kbpsi_str.orbital_index[ion].size();


        dgemm ("T", "N", &num_orb, &nh, &nh, &alpha, double_ptr, &nh, dnmI, &nh, &zero, temA, &num_orb);
        dgemm ("N", "N", &num_orb, &tot_orb, &nh, &one, temA, &num_orb, double_ptr, &nh, &zero, temB, &num_orb);




        for(idx1 = 0; idx1 < num_orb; idx1++)
        {
            st1 = int_ptr[idx1];
            for(idx2 = 0; idx2 < tot_orb; idx2++)
            {
                st2 = int_ptr[idx2];
                ist = (st1 - ct.state_begin) * ct.num_states + st2;

                Aij[ist] += temB[idx1 + idx2 * num_orb];
            }

        }

        dgemm ("T", "N", &num_orb, &nh, &nh, &alpha, double_ptr, &nh, qnmI, &nh, &zero, temA, &num_orb);
        dgemm ("N", "N", &num_orb, &tot_orb, &nh, &one, temA, &num_orb, double_ptr, &nh, &zero, temB, &num_orb);

        for(idx1 = 0; idx1 < num_orb; idx1++)
        {
            st1 = int_ptr[idx1];
            for(idx2 = 0; idx2 < tot_orb; idx2++)
            {
                st2 = int_ptr[idx2];
                ist = (st1 - ct.state_begin) * ct.num_states + st2;

                Bij[ist] += temB[idx1 + idx2 * num_orb];
            }

        }

    }

}
