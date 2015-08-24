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
    int nh, ion, ip1, ip2, st1, st2, ist;
    double alpha;
    MPI_Status mstatus;
    int ion1, ion2, ion1_global, ion2_global;
    int iip1, iip2, iip1a, iip2a;
    int size, proc, proc1, proc2, idx;
    double *dnmI, *qnmI;
    double temA, temB;
    double *double_ptr;
    int *int_ptr;
    unsigned int tot_orbital, idx1, idx2;




    double t1 = (double) (Rmg_G->get_NX_GRID(1) * Rmg_G->get_NY_GRID(1) * Rmg_G->get_NZ_GRID(1));
    alpha = t1/Rmg_L.get_omega();

    for (ion = 0; ion < pct.n_ion_center; ion++)
    {
        /* begin shuchun wang */
        ion1 = pct.ionidx[ion];
        nh = pct.prj_per_ion[ion1];
        dnmI = pct.dnmI[ion1];
        qnmI = pct.qqq[ion1];
        double_ptr = Kbpsi_str.kbpsi_ion[ion].data();
        int_ptr = Kbpsi_str.orbital_index[ion].data();
        tot_orbital = Kbpsi_str.orbital_index[ion].size();

        for(idx1 = 0; idx1 < Kbpsi_str.num_orbital_thision[ion]; idx1++)
        {
            st1 = int_ptr[idx1];
            for(idx2 = 0; idx2 < tot_orbital; idx2++)
            {
                st2 = int_ptr[idx2];
                ist = (st1 - ct.state_begin) * ct.num_states + st2;

       //     if(st1 == 0 && st2 == 0) 
       //     {  printf("\n aaa %d %d %d %f  %f", ion, idx1, idx2, double_ptr[idx1 * nh], double_ptr[idx1*nh + 1]);
       //     }

                for (ip1 = 0; ip1 < nh; ip1++)
                {
                    for (ip2 = 0; ip2 < nh; ip2++)
                    {
                        Aij[ist] +=
                            alpha * dnmI[ip2 * nh + ip1] * double_ptr[idx1*nh + ip1] * double_ptr[idx2*nh + ip2];
                        Bij[ist] +=
                            alpha * qnmI[ip2 * nh + ip1] * double_ptr[idx1*nh + ip1] * double_ptr[idx2*nh + ip2];


                    }
                }
            }                 
        }                    


    }

}
