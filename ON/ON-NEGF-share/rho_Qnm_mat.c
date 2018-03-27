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

#include "main.h"
#include "grid.h"
#include "prototypes_on.h"
#include "init_var.h"


void rho_Qnm_mat(double *Aij, double * global_mat_X, int
*state_begin, int *state_end, int *num_nonlocal_ion, double *kbpsi,
int max_ion_nonlocal, double *kbpsi_comm, int *ionidx_allproc)
{
    int ion, ip1, ip2, st1, st2, ist;
    MPI_Status mstatus;
    int ion1, ion2, ion1_global, ion2_global;
    int iip1, iip2, iip1a, iip2a;
    int size, proc, proc1, proc2, idx;
    int nh;
    int st11, index;
    double tem;



    /* Loop over states on this proce onle 
       (distribution of work AND Aij contributions) */
    proc = pct.gridpe;
    for (st1 = ct.state_begin; st1 < ct.state_end; st1++)
        for (st2 = state_begin[proc]; st2 < state_end[proc]; st2++)
        {
            st11 = st1 - ct.state_begin;
            iip1 = (st1 - state_begin[proc]) * num_nonlocal_ion[proc] * ct.max_nl;
            iip2 = (st2 - state_begin[proc]) * num_nonlocal_ion[proc] * ct.max_nl;

            for (ion = 0; ion < num_nonlocal_ion[proc]; ion++)
            {
                /* begin shuchun wang */
                ion1 = pct.ionidx[ion];
                nh = pct.prj_per_ion[ion1];
                for (ip1 = 0; ip1 < nh; ip1++)
                {
                    for (ip2 = 0; ip2 < nh; ip2++)
                    {
                        if (fabs(kbpsi[iip1 + ip1]) > 0.)
                        {
                            ist = ion1 * ct.max_nl * ct.max_nl + ip1 * ct.max_nl + ip2;
                            Aij[ist] +=
                                global_mat_X[st11 * ct.num_states + st2] *
                                kbpsi[iip1 + ip1] * kbpsi[iip2 + ip2];
                        }
                    }
                }
                iip1 += ct.max_nl;
                iip2 += ct.max_nl;
                /*end shuchun wang */
            }                   /* end for ion */
        }                       /* end for st1 and st2 */


    /* Now calculate the part that kbpsi is stored in other processors */

    size = ct.state_per_proc * max_ion_nonlocal * ct.max_nl;

    void *RT1 = BeginRmgTimer("3-get_new_rho: augmented_comm_loop ");
    for (idx = 0; idx < kbpsi_num_loop; idx++)
    {

        proc1 = kbpsi_comm_send[idx];
        proc2 = kbpsi_comm_recv[idx];

        int tag1 = idx * pct.grid_npes  + pct.gridpe;
        int tag2 = idx * pct.grid_npes  + proc2;
        MPI_Request request;

        if(proc1 >=0)
        {
            MPI_Isend(kbpsi, size, MPI_DOUBLE, proc1, tag1, pct.grid_comm, &request);
        }
        if(proc2 >=0)
            MPI_Recv(kbpsi_comm, size, MPI_DOUBLE, proc2, tag2, pct.grid_comm, &mstatus);
        if(proc1 >=0) MPI_Wait(&request, &mstatus);

        if(proc2 < 0) continue;

      //  dprintf("\n loopaaa %d  proc1 %d  proc2 %d thispe %d", idx, proc1, proc2,
//pct.gridpe);


    void *RT2 = BeginRmgTimer("3-get_new_rho: augmented_comm_loop_calc ");

        for (ion1 = 0; ion1 < num_nonlocal_ion[proc]; ion1++)
            for (ion2 = 0; ion2 < num_nonlocal_ion[proc2]; ion2++)
            {
                ion1_global = ionidx_allproc[proc * max_ion_nonlocal + ion1];
                ion2_global = ionidx_allproc[proc2 * max_ion_nonlocal + ion2];

                if (ion1_global == ion2_global)
                {
                    for (st1 = ct.state_begin; st1 < ct.state_end; st1++)
                    {
                        st11 = st1 - ct.state_begin;
                        index = st11 * ct.num_ions + ion1_global;
                        if(!ion_orbit_overlap_region_nl[index].flag) continue;
                        iip1 = (st1 - state_begin[proc]) * num_nonlocal_ion[proc] * ct.max_nl;

                        for (ip2 = 0; ip2 < ct.max_nl; ip2++)
                        {

                            tem = 0.0;
                            for (st2 = state_begin[proc2]; st2 < state_end[proc2]; st2++)
                            {

                                iip2 = (st2 - state_begin[proc2]) * num_nonlocal_ion[proc2] * ct.max_nl;
                                iip2a = iip2 + ion2 * ct.max_nl + ip2;
                                tem += global_mat_X[st11 * ct.num_states + st2] * kbpsi_comm[iip2a];
                            }

                            for (ip1 = 0; ip1 < ct.max_nl; ip1++)
                            {
                                iip1a = iip1 + ion1 * ct.max_nl + ip1;

                                ist = ion1_global * ct.max_nl *
                                    ct.max_nl + ip1 * ct.max_nl + ip2;
                                Aij[ist] += tem * kbpsi[iip1a];
                            }
                        }   /* end shuchun wang */
                    }       /* end if (ion1_glo... */
                }           /* end for ion1 and ion2 */

            }                   /* end for st1 and st2 */
        EndRmgTimer(RT2);
    }                           /* end for idx */

    EndRmgTimer(RT1);

}

