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

#include "main.h"
#include "prototypes_on.h"
#include "init_var.h"


void get_Hvnlij(double *Aij, double *Bij)
{
    int nh, ion, ip1, ip2, st1, st2, ist;
    double alpha;
    MPI_Status mstatus;
    int ion1, ion2, ion1_global, ion2_global;
    int iip1, iip2, iip1a, iip2a;
    int size, proc, proc1, proc2, idx;
    double *dnmI, *qnmI;
    double temA, temB;




    alpha = 1. / get_vel();


    /* Loop over states on this proce onle 
       (distribution of work AND Aij contributions) */
    proc = pct.gridpe;
    for (st1 = ct.state_begin; st1 < ct.state_end; st1++)
        for (st2 = ct.state_begin; st2 < ct.state_end; st2++)
        {

            ist = (st1 - ct.state_begin) * ct.num_states + st2;
            iip1 = (st1 - state_begin[proc]) * num_nonlocal_ion[proc] * ct.max_nl;
            iip2 = (st2 - state_begin[proc]) * num_nonlocal_ion[proc] * ct.max_nl;
            for (ion = 0; ion < num_nonlocal_ion[proc]; ion++)
            {
                /* begin shuchun wang */
                ion1 = pct.ionidx[ion];
                nh = pct.prj_per_ion[ion1];
                dnmI = pct.dnmI[ion1];
                qnmI = pct.qqq[ion1];
                for (ip1 = 0; ip1 < nh; ip1++)
                {
                    for (ip2 = 0; ip2 < nh; ip2++)
                    {
                        if (fabs(kbpsi[iip1 + ip1]) > 0.)
                        {
                            Aij[ist] +=
                                alpha * dnmI[ip2 * nh + ip1] * kbpsi[iip1 +
                                ip1] * kbpsi[iip2 + ip2];
                            Bij[ist] +=
                                alpha * qnmI[ip2 * nh + ip1] * kbpsi[iip1 +
                                ip1] * kbpsi[iip2 + ip2];
                        }


                    }
                }
                iip1 += ct.max_nl;
                iip2 += ct.max_nl;
                /*end shuchun wang */
            }                   /* end for ion */
        }                       /* end for st1 and st2 */

    /*  print_sum(ct.num_states* ct.num_states, Aij, "get_Hij vnnn");
     */


    /* Now calculate the part that kbpsi is stored in other processors */

    size = ct.state_per_proc * max_ion_nonlocal * ct.max_nl;

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


        for (ion1 = 0; ion1 < num_nonlocal_ion[proc]; ion1++)
            for (ion2 = 0; ion2 < num_nonlocal_ion[proc2]; ion2++)
            {
                ion1_global = ionidx_allproc[proc * max_ion_nonlocal + ion1];
                ion2_global = ionidx_allproc[proc2 * max_ion_nonlocal + ion2];

                if (ion1_global == ion2_global)
                {

                    /* begin shuchun wang */
                    nh = pct.prj_per_ion[ion1_global];
                    dnmI = pct.dnmI[ion1_global];
                    qnmI = pct.qqq[ion1_global];

                    for (st1 = ct.state_begin; st1 < ct.state_end; st1++)
                    {
                        iip1 = (st1 - state_begin[proc]) * num_nonlocal_ion[proc] * ct.max_nl;

                        for (ip2 = 0; ip2 < nh; ip2++)
                        {

                            temA = 0.0;
                            temB = 0.0;
                            for (ip1 = 0; ip1 < nh; ip1++)
                            {
                                iip1a = iip1 + ion1 * ct.max_nl + ip1;
                                temA += alpha * dnmI[ip2 * nh + ip1] * kbpsi[iip1a] ;
                                temB += alpha * qnmI[ip2 * nh + ip1] * kbpsi[iip1a] ;
                            }

                            for (st2 = state_begin[proc2]; st2 < state_end[proc2]; st2++)
                            {

                                ist = (st1 - ct.state_begin) * ct.num_states + st2;
                                iip2 = (st2 - state_begin[proc2]) * num_nonlocal_ion[proc2] * ct.max_nl;

                                iip2a = iip2 + ion2 * ct.max_nl + ip2;

                                Aij[ist] += temA * kbpsi_comm[iip2a];
                                Bij[ist] += temB * kbpsi_comm[iip2a];
                            }
                        }   /* end shuchun wang */
                    }       /* end if (ion1_glo... */
                }           /* end for ion1 and ion2 */

            }                   /* end for st1 and st2 */
    }                           /* end for idx */




}
