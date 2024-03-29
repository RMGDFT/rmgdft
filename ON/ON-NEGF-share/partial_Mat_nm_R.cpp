/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
/*
    partial_rho_nm_mat:
    
    Get the elements of partial_rho_nm/partial_R for each ion
    and add them into Aij.


*/
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include "main.h"
#include "prototypes_on.h"
#include "init_var.h"


void partial_Mat_nm_R(double *partial_x, double *partial_y, double *partial_z, double * global_mat_X)
{
    int ion, ip1, ip2, st1, st2, ist;
    MPI_Status mstatus;
    int ion1, ion2, ion1_global, ion2_global;
    int iip1, iip2, iip1a, iip2a;
    int size, proc, proc1, proc2, idx;
    int nh;
    double *kbpsi_comm_x, *kbpsi_comm_y, *kbpsi_comm_z;
    int st11;



    size = ct.state_per_proc * max_ion_nonlocal * ct.max_nl;
    my_malloc_init( kbpsi_comm_x, 3 * size, double );
    kbpsi_comm_y = kbpsi_comm_x + size;
    kbpsi_comm_z = kbpsi_comm_y + size;


    /* Loop over states on this proce onle */
    proc = pct.gridpe;
    for (st1 = ct.state_begin; st1 < ct.state_end; st1++)
        for (st2 = ct.state_begin; st2 < ct.state_end; st2++)
        {

            st11 = st1 - ct.state_begin;
            iip1 = (st1 - ct.state_begin) * num_nonlocal_ion[proc] * ct.max_nl;
            iip2 = (st2 - ct.state_begin) * num_nonlocal_ion[proc] * ct.max_nl;

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
                            partial_x[ist] +=  global_mat_X[st11 * ct.num_states + st2] *
                                (kbpsi[iip1 + ip1] * partial_kbpsi_x[iip2 + ip2] +
                                 partial_kbpsi_x[iip1 + ip1] * kbpsi[iip2 + ip2]);
                            partial_y[ist] +=  global_mat_X[st11 * ct.num_states + st2] *
                                (kbpsi[iip1 + ip1] * partial_kbpsi_y[iip2 + ip2] +
                                 partial_kbpsi_y[iip1 + ip1] * kbpsi[iip2 + ip2]);
                            partial_z[ist] +=  global_mat_X[st11 * ct.num_states + st2] *
                                (kbpsi[iip1 + ip1] * partial_kbpsi_z[iip2 + ip2] +
                                 partial_kbpsi_z[iip1 + ip1] * kbpsi[iip2 + ip2]);
                        }
                    }
                }
                iip1 += ct.max_nl;
                iip2 += ct.max_nl;
                /*end shuchun wang */
            }                   /* end for ion */
        }                       /* end for st1 and st2 */


    /* Now calculate the part that kbpsi is stored in other processors */


    for (idx = 1; idx < pct.grid_npes; idx++)
    {

        proc1 = pct.gridpe + idx;
        if (proc1 >= pct.grid_npes)
            proc1 = proc1 - pct.grid_npes;
        proc2 = pct.gridpe - idx;
        if (proc2 < 0)
            proc2 += pct.grid_npes;


        MPI_Sendrecv(kbpsi, size, MPI_DOUBLE, proc1, idx,
                     kbpsi_comm, size, MPI_DOUBLE, proc2, idx, pct.grid_comm, &mstatus);
        MPI_Sendrecv(partial_kbpsi_x, size, MPI_DOUBLE, proc1, idx,
                     kbpsi_comm_x, size, MPI_DOUBLE, proc2, idx, pct.grid_comm, &mstatus);
        MPI_Sendrecv(partial_kbpsi_y, size, MPI_DOUBLE, proc1, idx,
                     kbpsi_comm_y, size, MPI_DOUBLE, proc2, idx, pct.grid_comm, &mstatus);
        MPI_Sendrecv(partial_kbpsi_z, size, MPI_DOUBLE, proc1, idx,
                     kbpsi_comm_z, size, MPI_DOUBLE, proc2, idx, pct.grid_comm, &mstatus);


        for (st1 = ct.state_begin; st1 < ct.state_end; st1++)
            for (st2 = state_begin[proc2]; st2 < state_end[proc2]; st2++)
            {

                st11 = st1 - ct.state_begin;
                iip1 = (st1 - ct.state_begin) * num_nonlocal_ion[proc] * ct.max_nl;
                iip2 = (st2 - state_begin[proc2]) * num_nonlocal_ion[proc2] * ct.max_nl;

                for (ion1 = 0; ion1 < num_nonlocal_ion[proc]; ion1++)
                    for (ion2 = 0; ion2 < num_nonlocal_ion[proc2]; ion2++)
                    {
                        ion1_global = ionidx_allproc[proc * max_ion_nonlocal + ion1];
                        ion2_global = ionidx_allproc[proc2 * max_ion_nonlocal + ion2];

                        if (ion1_global == ion2_global)
                        {

                            /* begin shuchun wang */
                            for (ip1 = 0; ip1 < ct.max_nl; ip1++)
                            {
                                iip1a = iip1 + ion1 * ct.max_nl + ip1;
                                for (ip2 = 0; ip2 < ct.max_nl; ip2++)
                                {
                                    iip2a = iip2 + ion2 * ct.max_nl + ip2;
                                    if (fabs(kbpsi[iip1a]) > 0.)
                                    {
                                        ist =
                                            ion1_global * ct.max_nl * ct.max_nl + ip1 * ct.max_nl +
                                            ip2;
                                        partial_x[ist] +=
                                             global_mat_X[st11 * ct.num_states +
                                                                      st2] * (kbpsi[iip1a] *
                                                                              kbpsi_comm_x[iip2a] +
                                                                              partial_kbpsi_x[iip1a]
                                                                              * kbpsi_comm[iip2a]);

                                        partial_y[ist] +=
                                             global_mat_X[st11 * ct.num_states +
                                                                      st2] * (kbpsi[iip1a] *
                                                                              kbpsi_comm_y[iip2a] +
                                                                              partial_kbpsi_y[iip1a]
                                                                              * kbpsi_comm[iip2a]);

                                        partial_z[ist] +=
                                             global_mat_X[st11 * ct.num_states +
                                                                      st2] * (kbpsi[iip1a] *
                                                                              kbpsi_comm_z[iip2a] +
                                                                              partial_kbpsi_z[iip1a]
                                                                              * kbpsi_comm[iip2a]);
                                    }
                                }
                            }   /* end shuchun wang */
                        }
                    }
            }
    }

    my_free(kbpsi_comm_x);



}
