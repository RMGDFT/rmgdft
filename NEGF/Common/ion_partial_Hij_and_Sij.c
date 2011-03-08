/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include "method.h"
#include "md.h"


void ion_partial_Hij_and_Sij (int iion, int flag,  double *Hij, double *Sij)
{
    int nh, ion, ip1, ip2, st1, st2, ist;
    MPI_Status mstatus;
    int ion1, ion2, ion1_global, ion2_global;
    int iip1, iip2, iip1a, iip2a;
    int size, proc, proc1, proc2, idx, idx1, idx2;
    REAL *dnmI_R, *dnmI, *qqq, temp;
    REAL alpha, *partial_kbpsi, *partial_kbpsi_comm;


    double time1, time2;
    time1 = my_crtc ();

    for(idx = 0; idx < ct.num_states * ct.num_states; idx++) 
    {
        Hij[idx] = Sij[idx] = 0.0;
    }

    alpha = 1.0 / ct.vel;

    size = ct.state_per_proc * max_ion_nonlocal * ct.max_nl;
    my_malloc_init( partial_kbpsi_comm, size, REAL );

    switch (flag)
    {
        case 1:
            partial_kbpsi = partial_kbpsi_x;
            break;
        case 2:
            partial_kbpsi = partial_kbpsi_y;
            break;
        case 3:
            partial_kbpsi = partial_kbpsi_z;
            break;
        default:
            error_handler("Undefined flag");
    }


    /* Loop over states on this proce onle */
    proc = pct.thispe;
    for (st1 = state_begin[proc]; st1 < state_end[proc]; st1++)
        for (st2 = state_begin[proc]; st2 < state_end[proc]; st2++)
        {

            ist = st1 * ct.num_states + st2;
            iip1 = (st1 - state_begin[proc]) * num_nonlocal_ion[proc] * ct.max_nl;
            iip2 = (st2 - state_begin[proc]) * num_nonlocal_ion[proc] * ct.max_nl;

            for (ion = 0; ion < num_nonlocal_ion[proc]; ion++)
            {

                ion1 = pct.ionidx[ion];
                if (ion1 == iion)
                {
                    nh = pct.prj_per_ion[ion1];

                    qqq = pct.qqq[ion1];
                    dnmI = pct.dnmI[ion1];
                    switch (flag)
                    {
                        case 1:
                            dnmI_R = pct.dnmI_x[ion1];
                            break;
                        case 2:
                            dnmI_R = pct.dnmI_y[ion1];
                            break;
                        case 3:
                            dnmI_R = pct.dnmI_z[ion1];
                            break;
                        default:
                            error_handler("Undefined flag");
                    }
                    for (ip1 = 0; ip1 < nh; ip1++)
                    {
                        for (ip2 = 0; ip2 < nh; ip2++)
                        {

                            /****** from sum_nm ( par_Dnm/par_R * <phi|beta_n><beta_m|phi> ) ******/
                            if (fabs (kbpsi[iip1 + ip1]) > 0.)
                                Hij[ist] += alpha * dnmI_R[ip2 * nh + ip1] * kbpsi[iip1 + ip1] * kbpsi[iip2 + ip2];

                            /***** from sum_nm[Dnm*(<phi|par_beta_n><beta_m|phi> + <phi|beta_n><par_beta_m|phi>)] ****/
                            if ((fabs (kbpsi[iip1 + ip1]) > 0.) || (fabs(kbpsi[iip2 + ip2]) > 0.))
                            {
                                temp = partial_kbpsi[iip1 + ip1] * kbpsi[iip2 + ip2] 
                                            + kbpsi[iip1 + ip1] * partial_kbpsi[iip2 + ip2];
                                Hij[ist] += alpha * dnmI[ip1 * nh + ip2] * temp;
                                Sij[ist] += alpha * qqq[ip1 * nh + ip2] * temp;
                            }
                        }
                    }
                }
                iip1 += ct.max_nl;
                iip2 += ct.max_nl;

            }/* end for ion */

        }/* end for st1 and st2 */


    /* Now calculate the part that kbpsi is stored in other processors */

    for (idx = 1; idx < NPES; idx++)
    {

        proc1 = pct.thispe + idx;
        if (proc1 >= NPES)
            proc1 = proc1 - NPES;
        proc2 = pct.thispe - idx;
        if (proc2 < 0)
            proc2 += NPES;


        MPI_Sendrecv (kbpsi, size, MPI_DOUBLE, proc1, idx, kbBpsi_comm, size,
                      MPI_DOUBLE, proc2, idx, MPI_COMM_WORLD, &mstatus);

        MPI_Sendrecv(partial_kbpsi, size, MPI_DOUBLE, proc1, idx, partial_kbpsi_comm, size, 
                      MPI_DOUBLE, proc2, idx, MPI_COMM_WORLD, &mstatus);


        for (st1 = state_begin[proc]; st1 < state_end[proc]; st1++)
            for (st2 = state_begin[proc2]; st2 < state_end[proc2]; st2++)
            {

                ist = st1 * ct.num_states + st2;
                iip1 = (st1 - state_begin[proc]) * num_nonlocal_ion[proc] * ct.max_nl;
                iip2 = (st2 - state_begin[proc2]) * num_nonlocal_ion[proc2] * ct.max_nl;

                for (ion1 = 0; ion1 < num_nonlocal_ion[proc]; ion1++)
                    for (ion2 = 0; ion2 < num_nonlocal_ion[proc2]; ion2++)
                    {
                        ion1_global = ionidx_allproc[proc * max_ion_nonlocal + ion1];
                        ion2_global = ionidx_allproc[proc2 * max_ion_nonlocal + ion2];

                        if ((ion1_global == ion2_global) && (iion == ion1_global))
                        {

                            nh = pct.prj_per_ion[ion1_global];

                            qqq = pct.qqq[ion1_global];
                            dnmI = pct.dnmI[ion1_global];
                            switch (flag)
                            {
                                case 1:
                                    dnmI_R = pct.dnmI_x[ion1_global];
                                    break;
                                case 2:
                                    dnmI_R = pct.dnmI_y[ion1_global];
                                    break;
                                case 3:
                                    dnmI_R = pct.dnmI_z[ion1_global];
                                    break;
                                default:
                                    error_handler("Undefined flag");
                            }

                            for (ip1 = 0; ip1 < nh; ip1++)
                            {
                                iip1a = iip1 + ion1 * ct.max_nl + ip1;
                                for (ip2 = 0; ip2 < nh; ip2++)
                                {
                                    iip2a = iip2 + ion2 * ct.max_nl + ip2;

                                    /* from sum_nm ( par_Dnm/par_R * <phi|beta_n><beta_m|phi> ) */
                                    if (fabs (kbpsi[iip1a]) > 0.)
                                        Hij[ist] += alpha * dnmI_R[ip2 * nh + ip1] * kbpsi[iip1a] * kbBpsi_comm[iip2a];

                                    /* from sum_nm[Dnm*(<phi|par_beta_n><beta_m|phi> + <phi|beta_n><par_beta_m|phi>)] */
                                    if ((fabs (kbpsi[iip1a]) > 0.) || (fabs(partial_kbpsi[iip1a]) > 0.))
                                    {
                                        temp = partial_kbpsi[iip1a] * kbBpsi_comm[iip2a] 
                                                    + kbpsi[iip1a] * partial_kbpsi_comm[iip2a];
                                        Hij[ist] += alpha * dnmI[ip1 * nh + ip2] * temp;
                                        Sij[ist] += alpha * qqq[ip1 * nh + ip2] * temp;
                                    }

                                }
                            }
                        }
                    }
            }
    }

    /* symmetrize the Hij_R */
    for (st1 = 0; st1 < ct.num_states - 1; st1++)
        for (st2 = st1 + 1; st2 < ct.num_states; st2++)
        {
            idx1 = st1 + st2 * ct.num_states;
            idx2 = st2 + st1 * ct.num_states;
            Hij[idx1] = 0.5 * (Hij[idx1] + Hij[idx2]);
            Hij[idx2] = Hij[idx1];
            Sij[idx1] = 0.5 * (Sij[idx1] + Sij[idx2]);
            Sij[idx2] = Sij[idx1];
        }

    my_free(partial_kbpsi_comm);

    time2 = my_crtc ();
    md_timings (PAR_D_H_AND_S, time2 - time1);


}
