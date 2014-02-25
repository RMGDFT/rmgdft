/************************** SVN Revision Information **************************
 **    $Id$    **
 ******************************************************************************/

#include <float.h>
#include <stdio.h>
#include <assert.h>
#include "main.h"
#include "init_var.h"
#include "LCR.h"



#include "my_scalapack.h"



void nlforce_partial_H_part2 (STATE * states, STATE * states1, rmg_double_t *GHG, rmg_double_t *force)
{
    int i, ion, ion1, st1, st2, idx, idx1, idx2;
    rmg_double_t temp;
    rmg_double_t *psi, *psi1, *psi2, *old_psi, *vloc_psi;
    rmg_double_t *vloc_x, *vloc_y, *vloc_z;
    MPI_Status mstatus;
    int loop, proc1, proc2, size1, size2, state_per_proc;
    int num_sendrecv, num_send, num_recv;
    double time1, time2, time3, time4;
    int st11;

    state_per_proc = ct.state_per_proc + 2;
    time1 = my_crtc ();

    idx1 = lcr[1].num_states / 2;
    idx2 = ct.num_states - lcr[2].num_states / 2;
    for (st1 = ct.state_begin; st1 < ct.state_end; st1++)
    {

        st11 = st1 - ct.state_begin;
        for( st2 = st1; st2 < ct.state_end; st2++)
        {
            if (state_overlap_or_not[st1 * ct.num_states + st2] == 0) break;

            if (((st1 < idx1) && (st2 < idx1)) || 
                    ((st2 >= idx2) && (st1 < idx1)) 
                    || ((st1 >= idx2) && (st2 >= idx2))) break;

            vloc_x = vnuc_x;
            vloc_y = vnuc_y;
            vloc_z = vnuc_z;
            psi = states[st2].psiR;
            vloc_psi = states1[st2].psiR;
            for (ion = 0; ion < pct.n_ion_center_loc; ion++)
            {
                ion1 = pct.ionidx_loc[ion];
                if (ct.ions[ion1].movable)
                {
                    if(vloc_state_overlap_or_not[ion1 * ct.num_states + st2] == 1) 
                    {
                        for(idx = 0; idx < states1[st2].size; idx++) vloc_psi[idx] = 0.0;

                        /**************************** X direction ************************/
                        /* Generate partial_Vloc * psi and store it  in states1[st2].psiR*/
                        partial_vlocpsi (states[st2], ion1, psi, vloc_x, vloc_psi);

                        temp = dot_product_orbit_orbit (&states[st1], &states1[st2]);
                        if (st2 == st1) force[3*ion1] -= 1.0 * temp * GHG[st11 * ct.num_states + st2];
                        else force[3*ion1] -= 2.0 * temp * GHG[st11 * ct.num_states + st2];

                        /**************************** Y direction ************************/
                        partial_vlocpsi (states[st2], ion1, psi, vloc_y, vloc_psi);

                        temp = dot_product_orbit_orbit (&states[st1], &states1[st2]);
                        if (st2 == st1) force[3*ion1 + 1] -= 1.0 * temp * GHG[st11 * ct.num_states + st2];
                        else force[3*ion1 + 1] -= 2.0 * temp * GHG[st11 * ct.num_states + st2];

                        /**************************** Z direction ************************/
                        partial_vlocpsi (states[st2], ion1, psi, vloc_z, vloc_psi);

                        temp = dot_product_orbit_orbit (&states[st1], &states1[st2]);
                        if (st2 == st1) force[3*ion1 + 2] -= 1.0 * temp * GHG[st11 * ct.num_states + st2];
                        else force[3*ion1 + 2] -= 2.0 * temp * GHG[st11 * ct.num_states + st2];
                    }
                }
                vloc_x += ct.max_lpoints;
                vloc_y += ct.max_lpoints;
                vloc_z += ct.max_lpoints;
            }
        }/* end for st2 */
    }/* end for st1 */ 

    time2 = my_crtc();
    rmg_timings(PAR_D_VNUC_LOC, (time2-time1));

    time3 = my_crtc();
    psi2 = orbit_tem;
    for (loop = 0; loop < num_sendrecv_loop; loop++)
    {

        proc1 = send_to[loop * state_per_proc];
        proc2 = recv_from[loop * state_per_proc];
        num_send = send_to[loop * state_per_proc + 1];
        num_recv = recv_from[loop * state_per_proc + 1];
        num_sendrecv = min (num_send, num_recv);

        for (i = 0; i < num_sendrecv; i++)
        {
            st1 = send_to[loop * state_per_proc + i + 2];
            st2 = recv_from[loop * state_per_proc + i + 2];
            size1 = states[st1].size;
            size2 = states[st2].size;
            psi1 = states[st1].psiR;

            MPI_Sendrecv (psi1, size1, MPI_DOUBLE, proc1, i, psi2, size2,
                    MPI_DOUBLE, proc2, i, MPI_COMM_WORLD, &mstatus);

            old_psi = states[st2].psiR;
            states[st2].psiR = psi2;
            for (st1 = ct.state_begin; st1 < ct.state_end; st1++)
            {
                st11 = st1 - ct.state_begin;
                if (state_overlap_or_not[st1 * ct.num_states + st2] == 0) break;

                if (((st1 < idx1) && (st2 < idx1)) || 
                        ((st2 >= idx2) && (st1 < idx1)) || 
                        ((st1 >= idx2) && (st2 < idx1)) || 
                        ((st1 >= idx2) && (st2 >= idx2))) break;

                vloc_x = vnuc_x;
                vloc_y = vnuc_y;
                vloc_z = vnuc_z;
                psi = states[st1].psiR;
                vloc_psi = states1[st1].psiR;

                for (ion = 0; ion < pct.n_ion_center_loc; ion++)
                {
                    ion1 = pct.ionidx_loc[ion];
                    if (ct.ions[ion1].movable)
                    {

                        if(vloc_state_overlap_or_not[ion1 * ct.num_states + st1] == 1) 
                        {
                            for(idx = 0; idx < states1[st1].size; idx++) vloc_psi[idx] = 0.0;
                            /**************************** X direction ************************/
                            /* Generate partial_Vloc * psi and store it  in states1[st1].psiR*/
                            partial_vlocpsi (states[st1], ion1, psi, vloc_x, vloc_psi);
                            temp = dot_product_orbit_orbit (&states1[st1], &states[st2]);
                            force[3*ion1] -= 2.0 * temp * GHG[st11 * ct.num_states + st2];

                            /**************************** Y direction ************************/
                            partial_vlocpsi (states[st1], ion1, psi, vloc_y, vloc_psi);
                            temp = dot_product_orbit_orbit (&states1[st1], &states[st2]);
                            force[3*ion1 + 1] -= 2.0 * temp * GHG[st11 * ct.num_states + st2];

                            /**************************** Z direction ************************/
                            partial_vlocpsi (states[st1], ion1, psi, vloc_z, vloc_psi);
                            temp = dot_product_orbit_orbit (&states1[st1], &states[st2]);
                            force[3*ion1 + 2] -= 2.0 * temp * GHG[st11 * ct.num_states + st2];


                        }

                    } /* end for if */
                    vloc_x += ct.max_lpoints;
                    vloc_y += ct.max_lpoints;
                    vloc_z += ct.max_lpoints;
                } /* end for ion */
            } /* end for st1 */

            states[st2].psiR = old_psi;
        } /* end for i */

        if (num_send < num_recv)
            for (i = num_send; i < num_recv; i++)
            {
                st2 = recv_from[loop * state_per_proc + i + 2];
                size2 = states[st2].size;

                MPI_Recv (psi2, size2, MPI_DOUBLE, proc2, i, MPI_COMM_WORLD, &mstatus);

                old_psi = states[st2].psiR;
                states[st2].psiR = psi2;
                for (st1 = ct.state_begin; st1 < ct.state_end; st1++)
                {
                    st11 = st1 - ct.state_begin;
                    if (state_overlap_or_not[st1 * ct.num_states + st2] == 0) break;

                    if (((st1 < idx1) && (st2 < idx1)) || 
                            ((st2 >= idx2) && (st1 < idx1)) || 
                            ((st1 >= idx2) && (st2 < idx1)) || 
                            ((st1 >= idx2) && (st2 >= idx2))) break;

                    vloc_x = vnuc_x;
                    vloc_y = vnuc_y;
                    vloc_z = vnuc_z;
                    psi = states[st1].psiR;
                    vloc_psi = states1[st1].psiR;

                    for (ion = 0; ion < pct.n_ion_center_loc; ion++)
                    {
                        ion1 = pct.ionidx_loc[ion];
                        if (ct.ions[ion1].movable)
                        {

                            if(vloc_state_overlap_or_not[ion1 * ct.num_states + st1] == 1) 
                            {
                                for(idx = 0; idx < states1[st1].size; idx++) vloc_psi[idx] = 0.0;
                                /**************************** X direction ************************/
                                /* Generate partial_Vloc * psi and store it  in states1[st1].psiR*/
                                partial_vlocpsi (states[st1], ion1, psi, vloc_x, vloc_psi);
                                temp = dot_product_orbit_orbit (&states1[st1], &states[st2]);
                                force[3*ion1] -= 2.0 * temp * GHG[st11 * ct.num_states + st2];

                                /**************************** Y direction ************************/
                                partial_vlocpsi (states[st1], ion1, psi, vloc_y, vloc_psi);
                                temp = dot_product_orbit_orbit (&states1[st1], &states[st2]);
                                force[3*ion1 + 1] -= 2.0 * temp * GHG[st11 * ct.num_states + st2];

                                /**************************** Z direction ************************/
                                partial_vlocpsi (states[st1], ion1, psi, vloc_z, vloc_psi);
                                temp = dot_product_orbit_orbit (&states1[st1], &states[st2]);
                                force[3*ion1 + 2] -= 2.0 * temp * GHG[st11 * ct.num_states + st2];

                            }

                        } /* end for if */
                        vloc_x += ct.max_lpoints;
                        vloc_y += ct.max_lpoints;
                        vloc_z += ct.max_lpoints;
                    } /* enf for ion */
                } /* end for st1 */

                states[st2].psiR = old_psi;
            }

        if (num_send > num_recv)
            for (i = num_recv; i < num_send; i++)
            {
                st1 = send_to[loop * state_per_proc + i + 2];
                size1 = states[st1].size;
                psi1 = states[st1].psiR;
                MPI_Send (psi1, size1, MPI_DOUBLE, proc1, i, MPI_COMM_WORLD);
            }

        my_barrier ();
    }

    time4 = my_crtc();
    rmg_timings(PAR_D_VNUC_COMM, (time4-time3));

}
