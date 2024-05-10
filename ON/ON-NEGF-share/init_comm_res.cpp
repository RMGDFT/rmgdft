/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
/*

init_comm.c

(1) we build up a matrix pct.grid_npes * pct.grid_npes and see which pairs need communications
    and how many orbits we need send
(2) pick up the smallest number first, build up a sendrecv pairs.
(3) pick up the smallest number again, and build up another sendrecv pairs.
    until there is nothing left.
(4) this is used in calculate the residual part 
	Phi_i(r) = Theta_ij * phi_j(r)
	in mg_eig.c with function of get_nonortho_res(...)



num_sendrecv_loop1: 
send_to1[loop * (ct.num_state_per_proc +2) ]:
	for each loop, send_to1[loop][0]: the processor rank we send to
	               send_to1[loop][1]: number of state we need to send
		       send_to1[loop][2, 3, ..]: state index we send
recv_from1[loop * (ct.num_state_per_proc +2) ]:
	for each loop, recv_from1[loop][0]: the processor rank we recev from
	               recv_from1[loop][1]: number of state we need to recv
		       recv_from1[loop][2, 3, ..]: state index we recv
		       



*/

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "main.h"
#include "prototypes_on.h"
#include "init_var.h"

void init_comm_res(STATE * states)
{
    int proc1, proc2, proc3;
    int communication_pair, state_per_proc, loop;
    int *matrix_pairs, *proc_recv;
    int st1, st2, idx;
    int max_send_states;
    int send_order;
    int recv_order;
    int num_overlap;
    int recv_proc;
    int st11;

    state_per_proc = ct.state_per_proc + 2;
    my_calloc( matrix_pairs, pct.grid_npes * pct.grid_npes, int );
    my_calloc( proc_recv, pct.grid_npes, int );
    recv_proc = MPI_PROC_NULL;

    if (send_to1 != NULL)
        my_free(send_to1);
    my_calloc( send_to1, 2 * pct.grid_npes * state_per_proc, int );

    if (recv_from1 != NULL)
        my_free(recv_from1);
    my_calloc( recv_from1, 2 * pct.grid_npes * state_per_proc, int );

    for (loop = 0; loop < 2 * pct.grid_npes; loop++)
    {
        send_to1[loop * state_per_proc] = MPI_PROC_NULL;
        send_to1[loop * state_per_proc + 1] = 0;
        recv_from1[loop * state_per_proc] = MPI_PROC_NULL;
        recv_from1[loop * state_per_proc + 1] = 0;
    }

    for (proc1 = 0; proc1 < pct.grid_npes * pct.grid_npes; proc1++)
        matrix_pairs[proc1] = 0;

proc1 = pct.gridpe;
        for (proc2 = 0; proc2 < pct.grid_npes; proc2++)
        {
            num_overlap = 0;
            for (st1 = ct.state_begin; st1 < ct.state_end; st1++)
            {
                st11 = st1 - ct.state_begin;
                for (st2 = 0; st2 < ct.num_states; st2++)
                    if ( (states[st2].pe == proc2)
                            && (proc2 != proc1) && state_overlap_or_not[st11 * ct.num_states + st2] == 1)
                    {
                        num_overlap++;
                        break;
                    }
            }
            matrix_pairs[proc1 * pct.grid_npes + proc2] = num_overlap;
        }

        idx = pct.grid_npes * pct.grid_npes;
        global_sums_int(matrix_pairs, &idx);

        /*
         *if (pct.gridpe == 0)
         *{
         *    rmg_printf("\n initial communication matrix ");
         *    for (i = 0; i < pct.grid_npes; i++)
         *    {
         *        rmg_printf("\n");
         *        for (j = 0; j < pct.grid_npes; j++)
         *            rmg_printf(" %d ", matrix_pairs[i * pct.grid_npes + j]);
         *    }
         *}
         */

        for (loop = 0; loop < 2 * pct.grid_npes; loop++)
        {
            for (proc1 = 0; proc1 < pct.grid_npes; proc1++)
                proc_recv[proc1] = 0;
            communication_pair = 0;
            for (proc1 = 0; proc1 < pct.grid_npes; proc1++)
            {
                max_send_states = 0;
                for (proc3 = proc1 + 1; proc3 < proc1 + pct.grid_npes; proc3++)
                {
                    proc2 = proc3;
                    if (proc2 > pct.grid_npes - 1)
                        proc2 = proc2 - pct.grid_npes;
                    if (matrix_pairs[proc1 * pct.grid_npes + proc2] > max_send_states && proc_recv[proc2] == 0)
                    {
                        max_send_states = matrix_pairs[proc1 * pct.grid_npes + proc2];
                        recv_proc = proc2;
                    }
                }

                if (max_send_states > 0)
                {
                    proc_recv[recv_proc] = 1;
                    if (pct.gridpe == proc1)
                    {
                        send_to1[loop * state_per_proc] = recv_proc;
                        send_to1[loop * state_per_proc + 1] = max_send_states;
                    }
                    if (pct.gridpe == recv_proc)
                    {
                        recv_from1[loop * state_per_proc] = proc1;
                        recv_from1[loop * state_per_proc + 1] = max_send_states;
                    }
                    matrix_pairs[proc1 * pct.grid_npes + recv_proc] = 0;
                    communication_pair++;
                }
            }
            if (communication_pair == 0)
            {
                num_sendrecv_loop1 = loop;
                loop += 2 * pct.grid_npes;   /* break; */
            }
        }


//printf("\n  num-sendrecv_loop1 %d", num_sendrecv_loop1);
        for (loop = 0; loop < num_sendrecv_loop1; loop++)
        {
            proc1 = pct.gridpe;
            proc2 = send_to1[loop * state_per_proc];
            send_order = 0;
            if (proc2 != MPI_PROC_NULL)
                for (st1 = ct.state_begin; st1 < ct.state_end; st1++)
                {
                    st11 = st1 - ct.state_begin;
                    for (st2 = 0; st2 < ct.num_states; st2++)
                        if (states[st2].pe == proc2
                                && state_overlap_or_not[st11 * ct.num_states + st2] == 1)
                        {
                            send_to1[loop * state_per_proc + send_order + 2] = st1;
                            send_order++;
                            break;
                        }
                }
            assert(send_order == send_to1[loop * state_per_proc + 1]);

            proc1 = pct.gridpe;
            proc2 = recv_from1[loop * state_per_proc];
            recv_order = 0;
            if (proc2 != MPI_PROC_NULL)
                for (st2 = 0; st2 < ct.num_states; st2++)
                {
                    if (states[st2].pe == proc2)
                        for (st1 = ct.state_begin; st1 < ct.state_end; st1++)
                        {
                            st11 = st1 - ct.state_begin;
                            if (state_overlap_or_not[st11 * ct.num_states + st2] == 1)
                            {
                                recv_from1[loop * state_per_proc + recv_order + 2] = st2;
                                recv_order++;
                                break;
                            }
                        }
                }
            assert(recv_order == recv_from1[loop * state_per_proc + 1]);
        }


#if 	DEBUG
        rmg_printf("\n communication info in init_comm_res\n ");
        for (loop = 0; loop < num_sendrecv_loop1; loop++)
        {

            rmg_printf("\nLoop: %d  PE:%d send %d states to PE:%d ---- ", loop,
                    pct.gridpe, send_to1[loop * state_per_proc + 1], send_to1[loop * state_per_proc]);
            for (st2 = 0; st2 < send_to1[loop * state_per_proc + 1]; st2++)
                rmg_printf("  %d ", send_to1[loop * state_per_proc + 2 + st2]);

            rmg_printf("\nLoop: %d  PE:%d receive %d states from PE:%d ---- ", loop,
                    pct.gridpe, recv_from1[loop * state_per_proc + 1],
                    recv_from1[loop * state_per_proc]);
            for (st2 = 0; st2 < recv_from1[loop * state_per_proc + 1]; st2++)
                rmg_printf("  %d ", recv_from1[loop * state_per_proc + 2 + st2]);

        }
#endif

        my_free(matrix_pairs);
        my_free(proc_recv);

}
