/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
/*


orbit_dot_orbit.c

work_matrix(i,j) = <states[i].psiR | states1[j].psiR>


*/

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "main.h"
#include "init_var.h"
#include "LCR.h"
/* For Linux and MPICH 
 * 	#include "/usr/lib/mpich/include/mpi.h"
 */

rmg_double_t  dot_product_orbit_orbit (STATE *, STATE *);


void orbit_dot_orbit (STATE * states, STATE * states1, rmg_double_t * work_matrix)
{
    int i;
    int st1, st2;
    rmg_double_t temp;
    rmg_double_t *psi1;
    rmg_double_t *psi2;
    rmg_double_t *old_psi;
    MPI_Status mstatus;
    rmg_double_t time1;
    int loop, proc1, proc2, size1, size2, state_per_proc;
    int num_sendrecv, num_send, num_recv;
    rmg_double_t sum;
    int idx;
    rmg_double_t temp2, temp3, temp4;


    state_per_proc = ct.state_per_proc + 2;
    time1 = my_crtc ();


    for (st1 = ct.state_begin; st1 < ct.state_end; st1++)
        for (st2 = st1; st2 < ct.state_end; st2++)
            if (state_overlap_or_not[st1 * ct.num_states + st2] == 1)
            {
                temp = dot_product_orbit_orbit (&states1[st2], &states[st1]);
                work_matrix[st1 * ct.num_states + st2] = temp;
                work_matrix[st2 * ct.num_states + st1] = temp;

            }

    my_barrier ();
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
                if (state_overlap_or_not[st1 * ct.num_states + st2] == 1)
                {
                    temp = dot_product_orbit_orbit (&states1[st1], &states[st2]);
                    work_matrix[st2 * ct.num_states + st1] = temp;
                    work_matrix[st1 * ct.num_states + st2] = temp;
                }
            states[st2].psiR = old_psi;

        }

        if (num_send < num_recv)
            for (i = num_send; i < num_recv; i++)
            {
                st2 = recv_from[loop * state_per_proc + i + 2];
                size2 = states[st2].size;

                MPI_Recv (psi2, size2, MPI_DOUBLE, proc2, i, MPI_COMM_WORLD, &mstatus);

                old_psi = states[st2].psiR;
                states[st2].psiR = psi2;
                for (st1 = ct.state_begin; st1 < ct.state_end; st1++)
                    if (state_overlap_or_not[st1 * ct.num_states + st2] == 1)
                    {
                        temp = dot_product_orbit_orbit (&states1[st1], &states[st2]);
                        work_matrix[st2 * ct.num_states + st1] = temp;
                        work_matrix[st1 * ct.num_states + st2] = temp;

                    }
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

    time1 = my_crtc () - time1;
    rmg_timings (ORBIT_DOT_ORBIT, time1);


#if     DEBUG
    print_sum (ct.num_states * ct.num_states, work_matrix, "sum of work_matrix in orbit_dot_orbit");
    printf ("\n Leave orbit_dot_orbit \n ");
    fflush (NULL);
#endif

}
