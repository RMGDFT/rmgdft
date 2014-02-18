/************************** SVN Revision Information **************************
 **    $Id$    **
 ******************************************************************************/

/*


   orbit_dot_orbit.c

   work_matrix_row(i,j) = <states[i].psiR in this pe | states1[j].psiR>


 */

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "main_on.h"



void orbit_dot_orbit(STATE * states, STATE * states1, REAL *Aij, REAL *Bij)
{
    int i, ii;
    int st1, st2;
    REAL temp;
    REAL *psi1;
    REAL *psi2;
    REAL *psi3;
    MPI_Status mstatus;
    REAL time1;
    int loop, proc1, proc2, size1, size2, state_per_proc;
    int num_send, num_recv;
    REAL sum;
    int idx;
    REAL temp2, temp3, temp4;
    MPI_Request mr_send, *mr_recv;
    int st11;
    double H, S;
	double time0;


    my_barrier();
    state_per_proc = ct.state_per_proc + 2;
    time1 = my_crtc();


    for (st1 = ct.state_begin; st1 < ct.state_end; st1++)
        for (st2 = ct.state_begin; st2 < ct.state_end; st2++)
            if (state_overlap_or_not[st1 * ct.num_states + st2] == 1)
            {
                st11 = st1 - ct.state_begin;
                dot_product_orbit_orbit(&states1[st1], &states[st2], &states[st1],  &H, &S);
                Aij[st11 * ct.num_states + st2] = H;
                Bij[st11 * ct.num_states + st2] = S;

            }

    int max_ii = 0;
    for (loop = 0; loop < num_sendrecv_loop1; loop++)
    {
        num_send = send_to1[loop * state_per_proc + 1];
        num_recv = recv_from1[loop * state_per_proc + 1];
        max_ii = max(max_ii, num_send);
        max_ii = max(max_ii, num_recv);
    }

    max_ii = int_max_all(max_ii);

    ii = num_sendrecv_loop1 * (max_ii +10) +1;

    my_calloc(mr_recv, ii, MPI_Request);

    psi2 = orbit_tem;
    my_malloc_init(psi3, ct.max_orbit_size, REAL );


    for (loop = 0; loop < num_sendrecv_loop1; loop++)
    {

	my_barrier();
        proc1 = send_to1[loop * state_per_proc];
        proc2 = recv_from1[loop * state_per_proc];
        num_send = send_to1[loop * state_per_proc + 1];
        num_recv = recv_from1[loop * state_per_proc + 1];

        ii = loop * max_ii +1;

        if(num_recv >0)
        {
            i=0;
            ii++;

            st2 = recv_from1[loop * state_per_proc + i + 2];
            size2 = states[st2].size;

            if(ii%2 == 0) MPI_Irecv(psi2, size2, MPI_DOUBLE, proc2, ii, MPI_COMM_WORLD, &mr_recv[ii]);
            if(ii%2 == 1) MPI_Irecv(psi3, size2, MPI_DOUBLE, proc2, ii, MPI_COMM_WORLD, &mr_recv[ii]);
        }


        ii = loop * max_ii +1;

        for (i = 0; i < num_send; i++)
        {
            ii++;
            st1 = send_to1[loop * state_per_proc + i + 2];
            psi1 = states[st1].psiR;
            size1 = states[st1].size;

            MPI_Isend(psi1, size1, MPI_DOUBLE, proc1, ii, MPI_COMM_WORLD, &mr_send);
            MPI_Request_free(&mr_send);
        }   


        ii = loop * max_ii +2;
        
        for(i = 1; i< num_recv+1; i++)
        {
            ii++;
            MPI_Wait(&mr_recv[ii-1], &mstatus);

            st2 = recv_from1[loop * state_per_proc + i-1 + 2];
            size2= states[st2].size;

            if(i != num_recv)
            {
                if(ii%2 == 0) MPI_Irecv(psi2, size2, MPI_DOUBLE, proc2, ii, MPI_COMM_WORLD, &mr_recv[ii]);
                if(ii%2 == 1) MPI_Irecv(psi3, size2, MPI_DOUBLE, proc2, ii, MPI_COMM_WORLD, &mr_recv[ii]);
            }

            if(ii%2 == 1) states[st2].psiR = psi2;
            if(ii%2 == 0) states[st2].psiR = psi3;

            for (st1 = ct.state_begin; st1 < ct.state_end; st1++)
                if (state_overlap_or_not[st1 * ct.num_states + st2] == 1)
                {
                    st11 = st1 - ct.state_begin;
                    dot_product_orbit_orbit(&states1[st1], &states[st2], &states[st1], &H, &S);
                    Aij[st11 * ct.num_states + st2] = H;
                    Bij[st11 * ct.num_states + st2] = S;
                }

        }

    }

    time1 = my_crtc() - time1;
    rmg_timings(ORBIT_DOT_ORBIT, time1);

    my_free(psi3);
    my_free(mr_recv);


}
