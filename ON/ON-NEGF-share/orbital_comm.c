/************************** SVN Revision Information **************************
 **    $Id: orbit_dot_orbit.c 3140 2015-08-06 15:48:24Z luw $    **
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
#include "main.h"
#include "prototypes_on.h"
#include "init_var.h"


void orbital_comm(STATE * states)
{
    int i, ii;
    int st1, st2;
    double *psi1;
    double *psi2;
    double *psi3;
    MPI_Status mstatus;
    int loop, proc1, proc2, size1, size2, state_per_proc;
    int num_send, num_recv;
    MPI_Request mr_send, *mr_recv;
    int st11, ione = 1;


    my_barrier();
    state_per_proc = ct.state_per_proc + 2;



    int max_ii = 0;
    for (loop = 0; loop < num_sendrecv_loop1; loop++)
    {
        num_send = send_to1[loop * state_per_proc + 1];
        num_recv = recv_from1[loop * state_per_proc + 1];
        max_ii = rmg_max(max_ii, num_send);
        max_ii = rmg_max(max_ii, num_recv);
    }

    max_ii = int_max_all(max_ii);

    ii = num_sendrecv_loop1 * (max_ii +10) +1;

    my_calloc(mr_recv, ii, MPI_Request);

    psi2 = orbit_tem;
    my_malloc_init(psi3, ct.max_orbit_size, double );


    my_barrier();
    for (loop = 0; loop < num_sendrecv_loop1; loop++)
    {

        proc1 = send_to1[loop * state_per_proc];
        proc2 = recv_from1[loop * state_per_proc];
        num_send = send_to1[loop * state_per_proc + 1];
        num_recv = recv_from1[loop * state_per_proc + 1];

        ii = loop * max_ii +1;

        for (i = 0; i < num_send; i++)
        {
            ii++;
            st1 = send_to1[loop * state_per_proc + i + 2];
            psi1 = states[st1].psiR;
            size1 = states[st1].size;

            MPI_Isend(psi1, size1, MPI_DOUBLE, proc1, ii, pct.grid_comm, &mr_recv[ii]);
        }   


        ii = loop * max_ii +1;

        for(i = 0; i< num_recv; i++)
        {
            ii++;

            st2 = recv_from1[loop * state_per_proc + i + 2];
            size2= states[st2].size;

            psi2 = states[st2].psiR;
            MPI_Recv(psi2, size2, MPI_DOUBLE, proc2, ii, pct.grid_comm, &mstatus);


        }
        ii = loop * max_ii +1;

        for (i = 0; i < num_send; i++)
        {
            ii++;
            MPI_Wait(&mr_recv[ii], &mstatus);
        }

    }


    my_barrier();
    my_free(mr_recv);


}
