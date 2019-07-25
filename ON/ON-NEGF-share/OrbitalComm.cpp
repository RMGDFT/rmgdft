/************************** SVN Revision Information **************************
 **    $Id: orbit_dot_orbit.c 3140 2015-08-06 15:48:24Z luw $    **
 ******************************************************************************/

/*


   orbit_dot_orbit.c


 */

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "main.h"
#include "prototypes_on.h"
#include "init_var.h"


void OrbitalComm(STATE * states)
{
    int i, ii;
    int st1, st2;
    MPI_Status mstatus;
    int loop, proc1, proc2, size1, size2, state_per_proc;
    int num_send, num_recv;
    MPI_Request *mr_recv;

    int send_size, recv_size, position;

    MPI_Barrier(pct.img_comm);
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

    mr_recv = new MPI_Request[ii]();

    double *psi3 = new double[ct.max_orbit_size]();


    MPI_Barrier(pct.img_comm);
    for (loop = 0; loop < num_sendrecv_loop1; loop++)
    {

        proc1 = send_to1[loop * state_per_proc];
        proc2 = recv_from1[loop * state_per_proc];
        num_send = send_to1[loop * state_per_proc + 1];
        num_recv = recv_from1[loop * state_per_proc + 1];

        ii = loop * max_ii +1;

        send_size = 0; 
        for (i = 0; i < num_send; i++)
        {
            st1 = send_to1[loop * state_per_proc + i + 2];
            size1 = states[st1].size;
            send_size += size1 * sizeof(double);
        }


        if(send_size >0)
        {
            position = 0;
            for (i = 0; i < num_send; i++)
            {
                st1 = send_to1[loop * state_per_proc + i + 2];
                size1 = states[st1].size;
                MPI_Pack(states[st1].psiR, size1, MPI_DOUBLE, states1[ct.state_begin].psiR, 
                        send_size, &position, pct.grid_comm);
            }

            MPI_Isend(states1[ct.state_begin].psiR, send_size, MPI_BYTE, proc1, ii, pct.grid_comm, &mr_recv[ii]);
        }



        recv_size = 0;
        for(i = 0; i< num_recv; i++)
        {
            st2 = recv_from1[loop * state_per_proc + i + 2];
            recv_size += states[st2].size * sizeof(double);
        }

        if(recv_size >0 )
        {
            MPI_Recv(states_tem[ct.state_begin].psiR, recv_size, MPI_BYTE, proc2, ii, pct.grid_comm, &mstatus);
            position = 0;
            for(i = 0; i< num_recv; i++)
            {
                st2 = recv_from1[loop * state_per_proc + i + 2];
                size2 = states[st2].size;
                MPI_Unpack(states_tem[ct.state_begin].psiR, recv_size, &position, states[st2].psiR, size2, 
                        MPI_DOUBLE, pct.grid_comm);
            }
        }




        if(send_size >0) 
            MPI_Wait(&mr_recv[ii], &mstatus);

    }


    MPI_Barrier(pct.img_comm);
    delete [] psi3;
    delete [] mr_recv;


}
