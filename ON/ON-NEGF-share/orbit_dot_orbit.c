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
#include "main.h"
#include "prototypes_on.h"
#include "init_var.h"


void orbit_dot_orbit(STATE * states, STATE * states1, double *Aij, double *Bij)
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
    double H, S;


    state_per_proc = ct.state_per_proc + 2;


    for (st1 = ct.state_begin; st1 < ct.state_end; st1++)
    {
        st11 = st1 - ct.state_begin;
        for (st2 = ct.state_begin; st2 < ct.state_end; st2++)
        {
            if (state_overlap_or_not[st11 * ct.num_states + st2] == 1)
            {
                dot_product_orbit_orbit(&states1[st1], &states[st2], &states[st1],  &H, &S);
                Aij[st11 * ct.num_states + st2] = H;
                Bij[st11 * ct.num_states + st2] = S;
            }
        }
    }


    for (loop = 0; loop < num_sendrecv_loop1; loop++)
    {

        num_recv = recv_from1[loop * state_per_proc + 1];

        for(i = 0; i< num_recv; i++)
        {

            st2 = recv_from1[loop * state_per_proc + i + 2];
            size2= states[st2].size;

            for (st1 = ct.state_begin; st1 < ct.state_end; st1++)
            {
                st11 = st1 - ct.state_begin;
                if (state_overlap_or_not[st11 * ct.num_states + st2] == 1)
                {
                    dot_product_orbit_orbit(&states1[st1], &states[st2], &states[st1], &H, &S);
                    Aij[st11 * ct.num_states + st2] = H;
                    Bij[st11 * ct.num_states + st2] = S;
                }
            }

        }

    }



}
