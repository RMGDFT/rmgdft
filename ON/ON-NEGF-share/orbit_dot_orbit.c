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
    int state_per_proc;
    int num_recv;
    int st1;


    state_per_proc = ct.state_per_proc + 2;

#pragma omp parallel private(st1)
{
#pragma omp for schedule(dynamic) nowait
    for (st1 = ct.state_begin; st1 < ct.state_end; st1++)
    {
        int st11 = st1 - ct.state_begin;
        for (int st2 = ct.state_begin; st2 < ct.state_end; st2++)
        {
            if (state_overlap_or_not[st11 * ct.num_states + st2] == 1)
            {
                double H, S;
                dot_product_orbit_orbit(&states1[st1], &states[st2], &states[st1],  &H, &S);
                Aij[st11 * ct.num_states + st2] = H;
                Bij[st11 * ct.num_states + st2] = S;
            }
        }
    }
}

    for (int loop = 0; loop < num_sendrecv_loop1; loop++)
    {

        num_recv = recv_from1[loop * state_per_proc + 1];

        for(int i = 0; i< num_recv; i++)
        {

            int st2 = recv_from1[loop * state_per_proc + i + 2];
#pragma omp parallel private(st1)
{
#pragma omp for schedule(dynamic) nowait
            for (st1 = ct.state_begin; st1 < ct.state_end; st1++)
            {
                int st11 = st1 - ct.state_begin;
                if (state_overlap_or_not[st11 * ct.num_states + st2] == 1)
                {
                    double H, S;
                    dot_product_orbit_orbit(&states1[st1], &states[st2], &states[st1], &H, &S);
                    Aij[st11 * ct.num_states + st2] = H;
                    Bij[st11 * ct.num_states + st2] = S;
                }
            }
}
        }

    }



}
