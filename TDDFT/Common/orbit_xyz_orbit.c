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
#include "init_var.h"


void orbit_xyz_orbit(STATE * states, double *Xij, double *Yij, double *Zij)
{
    int i, ii;
    int st1, st2;
    double temp;
    double *psi1;
    double *psi2;
    double *psi3;
    MPI_Status mstatus;
    int loop, proc1, proc2, size1, size2, state_per_proc;
    int num_send, num_recv;
    double sum;
    int idx;
    double temp2, temp3, temp4;
    MPI_Request mr_send, *mr_recv;
    int st11;
    double X0, Y0, Z0;


    my_barrier();
    state_per_proc = ct.state_per_proc + 2;


    for (st1 = ct.state_begin; st1 < ct.state_end; st1++)
    {
        st11 = st1 - ct.state_begin;
        for (st2 = ct.state_begin; st2 < ct.state_end; st2++)
            if (state_overlap_or_not[st11 * ct.num_states + st2] == 1)
            {
                dot_product_orbit_xyz_orbit(&states[st1], &states[st2],  &X0, &Y0, &Z0);
                Xij[st11 * ct.num_states + st2] = X0;
                Yij[st11 * ct.num_states + st2] = Y0;
                Zij[st11 * ct.num_states + st2] = Z0;

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
                    dot_product_orbit_xyz_orbit(&states[st1], &states[st2],  &X0, &Y0, &Z0);
                    Xij[st11 * ct.num_states + st2] = X0;
                    Yij[st11 * ct.num_states + st2] = Y0;
                    Zij[st11 * ct.num_states + st2] = Z0;
                }

            }
        }

    }




}


