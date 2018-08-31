/************************** SVN Revision Information **************************
 **    $Id: orbit_dot_orbit.c 4341 2018-04-12 01:20:36Z emil $    **
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


void OrbitDotOrbit(STATE * states, STATE * states1, double *Aij, double *Bij)
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
        for (int st22 = 0; st22 < ct.num_orbitals_total; st22++)
        {
            int st2 = ct.orbitals_list[st22];
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

}
