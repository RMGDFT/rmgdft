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

extern std::vector<ORBITAL_PAIR> OrbitalPairs;

void OrbitDotOrbitBlock(int pair_start, int pair_end, double *Aij, double *Bij)
{

    for(int ipair = pair_start; ipair < pair_end; ipair++)
    {
        ORBITAL_PAIR onepair;
        onepair = OrbitalPairs[ipair];
        int st1 = onepair.orbital1;
        int st2 = onepair.orbital2;
        int st11 = st1 - ct.state_begin;

        double H, S;
        DotProductOrbitOrbit(&states1[st1], &states[st2], &states[st1],  &H, &S, &onepair);
        Aij[st11 * ct.num_states + st2] = H;
        Bij[st11 * ct.num_states + st2] = S;
    }
}


