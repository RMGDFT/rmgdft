/************************** SVN Revision Information **************************
 **    $Id: is_state_overlap.c 587 2006-08-18 22:57:56Z miro $    **
******************************************************************************/
 
/*

is_state_overlap.c

*/

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "md.h"


void is_state_overlap(STATE * states)
{

    int state1, state2;
    REAL r;
    REAL r1, r2;

    double hx4;
/*  add hx4 to account derivative */
    hx4 = 4.0* ct.hxgrid * ct.xside;
    

    for (state1 = 0; state1 < ct.num_states; state1++)
    {
        for (state2 = state1; state2 < ct.num_states; state2++)
        {
            r = minimage1(states[state1].crds, states[state2].crds);
            r1 = states[state1].radius;
            r2 = states[state2].radius;

            if (r < (r1 + r2 + hx4))
                state_overlap_or_not[state1 * ct.num_states + state2] = 1;
            else
                state_overlap_or_not[state1 * ct.num_states + state2] = 0;

            state_overlap_or_not[state1 + state2 * ct.num_states]
                = state_overlap_or_not[state1 * ct.num_states + state2];
        }
    }

}
