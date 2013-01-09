/************************** SVN Revision Information **************************
 **    $Id$    **
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

    double hx4, ax, ay, az, rmin, x, y, z;


    /* Get lattice vectors */
    ax = NX_GRID * ct.hxgrid * ct.xside;
    ay = NY_GRID * ct.hygrid * ct.yside;
    az = NZ_GRID * ct.hzgrid * ct.zside;

/*  add hx4 to account derivative */
    hx4 = 4.0* ct.hxgrid * ct.xside;
    

    for (state1 = 0; state1 < ct.num_states; state1++)
    {
        for (state2 = state1; state2 < ct.num_states; state2++)
        {
            x = fabs(states[state1].crds[0] - states[state2].crds[0]);
            y = fabs(states[state1].crds[1] - states[state2].crds[1]);
            z = fabs(states[state1].crds[2] - states[state2].crds[2]);

            r1 = states[state1].radius;
            r2 = states[state2].radius;

            if( x > ax/2.0 ) x = ax -x;
            if( y > ay/2.0 ) y = ay -y;
            if( z > az/2.0 ) z = az -z;

            if (x < (r1+r2+hx4) && y < (r1+r2+hx4) && z < (r1+r2+hx4))
                state_overlap_or_not[state1 * ct.num_states + state2] = 1;
            else
                state_overlap_or_not[state1 * ct.num_states + state2] = 0;	

            state_overlap_or_not[state1 + state2 * ct.num_states]
                = state_overlap_or_not[state1 * ct.num_states + state2];
        }
    }

}
