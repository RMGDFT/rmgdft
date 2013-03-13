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

    int ixcenter1, iycenter1, izcenter1;
    int ixcenter2, iycenter2, izcenter2;
    int ixhalf1, iyhalf1, izhalf1;
    int ixhalf2, iyhalf2, izhalf2;

    int xoverlap, yoverlap, zoverlap;

    for (state1 = 0; state1 < ct.num_states; state1++)
    {
        ixcenter1 = (states[state1].ixmin + states[state1].ixmax)/2;
        iycenter1 = (states[state1].iymin + states[state1].iymax)/2;
        izcenter1 = (states[state1].izmin + states[state1].izmax)/2;
        ixhalf1 = (states[state1].ixmax - states[state1].ixmin)/2;
        iyhalf1 = (states[state1].iymax - states[state1].iymin)/2;
        izhalf1 = (states[state1].izmax - states[state1].izmin)/2;

        for (state2 = state1; state2 < ct.num_states; state2++)
        {
            ixcenter2 = (states[state2].ixmin + states[state2].ixmax)/2;
            iycenter2 = (states[state2].iymin + states[state2].iymax)/2;
            izcenter2 = (states[state2].izmin + states[state2].izmax)/2;
            ixhalf2 = (states[state2].ixmax - states[state2].ixmin)/2;
            iyhalf2 = (states[state2].iymax - states[state2].iymin)/2;
            izhalf2 = (states[state2].izmax - states[state2].izmin)/2;

            xoverlap =abs(ixcenter1 - ixcenter2);
            yoverlap =abs(iycenter1 - iycenter2);
            zoverlap =abs(izcenter1 - izcenter2);


            if ((xoverlap <= ixhalf1 + ixhalf2 ) && 
                    (yoverlap <= iyhalf1 + iyhalf2 ) && 
                    (zoverlap <= izhalf1 + izhalf2 ))
                state_overlap_or_not[state1 * ct.num_states + state2] = 1;
            else
                state_overlap_or_not[state1 * ct.num_states + state2] = 0;	

            state_overlap_or_not[state1 + state2 * ct.num_states]
                = state_overlap_or_not[state1 * ct.num_states + state2];

        }
    }

}
