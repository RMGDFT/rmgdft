/************************** SVN Revision Information **************************
 **    $Id: duplicate_states_info.c 1242 2011-02-02 18:55:23Z luw $    **
******************************************************************************/
 
#include 	"md.h"


void duplicate_states_info (STATE * states, STATE * states1)
{
    int st;

    for (st = 0; st < ct.num_states; st++)
    {
        states1[st].pe = states[st].pe;
        states1[st].crds[0] = states[st].crds[0];
        states1[st].crds[1] = states[st].crds[1];
        states1[st].crds[2] = states[st].crds[2];
        states1[st].radius = states[st].radius;
        states1[st].movable = states[st].movable;
        states1[st].frozen = states[st].frozen;
        states1[st].index = states[st].index;
        states1[st].ixmin = states[st].ixmin;
        states1[st].iymin = states[st].iymin;
        states1[st].izmin = states[st].izmin;
        states1[st].ixmax = states[st].ixmax;
        states1[st].iymax = states[st].iymax;
        states1[st].izmax = states[st].izmax;
        states1[st].xfold = states[st].xfold;
        states1[st].yfold = states[st].yfold;
        states1[st].zfold = states[st].zfold;
        states1[st].ixstart = states[st].ixstart;
        states1[st].iystart = states[st].iystart;
        states1[st].izstart = states[st].izstart;
        states1[st].ixend = states[st].ixend;
        states1[st].iyend = states[st].iyend;
        states1[st].izend = states[st].izend;
        states1[st].orbit_nx = states[st].orbit_nx;
        states1[st].orbit_nx = states[st].orbit_nx;
        states1[st].orbit_nz = states[st].orbit_nz;
        states1[st].size = states[st].size;
    }

}
