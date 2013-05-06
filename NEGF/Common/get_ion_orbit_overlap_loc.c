/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
/*

each orbit or the non-local projector is a cubic grid 

for orbit, its corner grid is (ixmin, iymin, izmin) , 
                                 (ixmax, iymax, izmax)
for local projector on ion2
                       (ixstart_loc, iystart_loc, izstart_loc) 
                       (ixend_loc,   iyend_loc,   izend_loc)

Here we determine their overlap regions.

Here is the index for state st1 and ion ion2: 

index =  (st1 - ct.state_begin) * ct.num_ions + ion2; 

*/

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "main.h"


void get_ion_orbit_overlap_loc (STATE * states)
{

    int i1, i2, i3, i4, i5, i6;
    int st1;
    int ion2;
    int flagx, flagy, flagz;
    int num_states_gridpe, index;

    num_states_gridpe = ct.state_end - ct.state_begin;
    my_malloc( ion_orbit_overlap_region_loc, ct.num_ions * num_states_gridpe, ION_ORBIT_OVERLAP );

    for (st1 = ct.state_begin; st1 < ct.state_end; st1++)
        for (ion2 = 0; ion2 < ct.num_ions; ion2++)
        {
            index = (st1 - ct.state_begin) * ct.num_ions + ion2;
            flagx = 0;
            flagy = 0;
            flagz = 0;

            i1 = states[st1].ixmin;
            i2 = states[st1].ixmax;
            i3 = ct.ions[ion2].ixstart_loc;
            i4 = ct.ions[ion2].ixend_loc;
            i5 = max (i1, i3);
            i6 = min (i2, i4);
            if (i6 >= i5)
                flagx = 1;
            ion_orbit_overlap_region_loc[index].xlow1 = i5;
            ion_orbit_overlap_region_loc[index].xhigh1 = i6;

            i3 = ct.ions[ion2].ixstart_loc + NX_GRID;
            i4 = ct.ions[ion2].ixend_loc + NX_GRID;
            i5 = max (i1, i3);
            i6 = min (i2, i4);
            if (i6 >= i5)
            {
                flagx = 1;
                ion_orbit_overlap_region_loc[index].xlow2 = i5;
                ion_orbit_overlap_region_loc[index].xhigh2 = i6;
                ion_orbit_overlap_region_loc[index].xshift = NX_GRID;
            }
            else
            {
                i3 = ct.ions[ion2].ixstart_loc - NX_GRID;
                i4 = ct.ions[ion2].ixend_loc - NX_GRID;
                i5 = max (i1, i3);
                i6 = min (i2, i4);
                if (i6 >= i5)
                    flagx = 1;
                ion_orbit_overlap_region_loc[index].xlow2 = i5;
                ion_orbit_overlap_region_loc[index].xhigh2 = i6;
                ion_orbit_overlap_region_loc[index].xshift = -NX_GRID;
            }

            i1 = states[st1].iymin;
            i2 = states[st1].iymax;
            i3 = ct.ions[ion2].iystart_loc;
            i4 = ct.ions[ion2].iyend_loc;
            i5 = max (i1, i3);
            i6 = min (i2, i4);
            if (i6 >= i5)
                flagy = 1;
            ion_orbit_overlap_region_loc[index].ylow1 = i5;
            ion_orbit_overlap_region_loc[index].yhigh1 = i6;

            i3 = ct.ions[ion2].iystart_loc + NY_GRID;
            i4 = ct.ions[ion2].iyend_loc + NY_GRID;

            i5 = max (i1, i3);
            i6 = min (i2, i4);
            if (i6 >= i5)
            {
                flagy = 1;
                ion_orbit_overlap_region_loc[index].ylow2 = i5;
                ion_orbit_overlap_region_loc[index].yhigh2 = i6;
                ion_orbit_overlap_region_loc[index].yshift = NY_GRID;
            }
            else
            {
                i3 = ct.ions[ion2].iystart_loc - NY_GRID;
                i4 = ct.ions[ion2].iyend_loc - NY_GRID;
                i5 = max (i1, i3);
                i6 = min (i2, i4);
                if (i6 >= i5)
                    flagy = 1;
                ion_orbit_overlap_region_loc[index].ylow2 = i5;
                ion_orbit_overlap_region_loc[index].yhigh2 = i6;
                ion_orbit_overlap_region_loc[index].yshift = -NY_GRID;
            }

            i1 = states[st1].izmin;
            i2 = states[st1].izmax;
            i3 = ct.ions[ion2].izstart_loc;
            i4 = ct.ions[ion2].izend_loc;

            i5 = max (i1, i3);
            i6 = min (i2, i4);

            if (i6 >= i5)
                flagz = 1;
            ion_orbit_overlap_region_loc[index].zlow1 = i5;
            ion_orbit_overlap_region_loc[index].zhigh1 = i6;

            i3 = ct.ions[ion2].izstart_loc + NZ_GRID;
            i4 = ct.ions[ion2].izend_loc + NZ_GRID;

            i5 = max (i1, i3);
            i6 = min (i2, i4);

            if (i6 >= i5)
            {
                flagz = 1;
                ion_orbit_overlap_region_loc[index].zlow2 = i5;
                ion_orbit_overlap_region_loc[index].zhigh2 = i6;
                ion_orbit_overlap_region_loc[index].zshift = NZ_GRID;
            }
            else
            {
                i3 = ct.ions[ion2].izstart_loc - NZ_GRID;
                i4 = ct.ions[ion2].izend_loc - NZ_GRID;
                i5 = max (i1, i3);
                i6 = min (i2, i4);
                if (i6 >= i5)
                    flagz = 1;
                ion_orbit_overlap_region_loc[index].zlow2 = i5;
                ion_orbit_overlap_region_loc[index].zhigh2 = i6;
                ion_orbit_overlap_region_loc[index].zshift = -NZ_GRID;
            }

            ion_orbit_overlap_region_loc[index].flag = flagx * flagy * flagz;
        }

}
