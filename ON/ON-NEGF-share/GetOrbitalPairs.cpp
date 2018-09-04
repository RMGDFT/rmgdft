/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
/*
for each orbit, its corner grid is (ixmin, iymin, izmin) , 
                                 (ixmax, iymax, izmax)
Here we determine their overlap regions.
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
void GetOrbitalPairs(STATE * states)
{


    int i1, i2, i3, i4, i5, i6;
    int size;
    int index;
    int st1, st2;
    int num_state_gridpe;
    ORBITAL_PAIR onepair;

    int ixcenter1, iycenter1, izcenter1;
    int ixcenter2, iycenter2, izcenter2;
    int ixhalf1, iyhalf1, izhalf1;
    int ixhalf2, iyhalf2, izhalf2;

    int xoverlap, yoverlap, zoverlap;

    for (st1 = ct.state_begin; st1 < ct.state_end; st1++)
    {
        ixcenter1 = (states[st1].ixmin + states[st1].ixmax)/2;
        iycenter1 = (states[st1].iymin + states[st1].iymax)/2;
        izcenter1 = (states[st1].izmin + states[st1].izmax)/2;
        ixhalf1 = (states[st1].ixmax - states[st1].ixmin)/2;
        iyhalf1 = (states[st1].iymax - states[st1].iymin)/2;
        izhalf1 = (states[st1].izmax - states[st1].izmin)/2;

        for (st2 = 0; st2 < ct.num_states; st2++)
        {
            ixcenter2 = (states[st2].ixmin + states[st2].ixmax)/2;
            iycenter2 = (states[st2].iymin + states[st2].iymax)/2;
            izcenter2 = (states[st2].izmin + states[st2].izmax)/2;
            ixhalf2 = (states[st2].ixmax - states[st2].ixmin)/2;
            iyhalf2 = (states[st2].iymax - states[st2].iymin)/2;
            izhalf2 = (states[st2].izmax - states[st2].izmin)/2;

            xoverlap =abs(ixcenter1 - ixcenter2);
            yoverlap =abs(iycenter1 - iycenter2);
            zoverlap =abs(izcenter1 - izcenter2);

            if (xoverlap > get_NX_GRID()/2) xoverlap = get_NX_GRID() - xoverlap;
            if (yoverlap > get_NY_GRID()/2) yoverlap = get_NY_GRID() - yoverlap;
            if (zoverlap > get_NZ_GRID()/2) zoverlap = get_NZ_GRID() - zoverlap;


            if ((xoverlap <= ixhalf1 + ixhalf2 ) && 
                    (yoverlap <= iyhalf1 + iyhalf2 ) && 
                    (zoverlap <= izhalf1 + izhalf2 ))

            {
                onepair.orbital1 = st1;
                onepair.orbital2 = st2;
                i1 = states[st1].ixmin;
                i2 = states[st1].ixmax;
                i3 = states[st2].ixmin;
                i4 = states[st2].ixmax;
                i5 = rmg_max(i1, i3);
                i6 = rmg_min(i2, i4);

                onepair.xlow1 = i5;
                onepair.xhigh1 = i6;

                i3 = states[st2].ixmin + get_NX_GRID();
                i4 = states[st2].ixmax + get_NX_GRID();
                i5 = rmg_max(i1, i3);
                i6 = rmg_min(i2, i4);

                if (i6 >= i5)
                {
                    onepair.xlow2 = i5;
                    onepair.xhigh2 = i6;
                    onepair.xshift = get_NX_GRID();
                }
                else
                {
                    i3 = states[st2].ixmin - get_NX_GRID();
                    i4 = states[st2].ixmax - get_NX_GRID();
                    i5 = rmg_max(i1, i3);
                    i6 = rmg_min(i2, i4);
                    onepair.xlow2 = i5;
                    onepair.xhigh2 = i6;
                    onepair.xshift = -get_NX_GRID();
                }

                i1 = states[st1].iymin;
                i2 = states[st1].iymax;
                i3 = states[st2].iymin;
                i4 = states[st2].iymax;
                i5 = rmg_max(i1, i3);
                i6 = rmg_min(i2, i4);

                onepair.ylow1 = i5;
                onepair.yhigh1 = i6;

                i3 = states[st2].iymin + get_NY_GRID();
                i4 = states[st2].iymax + get_NY_GRID();
                i5 = rmg_max(i1, i3);
                i6 = rmg_min(i2, i4);

                if (i6 >= i5)
                {
                    onepair.ylow2 = i5;
                    onepair.yhigh2 = i6;
                    onepair.yshift = get_NY_GRID();
                }
                else
                {
                    i3 = states[st2].iymin - get_NY_GRID();
                    i4 = states[st2].iymax - get_NY_GRID();
                    i5 = rmg_max(i1, i3);
                    i6 = rmg_min(i2, i4);
                    onepair.ylow2 = i5;
                    onepair.yhigh2 = i6;
                    onepair.yshift = -get_NY_GRID();
                }


                i1 = states[st1].izmin;
                i2 = states[st1].izmax;
                i3 = states[st2].izmin;
                i4 = states[st2].izmax;
                i5 = rmg_max(i1, i3);
                i6 = rmg_min(i2, i4);

                onepair.zlow1 = i5;
                onepair.zhigh1 = i6;

                i3 = states[st2].izmin + get_NZ_GRID();
                i4 = states[st2].izmax + get_NZ_GRID();
                i5 = rmg_max(i1, i3);
                i6 = rmg_min(i2, i4);

                if (i6 >= i5)
                {
                    onepair.zlow2 = i5;
                    onepair.zhigh2 = i6;
                    onepair.zshift = get_NZ_GRID();
                }
                else
                {
                    i3 = states[st2].izmin - get_NZ_GRID();
                    i4 = states[st2].izmax - get_NZ_GRID();
                    i5 = rmg_max(i1, i3);
                    i6 = rmg_min(i2, i4);

                    onepair.zlow2 = i5;
                    onepair.zhigh2 = i6;
                    onepair.zshift = -get_NZ_GRID();
                }
                
                OrbitalPairs.push_back(onepair);
            }
        }
    }

}
