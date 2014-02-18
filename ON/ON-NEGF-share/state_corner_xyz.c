/************************** SVN Revision Information **************************
 **    $Id$    **
 ******************************************************************************/

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "main_on.h"


void state_corner_xyz(STATE * states)
{
    int ix, iy, iz, state;
    REAL hgrid[3];
    REAL tem;
    int ixx, iyy, izz, imin, imax;
    int item;
    int max_nx, max_ny, max_nz;

    hgrid[0] = ct.hxgrid * ct.xside;
    hgrid[1] = ct.hygrid * ct.yside;
    hgrid[2] = ct.hzgrid * ct.zside;



    /* Loop over states */
    for (state = 0; state < ct.num_states; state++)
    {

        /* find the minmum and maximum grid point of the orbit on ion */
        /* x direction */

        ixx = states[state].crds[0]/hgrid[0];
        states[state].ixmin = ixx - states[state].orbit_nx/2;
        states[state].ixmax = ixx + states[state].orbit_nx/2;

        iyy = states[state].crds[1]/hgrid[1];
        states[state].iymin = iyy - states[state].orbit_ny/2;
        states[state].iymax = iyy + states[state].orbit_ny/2;

        izz = states[state].crds[2]/hgrid[2];
        states[state].izmin = izz - states[state].orbit_nz/2;
        states[state].izmax = izz + states[state].orbit_nz/2;
    }


    if(pct.gridpe == 0)
    {
        printf("\n lagest orbital size: %d %d %d", ct.max_orbit_nx, ct.max_orbit_ny, ct.max_orbit_nz);
    }

}
