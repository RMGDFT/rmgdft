/************************** SVN Revision Information **************************
 **    $Id$    **
 ******************************************************************************/

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "main.h"
#include "prototypes_on.h"
#include "transition.h"


void state_corner_xyz(STATE * states)
{
    int state;
    double hgrid[3];
    int ixx, iyy, izz;

    hgrid[0] = get_hxgrid() * get_xside();
    hgrid[1] = get_hygrid() * get_yside();
    hgrid[2] = get_hzgrid() * get_zside();



    /* Loop over states */
    for (state = 0; state < ct.num_states; state++)
    {

        /* find the minmum and maximum grid point of the orbit on ion */
        /* x direction */

        ixx = (int)(states[state].crds[0]/hgrid[0]+0.5);

        states[state].ixmin = ixx - states[state].orbit_nx/2;
        states[state].ixmax = ixx + states[state].orbit_nx/2;

        iyy = (int)(states[state].crds[1]/hgrid[1]+0.5);
        states[state].iymin = iyy - states[state].orbit_ny/2;
        states[state].iymax = iyy + states[state].orbit_ny/2;

        izz = (int)(states[state].crds[2]/hgrid[2]+0.5);
        states[state].izmin = izz - states[state].orbit_nz/2;
        states[state].izmax = izz + states[state].orbit_nz/2;
    }


    if(pct.gridpe == 0)
    {
        rmg_printf("\n lagest orbital size: %d %d %d", ct.max_orbit_nx, ct.max_orbit_ny, ct.max_orbit_nz);
    }

}
