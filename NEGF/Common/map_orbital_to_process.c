/************************** SVN Revision Information **************************
 **    $Id$    **
 ******************************************************************************/

/*

    map the localized orbital to the processor domain
*/



#include <sys/types.h>
#include <unistd.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>
#include "main.h"


void map_orbital_to_process(int st2, STATE *states, STATE *states_distribute, double *psi_whole)
{

    int idx, idx2;
    int ixmin, ixmax, iymin, iymax, izmin, izmax;
    int incx, incy, ix, iy, iz;

    int st,ixx, iyy, izz; 

    int x_off, y_off, z_off;


    x_off = pct.PX_OFFSET;
    y_off = pct.PY_OFFSET;
    z_off = pct.PZ_OFFSET;

    st = states_distribute[st2].istate;

    ixmin = states[st].ixmin;
    ixmax = states[st].ixmax;
    iymin = states[st].iymin;
    iymax = states[st].iymax;
    izmin = states[st].izmin;
    izmax = states[st].izmax;

    for(idx = 0; idx < pct.P0_BASIS; idx++) states_distribute[st2].psiR[idx] = 0.0;

    for(ix = ixmin; ix <= ixmax; ix++)
        for(iy = iymin; iy <= iymax; iy++)
            for(iz = izmin; iz <= izmax; iz++)
            {
                idx = (ix-ixmin) * states[st].orbit_ny * states[st].orbit_nz  
                    +(iy-iymin) * states[st].orbit_nz + iz - izmin;
                ixx = ix;
                iyy = iy;
                izz = iz;

                if(cei.num_probe <= 2) 
                {
                    if(iyy < 0 ) iyy += NY_GRID;
                    if(iyy >= NY_GRID) iyy -= NY_GRID;
                }

                if(izz < 0 ) izz += NZ_GRID;
                if(izz >= NZ_GRID) izz -= NZ_GRID;

                if(   ixx >= x_off && ixx < x_off + pct.PX0_GRID 
                        && iyy >= y_off && iyy < y_off + pct.PY0_GRID 
                        && izz >= z_off && izz < z_off + pct.PZ0_GRID )
                {

                    idx2 = (ixx - x_off) * pct.PY0_GRID * pct.PZ0_GRID
                        + (iyy - y_off)  * pct.PZ0_GRID + izz - z_off;

                    states_distribute[st2].psiR[idx2] = psi_whole[idx];
                }
            }
}

