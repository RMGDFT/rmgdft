/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
#include <stdio.h>
#include <assert.h>
#include "main_on.h"


void app_mask(int istate, double *u, int level)
{
    char *maskptr;
    double time1, time2;
    int idx, istart;
    int size, nx, ny, nz, item;


    time1 = my_crtc();

    maskptr = ct.states[istate].lmask[level];
    
    item = 1<<level;
    nx = (ct.states[istate].orbit_nx + item -1)/item;
    ny = (ct.states[istate].orbit_ny + item -1)/item;
    nz = (ct.states[istate].orbit_nz + item -1)/item;
    size = nx * ny * nz;

    if (maskptr != NULL)
        for (idx = 0; idx < size; idx++)
        {
            u[idx] *= (double) maskptr[idx];
        }
    else
        for (idx = 0; idx < size; idx++)
            u[idx] = 0.;

    time2 = my_crtc();
    rmg_timings(MASK_TIME, (time2 - time1));

}
