/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
#include <stdio.h>
#include <assert.h>
#include "md.h"


void app_mask(int istate, double *u, int level)
{
    char *maskptr;
    double time1, time2;
    int idx, istart;
    int size;


    time1 = my_crtc();

    maskptr = ct.states[istate].lmask[level];
    
    size = ct.states[istate].size/(1<<level)/(1<<level)/(1<<level);

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
