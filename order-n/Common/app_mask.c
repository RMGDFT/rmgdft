/************************** SVN Revision Information **************************
 **    $Id: app_mask.c 587 2006-08-18 22:57:56Z miro $    **
******************************************************************************/
 
#include <stdio.h>
#include <assert.h>
#include "md.h"


void app_mask(int istate, double *u, int level)
{
    char *maskptr;
    double time1, time2;
    int idx, istart;


    time1 = my_crtc();

    maskptr = ct.states[istate].lmask[level];

    if (maskptr != NULL)
        for (idx = 0; idx < ct.states[istate].size; idx++)
        {
            u[idx] *= (double) maskptr[idx];
        }
    else
        for (idx = 0; idx < ct.states[istate].size; idx++)
            u[idx] = 0.;

    time2 = my_crtc();
    rmg_timings(MASK_TIME, (time2 - time1), 0);

}
