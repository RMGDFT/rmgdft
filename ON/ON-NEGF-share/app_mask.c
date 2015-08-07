/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
#include <stdio.h>
#include <assert.h>
#include "main.h"
#include "prototypes_on.h"


void app_mask(int istate, double *u, int level)
{
    char *maskptr;
    int idx;
    int size, nx, ny, nz, item;



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


}
