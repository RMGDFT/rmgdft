/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
/*
                            make_mask_grid.c

    make mask_grid for each state according their radius and mask level  

*/

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "md.h"


void make_mask_grid_state(int level, STATE * states)
{
    int idx, ix, iy, iz;
    int ii;
    char *ptr_mask, *maskptr;
    int dim[3], incx, incy;
    char myMask;
    int state;
    int cut_grid[5];
   
    cut_grid[0] = 2;
    cut_grid[1] = 2;
    cut_grid[2] = 2;
    cut_grid[3] = 2;
    cut_grid[4] = 2;
    


    assert(level < 5);

    for (state = ct.state_begin; state < ct.state_end; state++)
    {
        dim[0] = (states[state].orbit_nx + (1<<level) -1) / (1 << level);
        dim[1] = (states[state].orbit_ny + (1<<level) -1) / (1 << level);
        dim[2] = (states[state].orbit_nz + (1<<level) -1) / (1 << level);
        incy = dim[2];
        incx = dim[2] * dim[1];

        maskptr = states[state].lmask[level];

        for( idx = 0; idx < dim[0] * dim[1] * dim[2]; idx++)
        {
            maskptr[idx] = 1;
        }

        for (ix = 0; ix < dim[0]; ix++)
            for (iy = 0; iy < dim[1]; iy++)
                for (iz = 0; iz < dim[2]; iz++)
                {
                    if(ix <cut_grid[level] | ix >=dim[0]-cut_grid[level]| iy <cut_grid[level] | iy >=dim[1]-cut_grid[level]| iz <cut_grid[level] | iz >=dim[2]-cut_grid[level]) 
                    maskptr[ix * incx + iy * incy + iz] = 0;
                }
    }

        printf(" make_mask_grid_state is  done\n");
}
