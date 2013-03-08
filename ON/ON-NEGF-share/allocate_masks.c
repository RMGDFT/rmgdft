/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "md.h"


void allocate_masks(STATE * states)
{

    int radius_order[MAX_STATES];
    REAL radius[MAX_STATES];
    int num_masks;
    REAL mask_radius[MAX_STATES];
    int mask_size[MAX_STATES];
    int state;
    int temp_i;
    REAL temp_r;
    char *mask_ptr;
    int i, item;
    int idx, level, nx, ny, nz;
    int order;



    for (state = ct.state_begin; state < ct.state_end; state++)
    {
        idx = states[state].size;
        for (level = 1; level <= ct.eig_parm.levels; level++)
        {
            item = 1<<level;
            nx = (states[state].orbit_nx + item -1)/item;
            ny = (states[state].orbit_ny + item -1)/item;
            nz = (states[state].orbit_nz + item -1)/item;
            idx += nx * ny * nz;
        }

        my_malloc( mask_ptr, idx, char );

        idx = 0;
        for (level = 0; level <= ct.eig_parm.levels; level++)
        {
            states[state].lmask[level] = mask_ptr + idx;
            item = 1<<level;
            nx = (states[state].orbit_nx + item -1)/item;
            ny = (states[state].orbit_ny + item -1)/item;
            nz = (states[state].orbit_nz + item -1)/item;
            idx += nx * ny * nz;
        }
    }
}
