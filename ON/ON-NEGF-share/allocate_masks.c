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


void allocate_masks(STATE * states)
{

    int state;
    char *mask_ptr;
    int item;
    int idx, level, nx, ny, nz;


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
