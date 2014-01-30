/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "main.h"

void get_vtot_psi (rmg_double_t * vtot_psi, rmg_double_t * vtot, int grid_ratio)
{
    int idx, ione = 1;
    int ix, iy,iz, idx1;

    /*If the grids are the same, just copy the data */
    if (grid_ratio == 1)
    {
        idx = get_FPX0_GRID() * get_FPY0_GRID() * get_FPZ0_GRID();
        QMD_dcopy (idx, vtot, ione, vtot_psi, ione);
    }
    /* For different grids, restriction algorithm is used to obtain potential on coarse grid */
    else
    {
//        mg_restrict_6 (vtot, vtot_psi, get_FPX0_GRID, get_FPY0_GRID, get_FPZ0_GRID, grid_ratio);
        for(ix = 0; ix < get_FPX0_GRID()/2; ix++)
        for(iy = 0; iy < get_FPY0_GRID()/2; iy++)
        for(iz = 0; iz < get_FPZ0_GRID()/2; iz++)
        {
            idx = ix * get_FPY0_GRID()/2 * get_FPZ0_GRID()/2 + iy * get_FPZ0_GRID()/2 + iz;
            idx1 = 2 *ix * get_FPY0_GRID() * get_FPZ0_GRID() + 2 *iy * get_FPZ0_GRID() + 2*iz;
            vtot_psi[idx] = vtot[idx1];
        }

    }
}
