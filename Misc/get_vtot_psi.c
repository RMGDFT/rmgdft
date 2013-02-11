/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "main.h"

void get_vtot_psi (REAL * vtot_psi, REAL * vtot, int grid_ratio)
{
    int idx, ione = 1;
    int ix, iy,iz, idx1;

    /*If the grids are the same, just copy the data */
    if (grid_ratio == 1)
    {
        idx = pct.FPX0_GRID * pct.FPY0_GRID * pct.FPZ0_GRID;
        QMD_scopy (idx, vtot, ione, vtot_psi, ione);
    }
    /* For different grids, restriction algorithm is used to obtain potential on coarse grid */
    else
    {
//        mg_restrict_6 (vtot, vtot_psi, pct.FPX0_GRID, pct.FPY0_GRID, pct.FPZ0_GRID, grid_ratio);
    for(ix = 0; ix < pct.FPX0_GRID/2; ix++)
    for(iy = 0; iy < pct.FPY0_GRID/2; iy++)
    for(iz = 0; iz < pct.FPZ0_GRID/2; iz++)
    {
        idx = ix * pct.FPY0_GRID/2 * pct.FPZ0_GRID/2 + iy * pct.FPZ0_GRID/2 + iz;
        idx1 = 2 *ix * pct.FPY0_GRID * pct.FPZ0_GRID + 2 *iy * pct.FPZ0_GRID + 2*iz;
        vtot_psi[idx] = vtot[idx1];
    }

    }
}
