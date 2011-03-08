/************************** SVN Revision Information **************************
 **    $Id: get_vtot_psi.c 1166 2010-11-17 20:22:00Z btan $    **
******************************************************************************/

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "main.h"

void get_vtot_psi (REAL * vtot_psi, REAL * vtot, int grid_ratio)
{
    int idx, ione = 1;
    REAL *sg_vtot;

    /*If the grids are the same, just copy the data */
    if (grid_ratio == 1)
    {
        idx = FPX0_GRID * FPY0_GRID * FPZ0_GRID;
        dcopy (&idx, vtot, &ione, vtot_psi, &ione);
    }
    /* For different grids, restriction algorithm is used to obtain potential on coarse grid */
    else
    {
        my_malloc (sg_vtot, (FPX0_GRID + 10) * (FPY0_GRID + 10) * (FPZ0_GRID + 10), REAL);
        trade_imagesx (vtot, sg_vtot, FPX0_GRID, FPY0_GRID, FPZ0_GRID, 5);

        mg_restrict_6 (sg_vtot, vtot_psi, FPX0_GRID, FPY0_GRID, FPZ0_GRID, grid_ratio);

        my_free (sg_vtot);
    }
}
