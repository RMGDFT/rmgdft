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

    /*If the grids are the same, just copy the data */
    if (grid_ratio == 1)
    {
        idx = FPX0_GRID * FPY0_GRID * FPZ0_GRID;
        dcopy (&idx, vtot, &ione, vtot_psi, &ione);
    }
    /* For different grids, restriction algorithm is used to obtain potential on coarse grid */
    else
    {
        mg_restrict_6 (vtot, vtot_psi, FPX0_GRID, FPY0_GRID, FPZ0_GRID, grid_ratio);
    }
}
