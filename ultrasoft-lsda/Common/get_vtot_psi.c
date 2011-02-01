/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "main.h"

void get_vtot_psi (REAL * vtot_psi, REAL * vtot)
{
    int i, j, k;
    P0_GRID *ptr_vtot_psi;
    FP0_GRID *ptr_vtot;

    ptr_vtot_psi = (P0_GRID *) vtot_psi;
    ptr_vtot = (FP0_GRID *) vtot;


    for (i = 0; i < PX0_GRID; i++)
    {
        for (j = 0; j < PY0_GRID; j++)
        {
            for (k = 0; k < PZ0_GRID; k++)
            {
                ptr_vtot_psi->s1.b[i][j][k] = ptr_vtot->s1.b[FG_NX * i][FG_NY * j][FG_NZ * k];
            }
        }
    }
}
