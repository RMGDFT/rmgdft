/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include "main.h"

void debug_write_rho_z (REAL * rhoz)
{
    int k;
    int basis;
    FILE *ftpr;

    if (pct.thispe == 0)
    {
        my_fopen (ftpr, "rho_soft_separate.txt", "a+");

        basis = (PX0_GRID / 2) * PY0_GRID * PZ0_GRID + (PY0_GRID / 2) * PZ0_GRID;
        for (k = 0; k < PZ0_GRID; k++)
            fprintf (ftpr, "%d   %f \n", k, rhoz[basis + k]);

        fclose (ftpr);
    }
}
