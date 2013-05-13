/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include "main.h"

void debug_write_rho_z (rmg_double_t * rhoz)
{
    int k;
    int basis;
    FILE *ftpr;

    if (pct.gridpe == 0)
    {
        my_fopen (ftpr, "rho_soft_separate.txt", "a+");

        basis = (pct.PX0_GRID / 2) * pct.PY0_GRID * pct.PZ0_GRID + (pct.PY0_GRID / 2) * pct.PZ0_GRID;
        for (k = 0; k < pct.PZ0_GRID; k++)
            fprintf (ftpr, "%d   %f \n", k, rhoz[basis + k]);

        fclose (ftpr);
    }
}
