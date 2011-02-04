/************************** SVN Revision Information **************************
 **    $Id: pack_vtot_ftoc.c 1242 2011-02-02 18:55:23Z luw $    **
******************************************************************************/
 
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "md.h"

void pack_vtot_ftoc (REAL * vt, REAL * vt_c)
{
    int i, j, k;
    int ixs, iys, ix1, iy1;

    ixs = FPZ0_GRID * FPY0_GRID;
    iys = FPZ0_GRID;

    ix1 = PZ0_GRID * PY0_GRID;
    iy1 = PZ0_GRID;

    for (i = 0; i < PX0_GRID; i++)
    {
        for (j = 0; j < PY0_GRID; j++)
        {
            for (k = 0; k < PZ0_GRID; k++)
            {
                vt_c[i * ix1 + j * iy1 + k] =
                    vt[(RHO_NX * i) * ixs + (RHO_NY * j) * iys + RHO_NZ * k];
            }
        }
    }
}
