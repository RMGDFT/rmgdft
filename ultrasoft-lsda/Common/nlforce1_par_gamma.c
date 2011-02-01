/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include "main.h"

void nlforce1_par_gamma (REAL * par_gamma, int ion, int nh)
{
    int idx, idx1, size, n, m;
    REAL forces[3];
    REAL *gamma_x, *gamma_y, *gamma_z, *dnmI;
    ION *iptr;

    size = nh * (nh + 1) / 2;

    gamma_x = par_gamma;
    gamma_y = gamma_x + size;
    gamma_z = gamma_y + size;

    dnmI = pct.dnmI[ion];
    iptr = &ct.ions[ion];

    for (idx = 0; idx < 3; idx++)
        forces[idx] = 0.0;

    idx = 0;
    for (n = 0; n < nh; n++)
    {
        for (m = n; m < nh; m++)
        {
            idx1 = n * nh + m;
            if (n == m)
            {
                forces[0] += dnmI[idx1] * gamma_x[idx];
                forces[1] += dnmI[idx1] * gamma_y[idx];
                forces[2] += dnmI[idx1] * gamma_z[idx];
            }
            else
            {
                forces[0] += 2.0 * dnmI[idx1] * gamma_x[idx];
                forces[1] += 2.0 * dnmI[idx1] * gamma_y[idx];
                forces[2] += 2.0 * dnmI[idx1] * gamma_z[idx];
            }

            ++idx;
        }
    }

    iptr->force[ct.fpt[0]][0] += forces[0];
    iptr->force[ct.fpt[0]][1] += forces[1];
    iptr->force[ct.fpt[0]][2] += forces[2];

}
