/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include "main.h"

void nlforce_par_gamma (REAL * par_gamma, int ion, int nh, REAL *force)
{
    int idx, idx1, size, n, m, three = 3;
    REAL forces[3];
    REAL *gamma_x, *gamma_y, *gamma_z, *dnmI;

    size = nh * (nh + 1) / 2;

    gamma_x = par_gamma;
    gamma_y = gamma_x + size;
    gamma_z = gamma_y + size;

    dnmI = pct.dnmI[ion];

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

    force[0] += forces[0];
    force[1] += forces[1];
    force[2] += forces[2];

}
