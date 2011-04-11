/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include "main.h"

void nlforce1_par_omega (REAL * par_omega, int ion, ION * iptr, int nh)
{
    int idx, idx1, size, n, m, three = 3;
    REAL forces[3];
    REAL *omega_x, *omega_y, *omega_z, *qqq;

    size = nh * (nh + 1) / 2;

    omega_x = par_omega;
    omega_y = omega_x + size;
    omega_z = omega_y + size;

    qqq = pct.qqq[ion];

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
                forces[0] += qqq[idx1] * omega_x[idx];
                forces[1] += qqq[idx1] * omega_y[idx];
                forces[2] += qqq[idx1] * omega_z[idx];
            }
            else
            {
                forces[0] += 2.0 * qqq[idx1] * omega_x[idx];
                forces[1] += 2.0 * qqq[idx1] * omega_y[idx];
                forces[2] += 2.0 * qqq[idx1] * omega_z[idx];
            }

            ++idx;
        }
    }

    if (pct.spin_flag)
	    global_sums (forces, &three, pct.spin_comm);

    iptr->force[ct.fpt[0]][0] -= forces[0];
    iptr->force[ct.fpt[0]][1] -= forces[1];
    iptr->force[ct.fpt[0]][2] -= forces[2];

}
