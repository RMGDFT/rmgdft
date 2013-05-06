/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include "main.h"

void nlforce_par_rho(REAL * par_gamma_x, REAL * par_gamma_y, REAL * par_gamma_z, int ion, int nh)
{
    int idx, idx1, idx2, size, n, m, max_nl;
    REAL forces[3];
    REAL *gamma_x, *gamma_y, *gamma_z, *dnmI;
    ION *iptr;

    REAL time1, time2;
    time1 = my_crtc();

    size = ct.max_nl * ct.max_nl;
    max_nl = ct.max_nl;

    gamma_x = par_gamma_x + ion * size;
    gamma_y = par_gamma_y + ion * size;
    gamma_z = par_gamma_z + ion * size;

    dnmI = pct.dnmI[ion];
    iptr = &ct.ions[ion];

    for (idx = 0; idx < 3; idx++)
        forces[idx] = 0.0;

    for (n = 0; n < nh; n++)
    {
        for (m = n; m < nh; m++)
        {
            idx = n * nh + m;
            idx1 = n * max_nl + m;
            idx2 = m * max_nl + n;
            if (n == m)
            {
                forces[0] += dnmI[idx] * gamma_x[idx1];
                forces[1] += dnmI[idx] * gamma_y[idx1];
                forces[2] += dnmI[idx] * gamma_z[idx1];
            }
            else
            {
                forces[0] += dnmI[idx] * (gamma_x[idx1] + gamma_x[idx2]);
                forces[1] += dnmI[idx] * (gamma_y[idx1] + gamma_y[idx2]);
                forces[2] += dnmI[idx] * (gamma_z[idx1] + gamma_z[idx2]);
            }

        }
    }

    iptr->force[ct.fpt[0]][0] += forces[0];
    iptr->force[ct.fpt[0]][1] += forces[1];
    iptr->force[ct.fpt[0]][2] += forces[2];

    time2 = my_crtc();
    rmg_timings(NLFORCE_PAR_RHO, time2 - time1);

}
