/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include "md.h"

void nlforce_par_rho(REAL * par_gamma_x, REAL * par_gamma_y, REAL * par_gamma_z, int ion, int nh, REAL *forces)
{
    int idx, idx1, idx2, size, n, m, max_nl;
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

    for (n = 0; n < nh; n++)
    {
        for (m = n; m < nh; m++)
        {
            idx = n * nh + m;
            idx1 = n * max_nl + m;
            idx2 = m * max_nl + n;
            if (n == m)
            {
                forces[3 * ion] += dnmI[idx] * gamma_x[idx1];
                forces[3 * ion + 1] += dnmI[idx] * gamma_y[idx1];
                forces[3 * ion + 2] += dnmI[idx] * gamma_z[idx1];
            }
            else
            {
                forces[3 * ion] += dnmI[idx] * (gamma_x[idx1] + gamma_x[idx2]);
                forces[3 * ion + 1] += dnmI[idx] * (gamma_y[idx1] + gamma_y[idx2]);
                forces[3 * ion + 2] += dnmI[idx] * (gamma_z[idx1] + gamma_z[idx2]);
            }

        }
    }

    time2 = my_crtc();
    rmg_timings(NLFORCE_PAR_RHO, time2 - time1);

}
