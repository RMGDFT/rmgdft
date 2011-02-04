/************************** SVN Revision Information **************************
 **    $Id: nlforce_par_Q.c 1242 2011-02-02 18:55:23Z luw $    **
******************************************************************************/
 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "md.h"

void nlforce_par_Q(REAL * rho_nm, int ion, int nh, REAL *forces)
{
    int idx, idx1, idx2, n, m, size, max_nl;
    REAL *dnmI_x, *dnmI_y, *dnmI_z, *gamma;
    ION *iptr;

    REAL time1, time2;
    time1 = my_crtc();

    size = ct.max_nl * ct.max_nl;
    max_nl = ct.max_nl;

    gamma = rho_nm + ion * size;

    dnmI_x = pct.dnmI_x[ion];
    dnmI_y = pct.dnmI_y[ion];
    dnmI_z = pct.dnmI_z[ion];
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
                forces[3 * ion] += dnmI_x[idx] * gamma[idx1];
                forces[3 * ion + 1] += dnmI_y[idx] * gamma[idx1];
                forces[3 * ion + 2] += dnmI_z[idx] * gamma[idx1];
            }
            else
            {
                forces[3 * ion] += dnmI_x[idx] * (gamma[idx1] + gamma[idx2]);
                forces[3 * ion + 1] += dnmI_y[idx] * (gamma[idx1] + gamma[idx2]);
                forces[3 * ion + 2] += dnmI_z[idx] * (gamma[idx1] + gamma[idx2]);
            }

        }
    }

    time2 = my_crtc();
    md_timings(NLFORCE_PAR_Q, time2 - time1);

}
