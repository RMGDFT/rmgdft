/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include "main_on.h"

void nlforce_par_omega(rmg_double_t * par_omega_x, rmg_double_t * par_omega_y, rmg_double_t * par_omega_z, int ion, int nh)
{
    int idx, idx1, idx2, size, n, m, max_nl;
    rmg_double_t forces[3];
    rmg_double_t *omega_x, *omega_y, *omega_z, *qqq;
    ION *iptr;

    rmg_double_t time1, time2;
    time1 = my_crtc();

    size = ct.max_nl * ct.max_nl;
    max_nl = ct.max_nl;

    omega_x = par_omega_x + ion * size;
    omega_y = par_omega_y + ion * size;
    omega_z = par_omega_z + ion * size;

    qqq = pct.qqq[ion];
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
                forces[0] += qqq[idx] * omega_x[idx1];
                forces[1] += qqq[idx] * omega_y[idx1];
                forces[2] += qqq[idx] * omega_z[idx1];
            }
            else
            {
                forces[0] += qqq[idx] * (omega_x[idx1] + omega_x[idx2]);
                forces[1] += qqq[idx] * (omega_y[idx1] + omega_y[idx2]);
                forces[2] += qqq[idx] * (omega_z[idx1] + omega_z[idx2]);
            }

        }
    }

    iptr->force[ct.fpt[0]][0] -= forces[0];
    iptr->force[ct.fpt[0]][1] -= forces[1];
    iptr->force[ct.fpt[0]][2] -= forces[2];

    time2 = my_crtc();
    rmg_timings(NLFORCE_PAR_OMEGA, time2 - time1);

}
