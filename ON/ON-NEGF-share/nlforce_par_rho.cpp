/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include "main.h"
#include "prototypes_on.h"

void nlforce_par_rho(double * par_gamma_x, double * par_gamma_y, double * par_gamma_z, int ion, int nh)
{
    int idx, idx1, idx2, size, n, m, max_nl;
    double forces[3];
    double *gamma_x, *gamma_y, *gamma_z, *dnmI;
    ION *iptr;


    size = ct.max_nl * ct.max_nl;
    max_nl = ct.max_nl;

    gamma_x = par_gamma_x + ion * size;
    gamma_y = par_gamma_y + ion * size;
    gamma_z = par_gamma_z + ion * size;

    iptr = &Atoms[ion];
    dnmI = iptr->dnmI;

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


}
