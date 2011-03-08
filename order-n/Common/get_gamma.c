/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "md.h"


double get_gamma(double *vtot, double small_eig)
{
    int idx, nits;
    double vmax, diag = -1. / ct.Ac, gamma;

    /* Definition of parameter gamma */
    vmax = -100000.0;
    for (idx = 0; idx < P0_BASIS; idx++)
        if (vmax < vtot[idx])
            vmax = vtot[idx];
    vmax = real_max_all(vmax);

    if (pct.thispe == 0)
        printf("\n sssss %f %f %f ", vmax, small_eig, diag);
    nits = ct.eig_parm.gl_pre + ct.eig_parm.gl_pst;

    /* diag * 4^{N_level+1} */
    gamma = diag;
    if (nits > 0)
        for (idx = 0; idx <= ct.eig_parm.levels; idx++)
            gamma *= 4.;

    /* gamma = inverse of the largest eigenvalue for the low frequency error */
    gamma = 1.0 / (2.0 / gamma + vmax + fabs(small_eig));

    if (pct.thispe == 0)
        printf("\n get_gamma %22.16f ", gamma);
    return gamma;
}
