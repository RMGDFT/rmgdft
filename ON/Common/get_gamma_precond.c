/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "main.h"
#include "prototypes_on.h"


double get_gamma_precond(double *vtot, double small_eig)
{
    int idx, nits;
    double vmax, gamma;
    double Zfac = 2.0 * ct.max_zvalence;


    double ihx = 1. / (get_hxgrid() * get_hxgrid() * get_xside() * get_xside());
    double ihy = 1. / (get_hygrid() * get_hygrid() * get_yside() * get_yside());
    double ihz = 1. / (get_hzgrid() * get_hzgrid() * get_zside() * get_zside());
    double diag = (-4. / 3.) * (ihx + ihy + ihz);
    diag = -1. / diag;

    /* Definition of parameter gamma */
    vmax = -100000.0;
    for (idx = 0; idx < get_P0_BASIS(); idx++)
        if (vmax < vtot[idx])
            vmax = vtot[idx];
    vmax = real_max_all(vmax);

    if (pct.gridpe == 0)
        printf("\n sssss %f %f %f ", vmax, small_eig, diag);
    nits = ct.eig_parm.gl_pre + ct.eig_parm.gl_pst;

    /* diag * 4^{N_level+1} */
    gamma = diag;
    if (nits > 0)
        for (idx = 0; idx <= ct.eig_parm.levels; idx++)
            gamma *= 4.;

    /* gamma = inverse of the largest eigenvalue for the low frequency error */
    gamma = 1.0 / (2.0 / gamma + vmax + fabs(small_eig));

//double t5 = diag - Zfac;
//t5 = -1.0 / t5;
//gamma = ct.eig_parm.gl_step * t5;

    if (pct.gridpe == 0)
        printf("\n get_gamma %22.16f ", gamma);


    return gamma;
}
