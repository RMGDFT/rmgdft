/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
/*
                        get_all_partial_kbpsi.c

    Calculating <partial_beta_x|psi>, <partial_beta_y|psi>
    and <partial_beta_z|psi>

    stores result in work array partial_kbpsi_x,
    partial_kbpsi_y and partial_kbpsi_z.

*/
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "main_on.h"


void get_all_partial_kbpsi(STATE * states)
{
    int st1, idx1, idx2, idx;
    int size;
    int ion, ion1, ip;
    double *psi, *prjptr_x, *prjptr_y, *prjptr_z;

    double time1, time2;
    time1 = my_crtc();

    /* size of the <psi|kb> in each processor */
    size = ct.state_per_proc * pct.n_ion_center * ct.max_nl;

    /* Zero out the projection array */
    for (idx = 0; idx < size; idx++)
    {
        partial_kbpsi_x[idx] = 0.;
        partial_kbpsi_y[idx] = 0.;
        partial_kbpsi_z[idx] = 0.;
    }

    for (st1 = ct.state_begin; st1 < ct.state_end; st1++)
    {
        idx1 = (st1 - ct.state_begin) * pct.n_ion_center * ct.max_nl;
        psi = states[st1].psiR;

        prjptr_x = projectors_x;
        prjptr_y = projectors_y;
        prjptr_z = projectors_z;

        for (ion1 = 0; ion1 < pct.n_ion_center; ion1++)
        {
            idx2 = ion1 * ct.max_nl;
            ion = pct.ionidx[ion1];

            /* Loop over projectors for this ion and evaluate */
            /* <KB|psi>                                   */
            for (ip = 0; ip < pct.prj_per_ion[ion]; ip++)
            {
                idx = idx1 + idx2 + ip;
                partial_kbpsi_x[idx] =
                    ct.vel * dot_product_orbit_nl(&states[st1], ion, psi, prjptr_x);
                partial_kbpsi_y[idx] =
                    ct.vel * dot_product_orbit_nl(&states[st1], ion, psi, prjptr_y);
                partial_kbpsi_z[idx] =
                    ct.vel * dot_product_orbit_nl(&states[st1], ion, psi, prjptr_z);

                prjptr_x += ct.max_nlpoints;
                prjptr_y += ct.max_nlpoints;
                prjptr_z += ct.max_nlpoints;
            }
        }

    }

    time2 = my_crtc();
    rmg_timings(PARTIAL_KBPSI, time2 - time1);

}
