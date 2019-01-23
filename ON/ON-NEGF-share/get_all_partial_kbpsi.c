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
#include "main.h"
#include "prototypes_on.h"


void get_all_partial_kbpsi(STATE * states, ION_ORBIT_OVERLAP
        *ion_orbit_overlap_region_nl, double *projectors_x, 
        double *projectors_y, double *projectors_z, 
        double *partial_kbpsi_x,
        double *partial_kbpsi_y,double *partial_kbpsi_z)
{
    int st1, idx1, idx2, idx;
    int size;
    int ion, ion1, ip;
    double *psi, *prjptr_x, *prjptr_y, *prjptr_z;


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
            idx = idx1 + idx2;
            DotProductOrbitNl(&states[st1], ion, psi, prjptr_x, ion_orbit_overlap_region_nl, pct.prj_per_ion[ion], &partial_kbpsi_x[idx]);
            DotProductOrbitNl(&states[st1], ion, psi, prjptr_y, ion_orbit_overlap_region_nl, pct.prj_per_ion[ion], &partial_kbpsi_y[idx]);
            DotProductOrbitNl(&states[st1], ion, psi, prjptr_z, ion_orbit_overlap_region_nl, pct.prj_per_ion[ion], &partial_kbpsi_z[idx]);

            prjptr_x += ct.max_nlpoints * pct.prj_per_ion[ion];
            prjptr_y += ct.max_nlpoints * pct.prj_per_ion[ion];
            prjptr_z += ct.max_nlpoints * pct.prj_per_ion[ion];
        }

    }


}
