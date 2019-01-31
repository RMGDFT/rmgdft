/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
/*
                        get_all_kbpsi.c

    Apply non-local operator to a state on a given grid level and 
    stores result in the work array.
*/
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include "params.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "RmgTimer.h"
#include "transition.h"
#include "prototypes_on.h"
#include "init_var.h"
#include "blas.h"
#include "Kbpsi.h"


void get_local_kbpsi (STATE *st1, double *psi, double *kbpsi,
        ION_ORBIT_OVERLAP *, double *);
void get_kbpsi (STATE *sp1, double *kbpsi_one_state,
        ION_ORBIT_OVERLAP *ion_orbit_overlap_region_nl, double *projectors);


void GetAllKbpsi (STATE *states1, STATE * states, 
        ION_ORBIT_OVERLAP *ion_orbit_overlap_region_nl, 
        double *projectors, double *kbpsi)
{
    int st1, idx;
    int size;


    /* size of the <psi|kb> in each processor */
    size = ct.state_per_proc * pct.n_ion_center * ct.max_nl;

    /* Zero out the projection array */
    for (idx = 0; idx < size; idx++)
    {
        kbpsi[idx] = 0.;
    }

#pragma omp parallel private(st1)
{
#pragma omp for schedule(dynamic) nowait
    for (st1 = ct.state_begin; st1 < ct.state_end; st1++)
    {
        idx = (st1 - ct.state_begin) * pct.n_ion_center * ct.max_nl;
        get_kbpsi (&states1[st1], &kbpsi[idx], ion_orbit_overlap_region_nl, projectors);
    }
}

/*  print_sum(size, kbpsi, "kbpsi sum get_all_kbpsi "); 
*/



}


/*
    get_kbpsi:
    
       Loop over projectors for the ions and evaluate 
       <KB|psi> for the states represented in sp1.
       Store the result in the right element of the array kbpsi.
    
*/
void get_kbpsi (STATE *sp1, double *kbpsi_one_state,
        ION_ORBIT_OVERLAP *ion_orbit_overlap_region_nl, double *projectors)
{
    double *psi;


    psi = sp1->psiR;

    get_local_kbpsi (sp1, psi, kbpsi_one_state,
            ion_orbit_overlap_region_nl, projectors);


}


/*
    Get <KB|psi> for one ion and one function psi 
*/
void get_local_kbpsi (STATE *st1, double *psi, double *kbpsi_one_state,
        ION_ORBIT_OVERLAP *ion_orbit_overlap_region_nl, double *projectors)
{

    double *prjptr;
    int ion;


    prjptr = projectors;

    for (unsigned int ion2 = 0; ion2 < pct.n_ion_center; ion2++)
    {
        ion = pct.ionidx[ion2];
        /* Loop over projectors for this ion and evaluate */
        /* <KB|psi>                                   */
        DotProductOrbitNl (st1, ion, psi, prjptr, ion_orbit_overlap_region_nl, 
                pct.prj_per_ion[ion], &kbpsi_one_state[ion2 * ct.max_nl]);
        prjptr += pct.prj_per_ion[ion] * ct.max_nlpoints;
    }
}
