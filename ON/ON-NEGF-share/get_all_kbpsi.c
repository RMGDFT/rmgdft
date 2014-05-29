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
#include "main.h"
#include "prototypes_on.h"


void get_local_kbpsi (STATE *st1, double *psi, double *kbpsi,
        ION_ORBIT_OVERLAP *, rmg_double_t *);
void get_kbpsi (STATE *sp1, double *kbpsi_one_state,
        ION_ORBIT_OVERLAP *ion_orbit_overlap_region_nl, rmg_double_t *projectors);


void get_all_kbpsi (STATE *states1, STATE * states, 
        ION_ORBIT_OVERLAP *ion_orbit_overlap_region_nl, 
        rmg_double_t *projectors, rmg_double_t *kbpsi)
{
    int st1, idx;
    int size;
    int ion, ion1, ip;


    /* size of the <psi|kb> in each processor */
    size = ct.state_per_proc * pct.n_ion_center * ct.max_nl;

    /* Zero out the projection array */
    for (idx = 0; idx < size; idx++)
    {
        kbpsi[idx] = 0.;
    }

    for (st1 = ct.state_begin; st1 < ct.state_end; st1++)
    {
        idx = (st1 - ct.state_begin) * pct.n_ion_center * ct.max_nl;
        get_kbpsi (&states1[st1], &kbpsi[idx], ion_orbit_overlap_region_nl, projectors);
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
        ION_ORBIT_OVERLAP *ion_orbit_overlap_region_nl, rmg_double_t *projectors)
{
    int ixx, iyy, izz;
    double *psi;


    psi = sp1->psiR;

    get_local_kbpsi (sp1, psi, kbpsi_one_state,
            ion_orbit_overlap_region_nl, projectors);


}


/*
    Get <KB|psi> for one ion and one function psi 
*/
void get_local_kbpsi (STATE *st1, double *psi, double *kbpsi_one_state,
        ION_ORBIT_OVERLAP *ion_orbit_overlap_region_nl, rmg_double_t *projectors)
{

    int ip;
    double *prjptr;
    int ion, ion2;


    prjptr = projectors;

    for (ion2 = 0; ion2 < pct.n_ion_center; ion2++)
    {
        ion = pct.ionidx[ion2];
        /* Loop over projectors for this ion and evaluate */
        /* <KB|psi>                                   */
        for (ip = 0; ip < pct.prj_per_ion[ion]; ip++)
        {
            kbpsi_one_state[ion2 * ct.max_nl + ip]
                = get_vel() * dot_product_orbit_nl (st1, ion, psi, prjptr,
                        ion_orbit_overlap_region_nl);
            prjptr += ct.max_nlpoints;
        }
    }
}
