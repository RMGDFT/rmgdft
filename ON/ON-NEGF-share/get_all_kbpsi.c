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
#include "init_var.h"


void get_local_kbpsi (STATE *st1, double *psi, double *kbpsi);


void get_all_kbpsi (STATE *states1, STATE * states)
{
    int st1, idx;
    int size;
    int ion, ion1, ip;
    double time1, time2;

    time1 = my_crtc ();

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
        get_kbpsi (&states1[st1], &kbpsi[idx]);
    }

/*  print_sum(size, kbpsi, "kbpsi sum get_all_kbpsi "); 
*/

    /* scale the kbpsi by ct.ions[ion].pd[ip] */
/*	for(st1=ct.state_begin;st1 < ct.state_end;st1++) 
*	{ 
*		idx = (st1-ct.state_begin) * pct.n_ion_center * ct.max_nl; 
*		for(ion1 =0; ion1 < pct.n_ion_center; ion1++) 
*		{ 
*			ion = pct.ionidx[ion1]; 
*			for(ip = 0; ip < ct.max_nl; ip++) 
*			{ 
*				kbpsi[idx] *= ct.ions[ion].pd[ip]; 
*				idx++; 
*			} 
*		} 
*	} 
*/

    time2 = my_crtc ();
    rmg_timings (NL_TIME, (time2 - time1));

}


/*
    get_kbpsi:
    
       Loop over projectors for the ions and evaluate 
       <KB|psi> for the states represented in sp1.
       Store the result in the right element of the array kbpsi.
    
*/
void get_kbpsi (STATE *sp1, double *kbpsi_one_state)
{
    int ixx, iyy, izz;
    double time1, time2;
    double *psi;

    time1 = my_crtc ();

    psi = sp1->psiR;

    get_local_kbpsi (sp1, psi, kbpsi_one_state);

    time2 = my_crtc ();
    rmg_timings (NL_TIME, (time2 - time1));

}


/*
    Get <KB|psi> for one ion and one function psi 
*/
void get_local_kbpsi (STATE *st1, double *psi, double *kbpsi_one_state)
{

    int incx = 1, stop, ip, idx;
    double *prjptr;
    int *pidx;
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
                = ct.vel * dot_product_orbit_nl (st1, ion, psi, prjptr);
            prjptr += ct.max_nlpoints;
        }
    }
}
