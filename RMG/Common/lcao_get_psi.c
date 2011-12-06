/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "main.h"


void lcao_get_psi (STATE * states)
{

    int ion, idx, ip, l, m, state_count, nspin, st;
    SPECIES *sp;
    ION *iptr;
    REAL *psi, occupancy, occ_per_wave, charge_count;

    nspin = ct.spin_flag+1;

    /*The first index is due to k-point*/
    psi = &states[0].psiR[0];


    /* Loop over ions */
    state_count = 0;
    charge_count = 0.0;

    for (ion = 0; ion < ct.num_ions; ion++)
    {
        /* Generate ion pointer */
        iptr = &ct.ions[ion];

        /* Get species type */
        sp = &ct.sp[iptr->species];

	/*Make sure that the wavefunctions have been read*/
	if (!sp->num_atomic_waves)
	    error_handler("No initial wavefunctions for ion %d, most likely the PP file does not have them", ion);
        
	
	/*Loop over atomic wavefunctions for given ion*/
	for (ip = 0; ip < sp->num_atomic_waves; ip++)
	{
	    l = sp->atomic_wave_l[ip];
	    occupancy = sp->atomic_wave_oc[ip];

	    /*Make sure that occupancy will not be too high, may happen if system is positively charged*/
	    if ((charge_count + occupancy ) > ct.nel)
		occupancy = ct.nel - charge_count;

	    occ_per_wave = occupancy / ((REAL) (2*l+1));

	    /*Loop over all m values for given l and get wavefunctions */
	    for (m=0; m < 2*l+1; m++)
	    {
		for (idx = 0; idx < P0_BASIS; idx++)
		    psi[idx] = 0.0;

		lcao_get_awave(psi, iptr, ip, l, m);

		if (state_count == ct.num_states) 
		    error_handler("Error too few states for LCAO start");
		
		/*Assign occupation to that from PP file*/
		/*This needs to be changed from spin*/
		states[state_count].occupation[0] = occ_per_wave;

		psi += P0_BASIS;
		state_count ++;

		charge_count +=  occ_per_wave;
	    }
	    
	}

    }


    /*Initialize any additional states to random start*/
    if ( ct.num_states > state_count)
    {
	int ix, iy, iz;
	int xoff, yoff, zoff;
	REAL xrand[2 * NX_GRID], yrand[2 * NY_GRID], zrand[2 * NZ_GRID];
	long idum;
	STATE *state_p;
	
	/* Initialize the random number generator */
	idum = 3356;
	rand0 (&idum);

	occ_per_wave = (ct.nel - charge_count) / (ct.num_states - state_count);

	printf("\n occ_per_wave for unoccupied states is %f, ct.nel %f, charge_count %f, ct.num_states %d, state_count %d", occ_per_wave, ct.nel, charge_count, ct.num_states, state_count); 
	
	for (st = state_count; st < ct.num_states; st++)
	{

	    /* Generate x, y, z random number sequences */
	    for (idx = 0; idx < ct.psi_nxgrid; idx++)
		xrand[idx] = rand0 (&idum) - 0.5;
	    for (idx = 0; idx < ct.psi_nygrid; idx++)
		yrand[idx] = rand0 (&idum) - 0.5;
	    for (idx = 0; idx < ct.psi_nzgrid; idx++)
		zrand[idx] = rand0 (&idum) - 0.5;

#if !GAMMA_PT
	    for (idx = ct.psi_nxgrid; idx < 2 * ct.psi_nxgrid; idx++)
		xrand[idx] = rand0 (&idum) - 0.5;
	    for (idx = ct.psi_nygrid; idx < 2 * ct.psi_nygrid; idx++)
		yrand[idx] = rand0 (&idum) - 0.5;
	    for (idx = ct.psi_nzgrid; idx < 2 * ct.psi_nzgrid; idx++)
		zrand[idx] = rand0 (&idum) - 0.5;
#endif

	    state_p = &states[st];

	    /* If random start and Fermi occupation, start with
	       each state equally occupied  */

	    for (idx = 0; idx < nspin; idx++)
		state_p->occupation[idx] = occ_per_wave / nspin ;	


	    idx = 0;
	    for (ix = 0; ix < PX0_GRID; ix++)
	    {

		for (iy = 0; iy < PY0_GRID; iy++)
		{

		    for (iz = 0; iz < PZ0_GRID; iz++)
		    {

			state_p->psiR[idx] = xrand[xoff + ix] * yrand[yoff + iy] * zrand[zoff + iz];
			state_p->psiR[idx] = state_p->psiR[idx] * state_p->psiR[idx];
#if !GAMMA_PT
			state_p->psiI[idx] =
			    xrand[ct.psi_nxgrid + xoff + ix] * yrand[ct.psi_nygrid + yoff +
			    iy] * zrand[ct.psi_nzgrid + zoff +
			    iz];
			state_p->psiI[idx] =  state_p->psiI[idx] *  state_p->psiI[idx];

#endif
			idx++;

		    }               /* end for */
		}                   /* end for */
	    }                       /* end for */


	}                           /* end for */

    }


}
/******/
