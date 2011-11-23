/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "main.h"


void lcao_init_psi (STATE * states)
{

    int ion, idx, ip, l, m;
    SPECIES *sp;
    ION *iptr;
    REAL *psi, occupancy;

    /*The first index is due to k-point*/
    psi = &states[0].psiR[0];


    /* Loop over ions */
    for (ion = 0; ion < ct.num_ions; ion++)
    {
        /* Generate ion pointer */
        iptr = &ct.ions[ion];

        /* Get species type */
        sp = &ct.sp[iptr->species];
        
	
	/*Loop over atomic wavefunctions for given ion*/
	for (ip = 0; ip < sp->num_atomic_waves; ip++)
	{
	    l = sp->atomic_wave_l[ip];
	    occupancy = sp->atomic_wave_oc[ip];

	    /*Loop over all m values for given l and get wavefunctions */
	    for (m=0; m < 2*l+1; m++)
	    {
		for (idx = 0; idx < P0_BASIS; idx++)
		    psi[idx] = 0.0;

		get_awave(psi, iptr, ip, l, m);

		psi += P0_BASIS;

	    }
	    

	}


    }



}
/******/
