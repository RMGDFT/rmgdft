/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

/****f* QMD-MGDFT/scf.c *****
 * NAME
 *   Ab initio real space code with multigrid acceleration
 *   Quantum molecular dynamics package.
 *   Version: 2.1.5
 * COPYRIGHT
 *   Copyright (C) 1995  Emil Briggs
 *   Copyright (C) 1998  Emil Briggs, Charles Brabec, Mark Wensell, 
 *                       Dan Sullivan, Chris Rapcewicz, Jerzy Bernholc
 *   Copyright (C) 2001  Emil Briggs, Wenchang Lu,
 *                       Marco Buongiorno Nardelli,Charles Brabec, 
 *                       Mark Wensell,Dan Sullivan, Chris Rapcewicz,
 *                       Jerzy Bernholc
 * FUNCTION
 *   Obtaines occupations for wavefunctions with opposite spin
 * INPUTS
 *   states: points to orbital structure
 *   vxc: exchange correlation potential
 *   vh:  Hartree potential
 *   vnuc: pseudopotential
 *   rho: total valence charge density
 *   rhocore:  core charge density
 *   rhoc: Gaussian compensating charge density
 * OUTPUT
 *   states, vxc, vh, rho are updated
 *   CONVERGENCE: 1 converged, 0 not
 * PARENTS
 *   cdfastrlx.c fastrlx.c moldyn.c quench.c
 * CHILDREN
 *   
 * SOURCE
 */


#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "main.h"


void get_opposite_occupancies (STATE * states)
{

    MPI_Status status; 
    rmg_double_t *occ_sd, *occ_rv;   
    int st, st1, kpt;


    /* allocate memory for occupation send array and receive array */
    my_malloc (occ_sd, 2 * ct.num_kpts * ct.num_states, rmg_double_t);
    occ_rv = occ_sd + ct.num_kpts * ct.num_states;


    /*Prepare the sending buffer of occupancies */
    st = 0;
    for (kpt =0; kpt < ct.num_kpts; kpt++)
    {
	for (st1 = 0; st1 < ct.num_states; st1++)
	{	
	    occ_sd[st] = ct.kp[kpt].kstate[st1].occupation[0];
	    st += 1;
	}	
    }


    /*Communicate for spin up and spin down occupations*/    
    MPI_Sendrecv(occ_sd, st, MPI_DOUBLE, (pct.spinpe+1)%2, pct.gridpe,
	    occ_rv, st, MPI_DOUBLE, (pct.spinpe+1)%2, pct.gridpe, pct.spin_comm, &status);


    /* Unpack the received occupations to state structure */
    st = 0;
    for (kpt =0; kpt < ct.num_kpts; kpt++)
    {
	for (st1 = 0; st1 < ct.num_states; st1++)
	{	
	    ct.kp[kpt].kstate[st1].occupation[1] = occ_rv[st];
	    st += 1;

	}	
    } 


	my_free (occ_sd);
}                               /* end scf */


/******/
