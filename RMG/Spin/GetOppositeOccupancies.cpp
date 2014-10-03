/************************** SVN Revision Information **************************
 **    $Id: get_opposite_occupancies.c 2009 2013-05-13 16:55:41Z ebriggs $    **
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


#include <complex>
#include "const.h"
#include "params.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "rmg_error.h"
#include "State.h"
#include "Kpoint.h"
#include "transition.h"


template void GetOppositeOccupancies(Kpoint<double> **);
template void GetOppositeOccupancies(Kpoint<std::complex<double> > **);

template <typename KpointType>
void GetOppositeOccupancies (Kpoint<KpointType> ** Kptr)
{

    MPI_Status status; 
    double *occ_sd, *occ_rv;   
    int st, st1, kpt;


    /* allocate memory for occupation send array and receive array */
    occ_sd = new double[2 * ct.num_kpts * ct.num_states];
    occ_rv = occ_sd + ct.num_kpts * ct.num_states;


    /*Prepare the sending buffer of occupancies */
    st = 0;
    for (kpt =0; kpt < ct.num_kpts; kpt++)
    {
	for (st1 = 0; st1 < ct.num_states; st1++)
	{	
	    occ_sd[st] = Kptr[kpt]->Kstates[st1].occupation[0];
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
	    Kptr[kpt]->Kstates[st1].occupation[1] = occ_rv[st];
	    st += 1;

	}	
    } 


    delete [] occ_sd;
}                               /* end scf */


/******/
