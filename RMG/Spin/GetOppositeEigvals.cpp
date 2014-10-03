/************************** SVN Revision Information **************************
 **    $Id: get_opposite_eigvals.c 2009 2013-05-13 16:55:41Z ebriggs $    **
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
 *   Obtaines eigenvlaues for wavefunctions with opposite spin
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
 *   get_vxc.c get_vh.c mg_eig_state.c ortho_full.c fill.c get_rho.c
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


template void GetOppositeEigvals(Kpoint<double> **);
template void GetOppositeEigvals(Kpoint<std::complex<double> > **);

template <typename KpointType>
void GetOppositeEigvals (Kpoint<KpointType> ** Kptr)
{

    MPI_Status status; 
    double *eigval_sd, *eigval_rv;   
    int st, st1, kpt;


    /* allocate memory for eigenvalue send array and receive array */
    eigval_sd = new double[2 * ct.num_kpts * ct.num_states];
    eigval_rv = eigval_sd + ct.num_kpts * ct.num_states;


    /*Prepare the sending buffer of eigenvalues */
    st = 0;
    for (kpt =0; kpt < ct.num_kpts; kpt++)
    {
	for (st1 = 0; st1 < ct.num_states; st1++)
	{	
	    eigval_sd[st] = Kptr[kpt]->Kstates[st1].eig[0];
	    st += 1;
	}	
    }


    /*Communicate for spin up and spin down energy eigenvalues*/    
    MPI_Sendrecv(eigval_sd, st, MPI_DOUBLE, (pct.spinpe+1)%2, pct.gridpe,
	    eigval_rv, st, MPI_DOUBLE, (pct.spinpe+1)%2, pct.gridpe, pct.spin_comm, &status);


    /* Unpack the received eigenvalue to state structure */
    st = 0;
    for (kpt =0; kpt < ct.num_kpts; kpt++)
    {
	for (st1 = 0; st1 < ct.num_states; st1++)
	{	
	    Kptr[kpt]->Kstates[st1].eig[1] = eigval_rv[st];
	    st += 1;

	}	
    } 


    delete [] eigval_sd;
}                               /* end scf */


