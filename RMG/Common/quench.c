/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

/****f* QMD-MGDFT/quench.c *****
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
 *   void quench(STATE *states, REAL *vxc, REAL *vh, REAL *vnuc, 
 *               REAL *rho, REAL *rhocore, REAL *rhoc)
 *   For a fixed atomic configuration, quench the electrons to find 
 *   the minimum total energy 
 * INPUTS
 *   states: point to orbital structure (see main.h)
 *   vxc:    exchange correlation potential
 *   vh:     Hartree potential
 *   vnuc:   Pseudopotential 
 *   rho:    total valence charge density
 *   rhocore: core chare density only for non-linear core correction
 *   rhoc:   Gaussian compensating charge density
 * OUTPUT
 *   states, vxc, vh, rho are updated
 * PARENTS
 *   cdfastrlx.c fastrlx.c main.c
 * CHILDREN
 *   scf.c force.c get_te.c subdiag.c get_ke.c
 * SOURCE
 */



#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "main.h"

bool quench (STATE * states, REAL * vxc, REAL * vh, REAL * vnuc, REAL * rho,
             REAL * rho_oppo, REAL * rhocore, REAL * rhoc)
{

    bool CONVERGED;
    int numacc = 1, ic;
    /*int ist, ik;
       REAL KE; */

    /* ---------- begin scf loop ---------- */
    

    for (ct.scf_steps = 0, CONVERGED = false;
         ct.scf_steps < ct.max_scf_steps && !CONVERGED; ct.scf_steps++, ct.total_scf_steps++)
    {

        if (pct.imgpe == 0)
            printf ("\n\nquench: ------ [md: %d/%d  scf: %d/%d] ------\n",
                    ct.md_steps, ct.max_md_steps, ct.scf_steps, ct.max_scf_steps);


        /* perform a single self-consistent step */
        CONVERGED = scf (states, vxc, vh, vnuc, rho, rho_oppo, rhocore, rhoc);


	get_te (rho, rho_oppo, rhocore, rhoc, vh, vxc, states, !ct.scf_steps);

	/* output the eigenvalues with occupations */
	if (ct.write_eigvals_period)
	{
	    if (ct.scf_steps % ct.write_eigvals_period == 0)
	    {
		if (pct.imgpe == 0)
		{
		    output_eigenvalues (states, 0, ct.scf_steps);
		    printf ("\nTotal charge in supercell = %16.8f\n", ct.tcharge);
		}
	    }
	}


    }

    /* ---------- end scf loop ---------- */

    if (CONVERGED)
    {
	printf ("\n");
	progress_tag ();
	printf ("potential convergence has been achieved. stopping ...\n");
	    
	/*Write PDOS if converged*/
	if (ct.pdos_flag)
	    get_pdos (states, ct.Emin, ct.Emax, ct.E_POINTS);
    }

    printf ("\n");
    progress_tag ();
    printf ("final total energy = %14.7f Ha\n", ct.TOTAL);



    /* output final eigenvalues with occupations */

    output_eigenvalues (states, 0, ct.scf_steps);
    printf ("\nTotal charge in supercell = %16.8f\n", ct.tcharge);

    wvfn_residual (states);




    /*When running MD, force pointers need to be rotated before calculating new forces */
    if ((ct.forceflag == MD_CVE) || (ct.forceflag == MD_CVT) || (ct.forceflag == MD_CPT))
    {

	/* rotate the force pointers */
	switch (ct.mdorder)
	{
	    case ORDER_2:
		numacc = 1;
		break;
	    case ORDER_3:
		numacc = 2;
		break;
	    case ORDER_5:
		numacc = 4;
		break;
	}
	for (ic = (numacc - 1); ic > 0; ic--)
	{
	    ct.fpt[ic] = ct.fpt[ic - 1];
	}
	ct.fpt[0] = ct.fpt[numacc - 1] + 1;
	if (ct.fpt[0] > (numacc - 1) || numacc == 1)
	    ct.fpt[0] = 0;

    }


    /* compute the forces */
    /* Do not calculate forces for quenching when we are not converged */
    if (CONVERGED || (ct.forceflag != MD_QUENCH))
	force (rho, rho_oppo, rhoc, vh, vxc, vnuc, states);

    /* output the forces */
    if (pct.imgpe == 0)
	write_force ();

    return CONVERGED;


}                               /* end quench */




/******/
