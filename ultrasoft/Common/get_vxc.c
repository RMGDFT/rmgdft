/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

/****f* QMD-MGDFT/get_vxc.c *****
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
 *   void get_vxc(REAL *rho, REAL *rho_oppo, REAL *rhocore, REAL *vxc)
 *   Top-level driver routine that calls a subfunction to generate
 *   a specific exchange-corrleation potential.
 * INPUTS
 *   rho: Electronic charge density in spin-pairwised calculation, while in spin polarized 
 *        calculation, it's processor's own spin density
 *   rho_oppo: THe opposite spin's charge density in spin polarized calculation     
 *   rhocore: Core charge density if Non-linear core corrections are being used.
 *   vxc: The generated exchange-corrleation potential for the prosessor's own spin
 * OUTPUT
 *   vxc: The generated exchange-corrleation potential for the processor's own spin
 * PARENTS
 *   init.c scf.c
 * CHILDREN
 *   get_xc.c
 * SOURCE
 */


#include "main.h"
#include <float.h>
#include <math.h>



void get_vxc (REAL * rho, REAL * rho_oppo, REAL * rhocore, REAL * vxc)
{

    int idx;
    REAL *exc, *nrho, *nrho_oppo;

    
    if (pct.spin_flag)
    {
    	my_malloc (exc, 3 * FP0_BASIS, REAL);
	nrho_oppo = exc + 2 * FP0_BASIS;
    }
    else
    	my_malloc (exc, 2 * FP0_BASIS, REAL);
    
    nrho = exc + FP0_BASIS;


    /* Add in the nonlinear core charges correction from pseudopotential file */
    if (pct.spin_flag)
    {
        /*In spin polarized calculation,  add half of the nonlinear core charge for both 
	 * processor's own spin density and opposite spin density */
    	for (idx = 0; idx < FP0_BASIS; idx++)
	{
        	nrho[idx] = rhocore[idx] * 0.5 + rho[idx];         
		nrho_oppo[idx] = rhocore[idx] * 0.5 + rho_oppo[idx];
	}
    }
    else
    {
    	for (idx = 0; idx < FP0_BASIS; idx++)
        	nrho[idx] = rhocore[idx] + rho[idx];
    }

    get_xc(nrho, nrho_oppo, vxc, exc, ct.xctype);


    my_free (exc);


}                               /* end get_vxc */



/******/
