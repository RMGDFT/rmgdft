/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

/****f* QMD-MGDFT/get_rho.c *****
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
 *   void get_rho(STATE *states, REAL *rho)
 *   Generates new charge density and mix linearly with old one.
 * INPUTS
 *   states:  point to orbital structure (see main.h)
 *   rho:     old charge density
 * OUTPUT
 *   rho:    updated charge density
 * PARENTS
 *   scf.c
 * CHILDREN
 *   gather_psi.c symmetrize_rho.f
 * SOURCE
 */


#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "main.h"


void mix_rho (REAL * new_rho, REAL * rho, REAL *rhocore, int length, int length_x, int length_y, int length_z)
{
    REAL t1, min, min2, nspin = (ct.spin_flag + 1.0);
    int step, idx, inc = 1;
    static REAL **rhohist=NULL, **residhist=NULL;

    /*Linear Mixing*/
    if (verify("charge_mixing_type","Linear"))
    {
	
	/* Scale old charge density first*/
	t1 = 1.0 - ct.mix;
	QMD_dscal (length, t1, rho, inc); 

	/*Add the new density*/
	QMD_daxpy (length, ct.mix, new_rho, inc, rho, inc);
    }
    else {
	if (verify("charge_mixing_type","Pulay"))
	{
	    step = ct.scf_steps;

	    if (ct.charge_pulay_refresh)
		step = ct.scf_steps % ct.charge_pulay_refresh;

	    /*Use pulay mixing, result will be in rho*/
	    pulay_rho(step, length, length_x, length_y, length_z, new_rho, rho, ct.charge_pulay_order, &rhohist, &residhist, ct.charge_pulay_special_metrics, ct.charge_pulay_special_metrics_weight);
	    
	}
	    
    }


    /*Find charge minimum */
    min = ZERO;
    min2 = ZERO;
    for (idx = 0; idx < length; idx++)
    {
        if (rho[idx] < min)
            min = rho[idx];
        
	/*Here we check charge density with rhocore added*/
	if ((rho[idx] + rhocore[idx] / nspin) < min2)
         	min2 = rho[idx] + rhocore[idx] / nspin;
    }




    /*Find absolute minimum from all PEs */
    min = real_min_all (min, pct.img_comm);
    min2 = real_min_all (min2, pct.img_comm);

    if (min < ZERO)
    {
        printf ("\n\n Charge density is NEGATIVE after interpolation, minimum is %e", min);
        printf ("\n Minimum charge density with core charge added is %e", min2);
    }
}                               /* end get_rho */

/******/
