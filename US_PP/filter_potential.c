/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

/****f* QMD-MGDFT/init_psp.c *****
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
 *   void init_psp(void)
 *   Initializes radial Kleinman-Bylander projectors and
 *   calculates the normalization coefficient.
 *
 *   Also subtracts off the potential due to the compensating charges from
 *   the local potential.
 * INPUTS
 *   nothing
 * OUTPUT
 *   projectors are stored in structure SPECIES (see main.h)
 * PARENTS
 *   init.c
 * CHILDREN
 *   rft1.c radiff.c 
 * SOURCE
 */


#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "main.h"

#define SMALL 1.e-35


void filter_potential (REAL *potential, REAL *r, int rg_points, REAL rmax, REAL offset, REAL parm, REAL* potential_lgrid, 
	REAL *rab, int l_value, REAL dr, REAL  gwidth, int lgrid_points, REAL rcut, REAL rwidth, REAL * drpotential_lgrid)
{

    REAL *work, rdist, t1, exp_fac;
    int idx, der_flag = 0;

    if(drpotential_lgrid)
	der_flag = 1;

    if (ct.mask_function)
	apply_mask_function(potential, r, rg_points, rmax, offset);


    /* Transform to g-space and filter it */
    rft1 (parm, potential, r, potential_lgrid, rab, rg_points, l_value, dr,
	    gwidth, lgrid_points);

    if (ct.mask_function)
	backout_mask_function (potential_lgrid, dr, lgrid_points,  rmax);

    /*Evaluate radial derivative, if requested*/
    if (der_flag)
    {
	my_malloc(work, lgrid_points, REAL);

	for (idx = 0; idx < lgrid_points; idx++)
	    work[idx] = dr * ((REAL) idx);

	radiff (potential_lgrid, drpotential_lgrid, work, lgrid_points, 0.0);

	/* Fix up the first point */
	drpotential_lgrid[1] = 2.0 * drpotential_lgrid[2] - drpotential_lgrid[3];
	drpotential_lgrid[0] = 2.0 * drpotential_lgrid[1] - drpotential_lgrid[2];

	my_free (work);
    }

    /*Fix up first point in filtered potential*/
    potential_lgrid[0] = 2.0 * potential_lgrid[1] - potential_lgrid[2];


    /* Without mask function, result needs to be dampened by gaussian starting at rcut*/
    if (!ct.mask_function)
    {
	rdist = 0.0;
	for (idx = 0; idx < lgrid_points; idx++)
	{
	    if (rdist > rcut)
	    {
		t1 = (rdist - rcut) / rcut;
		exp_fac = exp (-rwidth * t1 * t1);

		potential_lgrid[idx] *= exp_fac;

		if (fabs (potential_lgrid[idx]) < SMALL)
		    potential_lgrid[idx] = 0.0;

		if (der_flag)
		{
		    drpotential_lgrid[idx] *= exp_fac;

		    if (fabs (drpotential_lgrid[idx]) < SMALL)
			drpotential_lgrid[idx] = 0.0;
		}

	    }               /* end if */

	    rdist += dr;

	}
    }



} 

/******/
