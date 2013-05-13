/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

/****f* QMD-MGDFT/rmg_fastrelax.c *****
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
 *   void rmg_fastrelax(void)
 *   Performs molecular dynamics using fast relax algorithm 
 * INPUTS
 *   no explict input
 * OUTPUT
 *   atomic coordinates are updated
 * PARENTS
 *   cdfastrlx.c fastrlx.c
 * CHILDREN
 *   to_crystal.c to_cartesian.c
 * SOURCE
 */



#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "main.h"

void fire (rmg_double_t *step, rmg_double_t step_max, rmg_double_t f_inc, rmg_double_t f_dec, int n_min, int *n_count )
{
    int ion, fpt, which = -1, count = 0;
    ION *iptr;
    rmg_double_t  magf = 0.0, magv = 0.0, dotfv, p=0.0, unitf[3], eff_mass = 6.0 * mu_me;
    rmg_double_t max_move = 0.0, avg_move = 0.0, rms_move = 0.0;
    rmg_double_t move_x, move_y, move_z, move_sq, move;

    /*Fire parameters*/
    rmg_double_t alpha_start = 0.1;
    rmg_double_t f_alpha = 0.99;
    static rmg_double_t alpha = 0.1;

    fpt = ct.fpt[0];
	

    /* Loop over ions, calculate p, total magnitude of F and v vectors*/ 
    for (ion = 0; ion < ct.num_ions; ion++)
    {

	/* Get ion pointer */
	iptr = &ct.ions[ion];

	/* Dot product of f and v */
	dotfv = iptr->force[fpt][0] * iptr->velocity[0] +
	    iptr->force[fpt][1] * iptr->velocity[1] + iptr->force[fpt][2] * iptr->velocity[2];
	
	p += dotfv; 
	    
	/* Magnitude of f  */
	magf += iptr->force[fpt][0] * iptr->force[fpt][0] +
		iptr->force[fpt][1] * iptr->force[fpt][1] + iptr->force[fpt][2] * iptr->force[fpt][2];
	
	/* Magnitude of v  */
	magv += iptr->velocity[0] * iptr->velocity[0] +
		iptr->velocity[1] * iptr->velocity[1] + iptr->velocity[2] * iptr->velocity[2];

    }

    magf = sqrt(magf);
    magv = sqrt(magv);
    
    /* For negative p, decrease time step, freeze the system and exit */
    if ((p <= 0.0) && (ct.md_steps > 0))
    {
	*step *= f_dec;
	*n_count = 0;
	alpha = alpha_start;

        printf ("\n");
        progress_tag ();
	printf ("p (%e) <= zero, freezing system, setting timestep to %10.5f", p, *step);
    }
    

    for (ion = 0; ion < ct.num_ions; ion++)
    {
	/* Get ion pointer */
	iptr = &ct.ions[ion];

	/*Unit vector in direction of force*/
	unitf[0] = iptr->force[fpt][0] / magf;
	unitf[1] = iptr->force[fpt][1] / magf;
	unitf[2] = iptr->force[fpt][2] / magf;


	if (p > 0.0)
	{
	    /*Fire velocity update*/
	    iptr->velocity[0] *= (1.0 - alpha); 
	    iptr->velocity[1] *= (1.0 - alpha); 
	    iptr->velocity[2] *= (1.0 - alpha); 

	    iptr->velocity[0] += alpha * unitf[0] * magv;
	    iptr->velocity[1] += alpha * unitf[1] * magv;
	    iptr->velocity[2] += alpha * unitf[2] * magv;
	}

	else
	{
	    iptr->velocity[0] = 0.0;
	    iptr->velocity[1] = 0.0;
	    iptr->velocity[2] = 0.0;
	}

	/*Euler velocity update*/
	iptr->velocity[0] += *step *  iptr->force[fpt][0] / eff_mass;
	iptr->velocity[1] += *step *  iptr->force[fpt][1] / eff_mass;
	iptr->velocity[2] += *step *  iptr->force[fpt][2] / eff_mass;
    }

    /*Move ions by dt*velocity*/
    move_ions(*step);
    

    if ((p > 0.0) && (*n_count > n_min))
    {
	*step *= f_inc;
	if (*step > step_max) *step = step_max;

	alpha *= f_alpha;

    }

    if (p > 0.0) 
    {
	*n_count ++;

        printf ("\n");
        progress_tag ();
	printf ("alpha:%7.5f, n_count %d, timestep:%10.5f", alpha,  *n_count, *step);
    }

}

/******/
