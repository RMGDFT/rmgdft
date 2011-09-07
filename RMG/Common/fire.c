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

void fire (REAL *step, REAL step_max, REAL f_inc, REAL f_dec, int n_min )
{
    int ion, fpt, which = -1, count = 0;
    ION *iptr;
    REAL  magf = 0.0, magv = 0.0, dotfv, p=0.0, unitf[3], eff_mass = 6.0 * mu_me;
    REAL max_move = 0.0, avg_move = 0.0, rms_move = 0.0;
    REAL move_x, move_y, move_z, move_sq, move;

    /*Fire parameters*/
    REAL alpha_start = 0.1;
    REAL f_alpha = 0.99;
    static REAL alpha = 0.1;
    static int n_count = 0;

    fpt = ct.fpt[0];
	
    //printf ("\n Starting FIRE relax: parameters are step:%f, alpha:%f, n_count is %d", *step, alpha,  n_count);

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
	n_count = 0;
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


	/* Move the ions */
	if (iptr->movable)
	{
	    move_x = *step * iptr->velocity[0];
	    move_y = *step * iptr->velocity[1];
	    move_z = *step * iptr->velocity[2];

	    /*Update coordinates*/
            iptr->crds[0] += move_x;
            iptr->crds[1] += move_y;
            iptr->crds[2] += move_z;

            /* enforce periodic boundary conditions on the ions */
            to_crystal (iptr->xtal, iptr->crds);
            to_cartesian (iptr->xtal, iptr->crds);

	    /*Find maximum, average and RMS displacement*/
	    move_sq = move_x*move_x + move_y*move_y + move_z*move_z;
	    move = sqrt (move_sq);

	    if (move > max_move)
	    {
		max_move = move;
		which = ion;
	    }

	    avg_move += move;
	    count ++;

	    rms_move += move_sq;
	}                       /* end if */
    }
    
    /*Write out displacement info*/
    printf ("\n");
    progress_tag ();
    printf ("Max displacement: %8.5f a0  (ion %d)", max_move, which + 1);
    printf ("\n");
    progress_tag ();
    printf ("Avg displacement: %8.5f a0", avg_move/count);
    printf ("\n");
    progress_tag ();
    printf ("RMS displacement: %8.5f a0", sqrt(rms_move/count));

    if ((p > 0.0) && (n_count > n_min))
    {
	*step *= f_inc;
	if (*step > step_max) *step = step_max;

	alpha *= f_alpha;

    }

    if (p > 0.0) 
    {
	n_count ++;

        printf ("\n");
        progress_tag ();
	printf ("alpha:%7.5f, n_count %d, timestep:%10.5f", alpha,  n_count, *step);
    }

}

/******/
