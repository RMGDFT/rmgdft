/*
 *
 * Copyright 2014 The RMG Project Developers. See the COPYRIGHT file 
 * at the top-level directory of this distribution or in the current
 * directory.
 * 
 * This file is part of RMG. 
 * RMG is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * any later version.
 *
 * RMG is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
*/




#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "main.h"

void fire (double *step, double step_max, double f_inc, double f_dec, int n_min, int *n_count )
{
    int ion, fpt;
    ION *iptr;
    double  magf = 0.0, magv = 0.0, dotfv, p=0.0, unitf[3], eff_mass = 6.0 * mu_me;

    /*Fire parameters*/
    double alpha_start = 0.1;
    double f_alpha = 0.99;
    static double alpha = 0.1;

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
