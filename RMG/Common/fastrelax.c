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

void fastrelax (REAL *dt, REAL dt_max, REAL dt_inc, REAL dt_dec, int n_min)
{
    int ion, fpt, which = 0, count = 0;
    ION *iptr;
    REAL mass, magf, dotfv;
    REAL max_move = 0.0, avg_move = 0.0, rms_move = 0.0;
    REAL move_x, move_y, move_z, move_sq, move, p = 0.0;
    static int n_count;


    fpt = ct.fpt[0];

    if (verify("relax_dynamic_timestep",&SET))
    {
    
        for (ion = 0; ion < ct.num_ions; ion++)
        {

        /* Get ion pointer */
        iptr = &ct.ions[ion];
        
        /* Dot product of f and v */
        p += iptr->force[fpt][0] * iptr->velocity[0] +
            iptr->force[fpt][1] * iptr->velocity[1] + iptr->force[fpt][2] * iptr->velocity[2];
        }
    
        if (p < 0) 
        {
            *dt *= dt_dec;
            n_count = 0;
    
            printf ("\n");
            progress_tag ();
            printf("p(%e) < 0, decreasing timestep to %f", p, *dt);  
        }
        else
        {
            if (n_count > n_min)
            {   
                *dt *= dt_inc;

                if (*dt > dt_max)
                    *dt = dt_max;
            }   

            n_count ++;
            
            printf ("\n");
            progress_tag ();
            printf("timestep:%f n_count:%d", *dt, n_count);  
        }



    }

    /* Loop over ions */
    for (ion = 0; ion < ct.num_ions; ion++)
    {

        /* Get ion pointer */
        iptr = &ct.ions[ion];


        /* Use either actual ionic mass or equal mass for all atoms*/
	if (ct.relax_mass == 0)
	    mass = ct.sp[iptr->species].atomic_mass * mu_me;
	else
	    mass =  12.0 * mu_me;
	


        /* Magnitudes of f and v */
        magf = iptr->force[fpt][0] * iptr->force[fpt][0] +
            iptr->force[fpt][1] * iptr->force[fpt][1] + iptr->force[fpt][2] * iptr->force[fpt][2];


        /* Dot product of f and v */
        dotfv = iptr->force[fpt][0] * iptr->velocity[0] +
            iptr->force[fpt][1] * iptr->velocity[1] + iptr->force[fpt][2] * iptr->velocity[2];
            
        iptr->velocity[0] = *dt * iptr->force[fpt][0] / mass;
        iptr->velocity[1] = *dt * iptr->force[fpt][1] / mass;
        iptr->velocity[2] = *dt * iptr->force[fpt][2] / mass;

        
        if (dotfv >= 1.0e-12)
        {

            iptr->velocity[0] += dotfv * iptr->force[fpt][0] / magf;
            iptr->velocity[1] += dotfv * iptr->force[fpt][1] / magf;
            iptr->velocity[2] += dotfv * iptr->force[fpt][2] / magf;

        }
    }                           /* end for */

    /*Move ions by dt*velocity*/
    move_ions(*dt);
	
}                               /* end rmg_fastrelax */

/******/
