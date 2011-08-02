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

void rmg_fastrelax (void)
{
    int ion, fpt;
    ION *iptr;
    REAL step, mass, magf, dotfv;


    fpt = ct.fpt[0];
    step = ct.iondt;

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

        if (dotfv >= 1.0e-12)
        {

            iptr->velocity[0] = dotfv * iptr->force[fpt][0] / magf +
                step * iptr->force[fpt][0] / mass;

            iptr->velocity[1] = dotfv * iptr->force[fpt][1] / magf +
                step * iptr->force[fpt][1] / mass;

            iptr->velocity[2] = dotfv * iptr->force[fpt][2] / magf +
                step * iptr->force[fpt][2] / mass;

        }
        else
        {

            iptr->velocity[0] = step * iptr->force[fpt][0] / mass;

            iptr->velocity[1] = step * iptr->force[fpt][1] / mass;

            iptr->velocity[2] = step * iptr->force[fpt][2] / mass;

        }                       /* end if */




        /* Move the ion */
        if (iptr->movable)
        {

            iptr->crds[0] += step * iptr->velocity[0];

            iptr->crds[1] += step * iptr->velocity[1];

            iptr->crds[2] += step * iptr->velocity[2];

            /* enforce periodic boundary conditions on the ions */

            to_crystal (iptr->xtal, iptr->crds);
            if (iptr->xtal[0] > ONE)
                iptr->xtal[0] -= ONE;
            if (iptr->xtal[0] < ZERO)
                iptr->xtal[0] += ONE;

            if (iptr->xtal[1] > ONE)
                iptr->xtal[1] -= ONE;
            if (iptr->xtal[1] < ZERO)
                iptr->xtal[1] += ONE;

            if (iptr->xtal[2] > ONE)
                iptr->xtal[2] -= ONE;
            if (iptr->xtal[2] < ZERO)
                iptr->xtal[2] += ONE;

            to_cartesian (iptr->xtal, iptr->crds);

        }                       /* end if */


    }                           /* end for */

}                               /* end rmg_fastrelax */

/******/
