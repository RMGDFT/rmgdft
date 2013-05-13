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

void rmg_lbfgs (void)
{
    int ion, fpt;
    ION *iptr;
    rmg_double_t step, mass, magf, dotfv;

    double *position, *force;

    my_malloc(position, 3*ct.num_ions, double);
    my_malloc(force, 3*ct.num_ions, double);
    fpt = ct.fpt[0];
    step = ct.iondt;

    /* Loop over ions */
    for (ion = 0; ion < ct.num_ions; ion++)
    {

        /* Get ion pointer */
        iptr = &ct.ions[ion];
        force[ion * 3 + 0] = iptr->force[fpt][0];
        force[ion * 3 + 1] = iptr->force[fpt][1];
        force[ion * 3 + 2] = iptr->force[fpt][2];

        position[ion * 3 + 0] = iptr->crds[0];
        position[ion * 3 + 1] = iptr->crds[1];
        position[ion * 3 + 2] = iptr->crds[2];
    }

    int num_images = pct.images;
    num_images = 1;
    lbfgs (position, force, ct.num_ions, num_images);

    for (ion = 0; ion < ct.num_ions; ion++)
    {

        iptr = &ct.ions[ion];
        if(iptr->movable )
        {
            iptr->crds[0] = position[ion * 3 + 0] ;
            iptr->crds[1] = position[ion * 3 + 1] ;
            iptr->crds[2] = position[ion * 3 + 2] ;

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
        }

    }


}                               /* end rmg_fastrelax */

