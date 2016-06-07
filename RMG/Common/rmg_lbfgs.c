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

void rmg_lbfgs (void)
{
    int ion, fpt;
    ION *iptr;
    double *position, *force;

    my_malloc(position, 3*ct.num_ions, double);
    my_malloc(force, 3*ct.num_ions, double);
    fpt = ct.fpt[0];

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

