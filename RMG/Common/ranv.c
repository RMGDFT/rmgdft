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
#include "main.h"


void ranv (void)
{
    int i, ii, N, ion;
    ION *iptr;
    double c, p, tmass, vtot[3], ek, scale;
    double kB, mass;

    /* init the random number generator */
    srand (123);

    /* count up moving atoms */
    N = 0;
    for (ion = 0; ion < ct.num_ions; ion++)
    {

        if (ct.ions[ion].movable)
            N++;

    }                           /* end for */

    /* define Boltzmann param */
    kB = 1.0 / (11605.0 * Ha_eV);

    if (ct.nose.temp > 0.0)
    {

        tmass = 0;
        vtot[0] = vtot[1] = vtot[2] = 0.0;

        /* randomize the velocity for each atom */
        for (ion = 0; ion < ct.num_ions; ion++)
        {

            /* Get ion pointer */
            iptr = &ct.ions[ion];

            /* Get ionic mass */
            mass = ct.sp[iptr->species].atomic_mass * mu_me;

            if (iptr->movable)
            {
                c = -kB * ct.nose.temp / mass;
                tmass += mass;
                for (i = 0; i < 3; i++)
                {

                    ii = rand ();
                    p = (double) ii / (double) RAND_MAX;

                    iptr->velocity[i] = sqrt (c * log (p));

                    if (p > 0.5)
                        iptr->velocity[i] = -iptr->velocity[i];
                    vtot[i] += mass * iptr->velocity[i];
                }
            }                   /* end of if */

        }                       /* end for ion */

        /* scale set velocity by total mass */
        vtot[0] = vtot[0] / tmass;
        vtot[1] = vtot[1] / tmass;
        vtot[2] = vtot[2] / tmass;

        /* subtract out the center of mass velocity */
        for (ion = 0; ion < ct.num_ions; ion++)
        {

            /* Get ion pointer */
            iptr = &ct.ions[ion];

            if (iptr->movable)
            {
                iptr->velocity[0] -= vtot[0];
                iptr->velocity[1] -= vtot[1];
                iptr->velocity[2] -= vtot[2];
            }
        }                       /* end of for ion */

        /* find kinetic energy */
        ek = 0.0;
        for (ion = 0; ion < ct.num_ions; ion++)
        {

            /* Get ion pointer */
            iptr = &ct.ions[ion];

            /* Get ionic mass */
            mass = ct.sp[iptr->species].atomic_mass * mu_me;

            if (iptr->movable)
            {
                ek += iptr->velocity[0] * iptr->velocity[0] * mass;
                ek += iptr->velocity[1] * iptr->velocity[1] * mass;
                ek += iptr->velocity[2] * iptr->velocity[2] * mass;
            }

        }                       /* end of for ion */

        ek = ek * 0.5;
        scale = sqrt (1.5 * N * kB * ct.nose.temp / ek);

        ek = 0.0;
        /* rescale to correct temperature */
        for (ion = 0; ion < ct.num_ions; ion++)
        {

            /* Get ion pointer */
            iptr = &ct.ions[ion];

            /* Get ionic mass */
            mass = ct.sp[iptr->species].atomic_mass * mu_me;

            if (iptr->movable)
            {
                iptr->velocity[0] *= scale;
                iptr->velocity[1] *= scale;
                iptr->velocity[2] *= scale;
                ek += iptr->velocity[0] * iptr->velocity[0] * mass;
                ek += iptr->velocity[1] * iptr->velocity[1] * mass;
                ek += iptr->velocity[2] * iptr->velocity[2] * mass;
            }

        }                       /* end of for ion */

        ct.ionke = 0.5 * ek;

    }                           /* end of if */

}                               /* end ranv */


/******/
