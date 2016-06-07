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


#include "portability.h"
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "main.h"
#include "input.h"

void fastrelax (double *dt, double dt_max, double dt_inc, double dt_dec, int n_min, int *n_count)
{
    int ion, fpt;
    ION *iptr;
    double mass, magf, dotfv, force[3], p = 0.0;
    bool set = true;


    fpt = ct.fpt[0];

    if (verify_boolean("relax_dynamic_timestep",&set))
    {

        for (ion = 0; ion < ct.num_ions; ion++)
        {

            /* Get ion pointer */
            iptr = &ct.ions[ion];

            force[0] = iptr->force[fpt][0];
            force[1] = iptr->force[fpt][1];
            force[2] = iptr->force[fpt][2];

            if (ct.constrainforces)
            {
                force[0] += iptr->constraint.forcemask[0];
                force[1] += iptr->constraint.forcemask[1];
                force[2] += iptr->constraint.forcemask[2];
            }
            
            /* Dot product of f and v */
            p += force[0] * iptr->velocity[0] +
                 force[1] * iptr->velocity[1] +
                 force[2] * iptr->velocity[2];
        }

        if (p < 0) 
        {
            *dt *= dt_dec;
            *n_count = 0;

            printf ("\n");
            progress_tag ();
            printf("p(%e) < 0, decreasing timestep to %f", p, *dt);  
        }
        else
        {
            if (*n_count > n_min)
            {   
                *dt *= dt_inc;

                if (*dt > dt_max)
                    *dt = dt_max;
            }   

            //*n_count ++;
            *n_count = (*n_count) + 1;

            printf ("\n");
            progress_tag ();
            printf("timestep:%f n_count:%d", *dt, *n_count);  
        }



    }

    /* Loop over ions */
    for (ion = 0; ion < ct.num_ions; ion++)
    {

        /* Get ion pointer */
        iptr = &ct.ions[ion];

        force[0] = iptr->force[fpt][0];
        force[1] = iptr->force[fpt][1];
        force[2] = iptr->force[fpt][2];

        if (ct.constrainforces)
        {
            force[0] += iptr->constraint.forcemask[0];
            force[1] += iptr->constraint.forcemask[1];
            force[2] += iptr->constraint.forcemask[2];
        }


        /* Use either actual ionic mass or equal mass for all atoms*/
        if (ct.relax_mass == 0)
            mass = ct.sp[iptr->species].atomic_mass * mu_me;
        else
            mass =  12.0 * mu_me;



        /* Magnitudes of f and v */
        magf = force[0] * force[0] +
               force[1] * force[1] +
               force[2] * force[2];


        /* Dot product of f and v */
        dotfv = force[0] * iptr->velocity[0] +
                force[1] * iptr->velocity[1] +
                force[2] * iptr->velocity[2];

        iptr->velocity[0] = *dt * force[0] / mass;
        iptr->velocity[1] = *dt * force[1] / mass;
        iptr->velocity[2] = *dt * force[2] / mass;


        if (dotfv >= 1.0e-12)
        {

            iptr->velocity[0] += dotfv * force[0] / magf;
            iptr->velocity[1] += dotfv * force[1] / magf;
            iptr->velocity[2] += dotfv * force[2] / magf;

        }
    }                           /* end for */

    /*Move ions by dt*velocity*/
    move_ions(*dt);

}                               /* end rmg_fastrelax */

/******/
