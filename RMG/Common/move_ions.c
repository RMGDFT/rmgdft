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

void move_ions (double dt)
{
    int ion, which = 0, count = 0;
    ION *iptr;
    double max_move = 0.0, avg_move = 0.0, rms_move = 0.0;
    double move_x, move_y, move_z, move_sq, move;


    /* Loop over ions */
    for (ion = 0; ion < ct.num_ions; ion++)
    {

        /* Get ion pointer */
        iptr = &ct.ions[ion];

        /*Save previous coordinates, needed for wavefunction extrapolation*/
        iptr->ocrds3[0] = iptr->ocrds2[0]; 
        iptr->ocrds3[1] = iptr->ocrds2[1]; 
        iptr->ocrds3[2] = iptr->ocrds2[2]; 

        iptr->ocrds2[0] = iptr->ocrds1[0]; 
        iptr->ocrds2[1] = iptr->ocrds1[1]; 
        iptr->ocrds2[2] = iptr->ocrds1[2]; 

        iptr->ocrds1[0] = iptr->crds[0]; 
        iptr->ocrds1[1] = iptr->crds[1]; 
        iptr->ocrds1[2] = iptr->crds[2]; 

        /* Move the ion */
        if (iptr->movable)
        {
            move_x = dt * iptr->velocity[0];
            move_y = dt * iptr->velocity[1];
            move_z = dt * iptr->velocity[2];


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


    }                           /* end for */

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

}                               /* end rmg_fastrelax */

/******/
