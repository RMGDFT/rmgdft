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
    int which = 0, count = 0;
    double max_move = 0.0, avg_move = 0.0, rms_move = 0.0;
    double move_sq, move;


    /* Loop over ions */
    for (size_t ion = 0, i_end = Atoms.size(); ion < i_end; ++ion)
    {

        /* Get ion pointer */
        ION &Atom = Atoms[ion];

        /*Save previous coordinates, needed for wavefunction extrapolation*/
        Atom.RotateCoordinates();

        /* Move the ion */
        int movable = Atom.movable[0] + Atom.movable[1] + Atom.movable[2];
        if (movable)
        {
            double move_x = dt * Atom.velocity[0] * Atom.movable[0];
            double move_y = dt * Atom.velocity[1] * Atom.movable[1];
            double move_z = dt * Atom.velocity[2] * Atom.movable[2];


            /*Update coordinates*/
            Atom.crds[0] += move_x;
            Atom.crds[1] += move_y;
            Atom.crds[2] += move_z;

            /* enforce periodic boundary conditions on the ions */
            to_crystal (Atom.xtal, Atom.crds);
            to_cartesian (Atom.xtal, Atom.crds);

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
    if(ct.verbose)
    {
        printf ("\n\n");
        progress_tag ();
        printf ("Max displacement: %8.5f a0  (ion %d)", max_move, which + 1);
        printf ("\n");
        progress_tag ();
        printf ("Avg displacement: %8.5f a0", avg_move/count);
        printf ("\n");
        progress_tag ();
        printf ("RMS displacement: %8.5f a0", sqrt(rms_move/count));
        printf ("\n");
    }
}                               /* end move_ions */

/******/
