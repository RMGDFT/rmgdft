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
#include <string.h>
#include "main.h"

/*Functions here handle generating xbs movie for dynamic simulations*/
static void init_xbsmovie (FILE * movie);


/*This writes .bs file (initial positions, bond info etc.) and returns pointer to 
 * the .mv file where further frames will go*/
FILE *open_xbs_movie (char *filename)
{
    FILE *xbsfp, *xbsfp1;
    char filename_bs[80];
    char filename_mv[80];


    /*file name where initial frame will go */
    sprintf (filename_bs, "%s%s", filename, ".bs");
    /*file name where other frames will go */
    sprintf (filename_mv, "%s%s", filename, ".mv");

    my_fopen (xbsfp, filename_bs, "w");
    my_fopen (xbsfp1, filename_mv, "w");

    /*This check is taken from rmmovie opening, I am not exactly sure what it does */
    if (setvbuf (xbsfp, (char *) NULL, _IOFBF, 4096 * 16) != 0)
        printf ("\n Warning: cant allocate XBS movie io buffer size\n");
    if (setvbuf (xbsfp1, (char *) NULL, _IOFBF, 4096 * 16) != 0)
        printf ("\n Warning: cant allocate XBS movie io buffer size\n");

    /*output initial info into xbs file */
    init_xbsmovie (xbsfp);

    /*This file can be closed, all other information should go into the .mv file */
    fclose (xbsfp);

    return (xbsfp1);
}




static void init_xbsmovie (FILE * movie)
{

    int ion, species, species2;
    ION *iptr;
    SPECIES *sp, *sp2;
    int max_atomic_number, min_atomic_number;
    double radius, color;

    max_atomic_number = 0;
    min_atomic_number = 1000;

    for (ion = 0; ion < ct.num_ions; ion++)
    {

        /* Get ion pointer */
        iptr = &ct.ions[ion];
        sp = &ct.sp[iptr->species];

        fprintf (movie, "atom %2s	%f %f %f\n",
                 sp->atomic_symbol, iptr->crds[0], iptr->crds[1], iptr->crds[2]);

        /*determine max and min atomic number */
        if (sp->atomic_number > max_atomic_number)
            max_atomic_number = sp->atomic_number;
        if (sp->atomic_number < min_atomic_number)
            min_atomic_number = sp->atomic_number;

    }                           /*end for(ion = 0;ion < ct.num_ions;ion++) */

    fprintf (movie, "\n\n");

    /*Loop over all species */
    for (species = 0; species < ct.num_species; species++)
    {

        sp = &ct.sp[species];

        /*Assign radii to species (from 0.5 to 1 according to atomic number) */
        radius =
            0.5 + 0.5 * (sp->atomic_number - min_atomic_number) / (max_atomic_number -
                                                                   min_atomic_number);

        /*Assign color for species */
        color = 1.0 - radius;

        fprintf (movie, "spec      %2s      %3f   %3f %3f  %3f \n", sp->atomic_symbol, radius,
                 color, color, color);

    }                           /*end for(species=0; species < ct.num_species; species++) */

    fprintf (movie, "\n\n");

    /*Now bonds specification for all possible combinations */
    for (species = 0; species < ct.num_species; species++)
    {
        sp = &ct.sp[species];
        for (species2 = 0; species2 < ct.num_species; species2++)
        {
            sp2 = &ct.sp[species2];


            /*Format should be: name1 name2 min_length max_length radius color */
            fprintf (movie, "bonds   %2s    %2s  0.00  3.5  0.10  0.95\n",
                     sp->atomic_symbol, sp2->atomic_symbol);

        }
    }


    fprintf (movie, "\n\n");

    fprintf (movie, "inc 10\n");
    fprintf (movie, "step 1\n");

}

void xbsmovie (FILE * movie)
{
    int ion;
    ION *iptr;

    fprintf (movie, "frame 	This is a frame\n");


    for (ion = 0; ion < ct.num_ions; ion++)
    {

        /* Get ion pointer */
        iptr = &ct.ions[ion];

        fprintf (movie, "%f %f %f ", iptr->crds[0], iptr->crds[1], iptr->crds[2]);
    }

    fprintf (movie, "\n\n");
}
