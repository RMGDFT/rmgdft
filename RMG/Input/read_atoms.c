/************************** SVN Revision Information **************************
 **    $Id$    **
 ******************************************************************************/

/****f* QMD-MGDFT/read_control.c *****
 * NAME
 *   Ab initio real space code with multigrid acceleration
 *   Quantum molecular dynamics package.
 * COPYRIGHT
 *   Copyright (C) 2005  Frisco Rose
 *                       Jerzy Bernholc
 *   
 * FUNCTION
 *	void read_control(CONTROL *c)
 * 		read all information from control input file
 *    
 * INPUTS
 *   main control structure
 * OUTPUT
 *   variables in structure CONTROL c are updated
 *   in most of other file, the name is ct.... 
 *   see main.h for structure CONTROL
 * PARENTS
 *   main.c
 * CHILDREN
 * 
 * SOURCE */

#include "main.h"

void read_atoms (void)
{
    int ion, args;
    char tbuf[MAX_PATH];
    char species[32];


    /* read and count coordinates for each ion */
    ion = 0;

    while (get_data ("atoms", tbuf, ITEM | STR, NULL))
    {
        args =
            sscanf (tbuf, "%s %lf %lf %lf %d", species, &ct.ions[ion].crds[0],
                    &ct.ions[ion].crds[1], &ct.ions[ion].crds[2], &ct.ions[ion].movable);

        if (args == 4)
        {
            get_data ("atoms", NULL, INIT | INFO, "#DEFAULT 1 on movable");
            ct.ions[ion].movable = 1;
        }
        else if (args < 4)
        {
            error_handler ("Error reading ion %d, not enough arguments in \"%s\"!", ion, tbuf );
        }

        /*Find the species of the atom */
        /* In case of pdb input verify species */
        if (verify ("pdb_atoms", NULL))
        {
            if (ct.ions[ion].species != assign_species (&ct, species))
                error_handler ("Mismatch between species in RMG and PDB input,"\
                               " expected %d - found %d for ion %d",\
                               ct.ions[ion].species, assign_species (&ct, species), ion);
        }
        else
        {
            ct.ions[ion].species = assign_species (&ct, species);
        }


        /* Handle lattice and angstrom units */
        if (verify ("atomic_coordinate_type", "Cell Relative"))
        {
            ct.ions[ion].xtal[0] = ct.ions[ion].crds[0];
            ct.ions[ion].xtal[1] = ct.ions[ion].crds[1];
            ct.ions[ion].xtal[2] = ct.ions[ion].crds[2];
        }
        else if (verify ("crds_units", "Angstrom"))
        {
            ct.ions[ion].crds[0] *= A_a0;
            ct.ions[ion].crds[1] *= A_a0;
            ct.ions[ion].crds[2] *= A_a0;
        }

        ion++;

    }                           /* end while (get_data( "atoms", tbuf, ITEM | STR, NULL) != NULL) */

    if (ion != ct.num_ions)
    {
        printf ("ion = %d != %d = ct.num_ions", ion, ct.num_ions);
        error_handler ("Mismatch in number of ions");
    }

}
