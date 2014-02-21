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
#include "init_var.h"
void read_atoms (void)
{
    int st1,ist, ion, args;
    char tbuf[MAX_PATH];
    char species[32];

    /*Count number of ions in input file */
    get_data ("atoms", &ion, INIT | LIST, NULL);

    if (verify ("pdb_atoms", NULL))
    {
        if (ion != ct.num_ions)
            error_handler ("Mismatch between PDB and regular atons");
    }
    else
    {
        ct.num_ions = ion;
    }



    //    if (pct.gridpe == 0)
    //       printf ("\n Number of ions is %d", ct.num_ions);

    /*Allocate memory for ions */
    my_calloc (ct.ions, ct.num_ions, ION);

    /* read and count coordinates for each ion */
    ion = 0;
    ist = 0;

    while (get_data ("atoms", tbuf, ITEM | STR, NULL))
    {


      int args = 0;
        args = sscanf (tbuf, "%s %lf %lf %lf %d %d %d",
                species,
                &ct.ions[ion].crds[0], &ct.ions[ion].crds[1], &ct.ions[ion].crds[2],
                &ct.ions[ion].movable, &ct.ions[ion].frozen, &ct.ions[ion].n_loc_states
                );

        if (args < 7)
        {
            printf ("Error reading ion %d args %d\n", ion, args);
            error_handler ("Not enough arguments in atoms list!");
        }

        /*Find the species of the atom */
        /* In case of pdb input verify species */
        if (verify ("pdb_atoms", NULL))
        {
            if (ct.ions[ion].species != assign_species (&ct, species))
                error_handler ("Mismatch between species in RMG and PDB input");
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
        else if (verify ("length_units", "Angstrom"))
        {
            ct.ions[ion].crds[0] *= A_a0;
            ct.ions[ion].crds[1] *= A_a0;
            ct.ions[ion].crds[2] *= A_a0;
        }



        for(st1 = ist; st1 < ist+ct.ions[ion].n_loc_states; st1++)
        {
            state_to_ion[st1] = ion;
            states[st1].atomic_orbital_index = st1-ist;
        }

        ist += ct.ions[ion].n_loc_states;


        ion++;

    }                           /* end while (get_data( "atoms", tbuf, ITEM | STR, NULL) != NULL) */

    if (ion != ct.num_ions)
    {
        printf ("ion = %d != %d = ct.num_ions", ion, ct.num_ions);
        error_handler ("Mismatch in number of ions");
    }

}
