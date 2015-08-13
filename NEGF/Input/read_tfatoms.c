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
#include "LCR.h"

void read_tfatoms (void)
{
    int st1,ist, ion, args;
    char tbuf[MAX_PATH];
    char species[32];

    ct.num_tfions = 0;
        
    if(pct.gridpe == 0) 
	printf("\n  *Starting read_tfatoms \n");
    
    /*Count number of ions in input file */
    get_data ("tf_atoms", &ion, INIT | LIST, NULL);
    
    if(pct.gridpe == 0) 
	printf("\n  Number TF ions determined as %d \n, ion");

    if (verify ("tf_atoms", NULL))
    {
	ct.num_tfions = ion;
	//    if (pct.gridpe == 0)
	//       printf ("\n Number of ions is %d", ct.num_ions);

	my_malloc (ct.tf_ions, ct.num_tfions, TF_ION);


	/* read and count coordinates for each ion */
	ion = 0;
	ist = 0;

	while (get_data ("tf_atoms", tbuf, ITEM | STR, NULL))
	{

	    int args = 0;
	    args = sscanf (tbuf, " %lf %lf %lf %lf %lf %lf %lf",
		    &ct.tfPions[ion].crds[0], &ct.tf_ions[ion].crds[1], &ct.tf_ions[ion].crds[2],
		    &ct.tf_ions[ion].q, &ct.tf_ions[ion].alpha, &ct.tf_ions[ion].q0, &ct.tf_ions[ion].alpha0
		    );

	    if (args < 7)
	    {
		printf ("Error reading tf_ion %d: %d args provided but 7 expected\n", ion, args);
		error_handler ("Not enough arguments in tf_atoms list!");
	    }

        /* Handle lattice and angstrom units */
        if (verify ("atomic_coordinate_type", "Cell Relative"))
        {
	    printf ("!!! WARNING: Cell Relative is specifies, but TF ions are not being converted, 
		    because they are expected to be in absolute coordinates");
        }
	else if (verify ("length_units", "Angstrom"))
        {
            ct.tf_ions[ion].crds[0] *= A_a0;
            ct.tf_ions[ion].crds[1] *= A_a0;
            ct.tf_ions[ion].crds[2] *= A_a0;
        }


        ion++;

	}                           /* end while (get_data( "atoms", tbuf, ITEM | STR, NULL) != NULL) */

    }
    
    if(pct.gridpe == 0) 
	printf("\n  Read %d TF ions\n, ion");

}
