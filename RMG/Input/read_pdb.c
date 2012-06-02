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

#include <string.h>
#include <ctype.h>
#include "main.h"


void read_parse_pdb_line (ION * iptr, int ion_index)
{
    char temp_str[10];
    char tbuf[256];
    char *temp_ptr;


    get_data ("pdb_atoms", &tbuf, ITEM | STR, NULL);



#if 1
    /*A line should have at least 78 characters, element symbol is 77-78, right justified */
    if (strlen (tbuf) < 78)
    {
        printf
            ("\n Line %d  in PDB input only has %d characters, but minimum of 78 is required",
             ion_index, (int) strlen (tbuf));
        error_handler ("Not enough characters in PDB line");
    }


    /* 1 -  6  Record name */
    my_strncpy (iptr->pdb.record_name, &tbuf[0], 6);


    /* 7 - 11 Atom serial number */
    my_strncpy (temp_str, &tbuf[6], 5);
    iptr->pdb.serial_num = strtol (temp_str, NULL, 10);


    /*13 - 16  Atom name */
    my_strncpy (iptr->pdb.name, &tbuf[12], 4);


    /* 17 Alternate location indicator. */
    my_strncpy (iptr->pdb.altLoc, &tbuf[16], 1);


    /* 18 - 20 Residue name */
    my_strncpy (iptr->pdb.resName, &tbuf[17], 3);


    /* 22 Chain identifier */
    my_strncpy (iptr->pdb.chainID, &tbuf[21], 1);


    /* 23 - 26 Residue sequence number */
    my_strncpy (temp_str, &tbuf[22], 4);
    iptr->pdb.resSeq = strtol (temp_str, NULL, 10);


    /* 27 Code for insertion of residues */
    my_strncpy (iptr->pdb.iCode, &tbuf[26], 1);



    /*Read crds */
    my_strncpy (temp_str, &tbuf[30], 8);
    iptr->crds[0] = A_a0 * strtod (temp_str, &temp_ptr);

    my_strncpy (temp_str, &tbuf[38], 8);
    iptr->crds[1] = A_a0 * strtod (temp_str, &temp_ptr);

    my_strncpy (temp_str, &tbuf[46], 8);
    iptr->crds[2] = A_a0 * strtod (temp_str, &temp_ptr);


    /*In PDB file, movable flags cannot be specified, therefore setup
     * everything movable by default*/
    iptr->movable = 1;



    /* 55 - 60 Occupancy */
    my_strncpy (temp_str, &tbuf[54], 6);
    iptr->pdb.occupancy = strtod (temp_str, &temp_ptr);


    /* 61 - 66 Temperature factor */
    my_strncpy (temp_str, &tbuf[60], 6);
    iptr->pdb.tempFactor = strtod (temp_str, &temp_ptr);

    /* 73 - 76  Element symbol, right-justified. */
    //my_strncpy (iptr->pdb.segID, &tbuf[72], 4);


    /* 77 - 78  Element symbol, right-justified. */
    my_strncpy (iptr->pdb.element, &tbuf[76], 2);

    /*Species symbol only, stripped of possible leading space */
    if (tbuf[76] == ' ')
        my_strncpy (temp_str, &tbuf[77], 1);
    else
        my_strncpy (temp_str, &tbuf[76], 2);

    /*Assign species number according to element symbol */
    iptr->species = assign_species (&ct, temp_str);

    Dprintf("found pdb species %s, %d: from ion %d", temp_str, iptr->species, ion_index);


    /*79 - 80  Charge on the atom. */
    if (strlen (tbuf) >= 79)
        my_strncpy (iptr->pdb.charge, &tbuf[78], 2);
    else
        my_strncpy (iptr->pdb.charge, "  ", 2);






#endif


}





void read_pdb (void)
{
    int ion;


    /*Absolute coordinates have to be used with PDB input */
    if (!verify ("atomic_coordinate_type", "Absolute"))
        error_handler ("atomic_coordinate_type has to be set to absolute when PDB input is used");



    /* read and count coordinates for each ion */
    for (ion = 0; ion < ct.num_ions; ion++)
    {
        read_parse_pdb_line (&ct.ions[ion], ion);
    }

}
