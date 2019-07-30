/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
/****f* QMD-MGDFT/write_pos.c *****
 * NAME
 *   Ab initio real space code with multigrid acceleration
 *   Quantum molecular dynamics package.
 *   Version: 2.1.5
 * COPYRIGHT
 *   Copyright (C) 1995  Emil Briggs
 *   Copyright (C) 1998  Emil Briggs, Charles Brabec, Mark Wensell, 
 *                       Dan Sullivan, Chris Rapcewicz, Jerzy Bernholc
 *   Copyright (C) 2001  Emil Briggs, Wenchang Lu,
 *                       Marco Buongiorno Nardelli,Charles Brabec, 
 *                       Mark Wensell,Dan Sullivan, Chris Rapcewicz,
 *                       Jerzy Bernholc
 * FUNCTION
 *   void write_pos(void)
 *   Writes out the postions of the ions and their displacements from their 
 *   initial postions. If the initial positions are negative values, they will
 *   add a lattice vector to make them positive							 
 * INPUTS
 *   nothing
 * OUTPUT
 *   print out atomic positions
 * PARENTS
 *   write_header.c
 * CHILDREN
 *   nothing
 * SOURCE
 */

#include <stdio.h>
#include <time.h>
#include <math.h>
#include "main.h"
#include "prototypes_on.h"

void write_pos(void)
{
    int ion;
    ION *iptr;


    printf("\n\n\n  IONIC POSITIONS AND DISPLACEMENTS:\n");


    printf("\nSpecies   X           Y           Z           dX          dY          dZ");

    for (ion = 0; ion < ct.num_ions; ion++)
    {

        iptr = &Atoms[ion];

        printf("\n  %d   %10.4f  %10.4f  %10.4f  %10.4f  %10.4f  %10.4f",
               iptr->species + 1,
               iptr->crds[0], iptr->crds[1], iptr->crds[2],
               iptr->crds[0] - iptr->icrds[0],
               iptr->crds[1] - iptr->icrds[1], iptr->crds[2] - iptr->icrds[2]);

    }                           /* end for */

    printf("\n");

}                               /* end write_pos */

/******/
