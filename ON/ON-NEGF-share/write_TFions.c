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

void write_tfions(void)
{
    int ion;
    TF_ION *iptr;


    printf("\n\n\n  Simplified water model:\n");
    printf("\n\n\n  TF IONIC POSITIONS AND PARAMETERS:\n");


    printf("\n          X           Y           Z           q          alpha       q0        alpha0");

    for (ion = 0; ion < ct.num_tfions; ion++)
    {

        iptr = &ct.tf_ions[ion];

        printf("\n   %10.4f  %10.4f  %10.4f  %10.4f  %10.4f  %10.4f  %10.4f",
               iptr->crds[0], iptr->crds[1], iptr->crds[2],
               iptr->q, iptr->alpha,
               iptr->q0, iptr->alpha0);

    }                           /* end for */

    printf("\n");

}                               /* end write_pos */

/******/
