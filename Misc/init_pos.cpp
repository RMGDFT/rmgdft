/****f* QMD-MGDFT/init_pos.c *****
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
 *   void init_pos()
 *   if positions are read in cartesian coordinates, get crystal coordinates
 *   if positions are read in crystal coordinates, get cartesian coordinates
 * INPUTS
 *   nothing
 * OUTPUT
 *   coordinates are stored in ct.ion 
 * PARENTS
 *   init.c
 * CHILDREN
 *   to_cartesian.c
 * SOURCE
 */

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "main.h"
#include "common_prototypes.h"

#define     SR2         1.414213562373
#define     SR3         1.732050807569

void init_pos ()
{
    int i, ir;
    ION *iptr;


    /* Scale the ion-coordinates by the lattice vectors (if required) */

    if (ct.crd_flag)
    {
        /* if positions in cartesian coordinates get crystal coordinates */
        for (i = 0; i < ct.num_ions; i++)
        {
            iptr = &Atoms[i];
            to_crystal(iptr->xtal, iptr->crds);
            to_cartesian (iptr->xtal, iptr->crds);
        }                       /* end for i */
    }
    else
    {
        /* If the positions are in crystal coordinates get cartesian coordinates */
        for (i = 0; i < ct.num_ions; i++)
        {
            iptr = &Atoms[i];
            to_cartesian (iptr->xtal, iptr->crds);
        }                       /* end for i */
    }                           /* if(ct.crd_flag) */

}                               /* end init_pos */

/******/
