/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

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
    int i, ir, ibrav;
    ION *iptr;

    ibrav = get_ibrav_type();

    /* Scale the ion-coordinates by the lattice vectors (if required) */

    if (ct.crd_flag)
    {

        /* if positions in cartesian coordinates
           get crystal coordinates */

        /* add part to bring into supercell (just in case) */

        for (i = 0; i < ct.num_ions; i++)
        {

            iptr = &ct.ions[i];

            for (ir = 0; ir < 3; ir++)
            {
                iptr->xtal[ir] = iptr->crds[0] * get_b0(ir)
                    + iptr->crds[1] * get_b1(ir) + iptr->crds[2] * get_b2(ir);
            }                   /* end for ir */

            if (ibrav == HEXAGONAL)
            {

                iptr->xtal[0] = (iptr->crds[0] - iptr->crds[1] / SR3) / get_celldm(0);
                iptr->xtal[1] = iptr->crds[1] / (SR3 / 2.0) / get_celldm(0);
                iptr->xtal[2] = iptr->crds[2] * get_b2(2);


            }                   /* end if */

            if (iptr->xtal[0] < 0.0)
                iptr->xtal[0] += 1.0;
            if (iptr->xtal[1] < 0.0)
                iptr->xtal[1] += 1.0;
            if (iptr->xtal[2] < 0.0)
                iptr->xtal[2] += 1.0;
            if (iptr->xtal[0] > 1.0)
                iptr->xtal[0] -= 1.0;
            if (iptr->xtal[1] > 1.0)
                iptr->xtal[1] -= 1.0;
            if (iptr->xtal[2] > 1.0)
                iptr->xtal[2] -= 1.0;
            to_cartesian (iptr->xtal, iptr->crds);

        }                       /* end for i */
    }
    else
    {
        /* If the positions are in crystal coordinates
           get cartesian coordinates */

        /* add part to bring into supercell (just in case) */
        for (i = 0; i < ct.num_ions; i++)
        {
            iptr = &ct.ions[i];
            for (ir = 0; ir < 3; ir++)
            {
                iptr->crds[ir] = iptr->xtal[0] * get_a0(ir)
                    + iptr->xtal[1] * get_a1(ir) + iptr->xtal[2] * get_a2(ir);
            }                   /* end for ir */

        }                       /* end for i */
    }                           /* if(ct.crd_flag) */

}                               /* end init_pos */

/******/
