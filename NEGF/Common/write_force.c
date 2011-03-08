/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
/****f* QMD-MGDFT/write_force.c *****
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
 *   void write_force(void)
 *   Writes out the postions of the ions and the current forces on them 
 * INPUTS
 *   forces are from ct.ions[]...
 * OUTPUT
 *   print out atomic coordinate and forces
 * PARENTS
 *   cdfastrlx.c fastrlx.c moldyn.c quench.c
 * CHILDREN
 *   nothing
 * SOURCE
 */


#include <stdio.h>
#include <time.h>
#include <math.h>
#include "md.h"



/* Writes out the postions of the ions and the current forces on them */
void write_force (void)
{
    int ion;
    ION *iptr;
    int num_movable = 0;
    REAL avfx = 0.0, avfy = 0.0, avfz = 0.0, maxfx = 0.0, maxfy = 0.0, maxfz = 0.0;


    printf("\n\n\n  IONIC POSITIONS [a0] AND FORCES [Ha/a0]:\n\n");

    printf
        ("@ION Ion Species    X           Y           Z          FX          FY          FZ          movable\n");

    for (ion = 0; ion < ct.num_ions; ion++)
    {


        iptr = &ct.ions[ion];

        printf ("@ION %3d %3d    %10.7f  %10.7f  %10.7f  %10.7f  %10.7f  %10.7f  %5d\n",
                ion + 1,
                iptr->species + 1,
                iptr->crds[0], iptr->crds[1], iptr->crds[2],
                iptr->force[ct.fpt[0]][0], iptr->force[ct.fpt[0]][1],
                iptr->force[ct.fpt[0]][2], iptr->movable);

        if (iptr->movable)
        {
            num_movable++;
            avfx += fabs (iptr->force[ct.fpt[0]][0]);
            avfy += fabs (iptr->force[ct.fpt[0]][1]);
            avfz += fabs (iptr->force[ct.fpt[0]][2]);

            maxfx = max (maxfx, fabs (iptr->force[ct.fpt[0]][0]));
            maxfy = max (maxfy, fabs (iptr->force[ct.fpt[0]][1]));
            maxfz = max (maxfz, fabs (iptr->force[ct.fpt[0]][2]));
        }                       /* if(iptr->movable) */
    }                           /* end for */

    if (num_movable != 0)
    {
        avfx = avfx / num_movable;
        avfy = avfy / num_movable;
        avfz = avfz / num_movable;
    }                           /* if (num_movable!=0) */
    printf ("\n mean FX = %12.8f", avfx);
    printf ("\n mean FY = %12.8f", avfy);
    printf ("\n mean FZ = %12.8f", avfz);

    printf ("\n max FX = %f", maxfx);
    printf ("\n max FY = %f", maxfy);
    printf ("\n max FZ = %f\n", maxfz);


}                               /* end write_force */

/******/
