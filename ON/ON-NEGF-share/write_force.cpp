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
 *   forces are from Atoms[]...
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
#include "main.h"
#include "prototypes_on.h"







/* Writes out the postions of the ions and the current forces on them */
void write_force(void)
{
    int ion;
    ION *iptr;
    int num_movable = 0;
    double avfx = 0.0, avfy = 0.0, avfz = 0.0, maxfx = 0.0, maxfy = 0.0, maxfz = 0.0;
    double sumx = 0.0, sumy = 0.0, sumz = 0.0;
    double avf = 0.0;
    double maxf = 0.0, max_all_f = 0.0;
    double f2;

    printf("\n\n\n  IONIC POSITIONS [a0] AND FORCES [Ha/a0]:\n\n");

    printf("@ION Ion Species           X           Y           Z          FX          FY          FZ movable\n");

    for (ion = 0; ion < ct.num_ions; ion++)
    {
        SPECIES *sp;
        double *fp;
      
        iptr = &Atoms[ion];

        fp = iptr->force[ct.fpt[0]];
        sp = &Species[iptr->species];

        printf ( "@ION %3d   %2s %2d  %10.7f  %10.7f  %10.7f  %10.7f  %10.7f  %10.7f %7d\n",
                 ion + 1,
                 sp->atomic_symbol,
                 iptr->species + 1,
                 iptr->crds[0], iptr->crds[1], iptr->crds[2], 
                 fp[0], fp[1], fp[2],
                 iptr->movable );

        if (iptr->movable)
        {
            num_movable++;

            avfx += fabs( fp[0] );
            avfy += fabs( fp[1] );
            avfz += fabs( fp[2] );

            sumx += fp[0];
            sumy += fp[1];
            sumz += fp[2];

            maxfx = std::max( maxfx, fabs (fp[0]) );
            maxfy = std::max( maxfy, fabs (fp[1]) );
            maxfz = std::max( maxfz, fabs (fp[2]) );


            f2 = fp[0] * fp[0] + fp[1] * fp[1] + fp[2] * fp[2];
            maxf = std::max( maxf, f2 );
            avf += f2;
        }
    }



    if (num_movable != 0)
    {
        avfx = avfx / num_movable;
        avfy = avfy / num_movable;
        avfz = avfz / num_movable;

        maxf = sqrt( maxf );
        avf = sqrt( avf / num_movable );
        max_all_f = std::max( maxfx, maxfy );
        max_all_f = std::max( max_all_f, maxfz );


        printf ("\n");
        progress_tag ();
        printf (" mean FX      = %12.8f Ha/a0\n", avfx);
        progress_tag ();
        printf (" mean FY      = %12.8f Ha/a0\n", avfy);
        progress_tag ();
        printf (" mean FZ      = %12.8f Ha/a0\n", avfz);

        printf ("\n");
        progress_tag ();
        printf (" max FX       = %12.8f Ha/a0\n", maxfx);
        progress_tag ();
        printf (" max FY       = %12.8f Ha/a0\n", maxfy);
        progress_tag ();
        printf (" max FZ       = %12.8f Ha/a0\n", maxfz);
        progress_tag ();
        printf (" max F        = %12.8f Ha/a0\n", max_all_f);
        if ((ct.forceflag == MD_FASTRLX) )
        {
            progress_tag ();
            printf (" tolerance    = %12.8f Ha/a0\n", ct.thr_frc);
        }

        printf ("\n");
        progress_tag ();
        printf (" sum FX       = %12.8f Ha/a0\n", sumx);
        progress_tag ();
        printf (" sum FY       = %12.8f Ha/a0\n", sumy);
        progress_tag ();
        printf (" sum FZ       = %12.8f Ha/a0\n", sumz);
        progress_tag ();
        printf (" Average      = %12.8f Ha/a0\n", (fabs (sumx) + fabs (sumy) + fabs (sumz)) / 3.0);

        printf ("\n");
        progress_tag ();
        printf (" sqrt < F^2 > = %12.8f Ha/a0\n", avf);
        progress_tag ();
        printf (" max | F |    = %12.8f Ha/a0\n", maxf);
        printf ("\n");
    }

}                               /* end write_force */

/******/
