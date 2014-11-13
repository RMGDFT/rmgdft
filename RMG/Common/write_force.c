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
#include "main.h"







/* Writes out the positions of the ions and the current forces on them */
void write_force (void)
{
    int ion;
    ION *iptr;
    int num_movable = 0, maxfx_ion = 0, maxfy_ion = 0, maxfz_ion = 0, maxf_ion = 0;
    double avfx = 0.0, avfy = 0.0, avfz = 0.0, maxfx = 0.0, maxfy = 0.0, maxfz = 0.0;
    double sumx = 0.0, sumy = 0.0, sumz = 0.0;
    double avf = 0.0;
    double maxf = 0.0, max_all_f = 0.0;
    double f2;
    double *fp;

    printf ("\n\n\n  IONIC POSITIONS [a0] AND FORCES [Ha/a0]:\n\n");

//    if (verify ("atom_constraints", NULL))
//    {
//        printf ("  CONSTRAINTS ON FORCES HAVE BEEN IMPOSED:\n\n");
//    }
    printf
        ("@ION Ion Species           X           Y           Z          FX          FY          FZ movable\n");

    for (ion = 0; ion < ct.num_ions; ion++)
    {
        SPECIES *sp;

        iptr = &ct.ions[ion];

        fp = iptr->force[ct.fpt[0]];
        sp = &ct.sp[iptr->species];

        printf ("@ION %3d   %2s %2d  %10.7f  %10.7f  %10.7f  %10.7f  %10.7f  %10.7f %7d\n",
                ion + 1,
                sp->atomic_symbol,
                iptr->species + 1,
                iptr->crds[0], iptr->crds[1], iptr->crds[2], fp[0], fp[1], fp[2], iptr->movable);

        if (iptr->movable)
        {

            num_movable++;

            avfx += fabs (fp[0]);
            avfy += fabs (fp[1]);
            avfz += fabs (fp[2]);

            sumx += fp[0];
            sumy += fp[1];
            sumz += fp[2];

            if (fabs (fp[0]) > maxfx)
            {
                maxfx = fabs (fp[0]);
                maxfx_ion = ion;
            }
            
            if (fabs (fp[1]) > maxfy)
            {
                maxfy = fabs(fp[1]);
                maxfy_ion = ion;
            }
            
            if (fabs (fp[2]) > maxfz)
            {
                maxfz = fabs(fp[2]);
                maxfz_ion = ion;
            }
            
            
            

            f2 = fp[0] * fp[0] + fp[1] * fp[1] + fp[2] * fp[2];

            if (f2 > maxf)
            {
                maxf = f2;
                maxf_ion = ion;
            } 

            avf += f2;
        }
    }



    if (num_movable != 0)
    {
        avfx /= num_movable;
        avfy /= num_movable;
        avfz /= num_movable;

        maxf = sqrt (maxf);
        avf = sqrt (avf / num_movable);
        max_all_f = max (maxfx, maxfy);
        max_all_f = max (max_all_f, maxfz);


        printf ("\n");
        progress_tag ();
        printf (" mean FX      = %12.8f Ha/a0\n", avfx);
        progress_tag ();
        printf (" mean FY      = %12.8f Ha/a0\n", avfy);
        progress_tag ();
        printf (" mean FZ      = %12.8f Ha/a0\n", avfz);

        printf ("\n");
        progress_tag ();
        printf (" max FX       = %12.8f Ha/a0   (ion %d)\n", maxfx, maxfx_ion + 1);
        progress_tag ();
        printf (" max FY       = %12.8f Ha/a0   (ion %d)\n", maxfy, maxfy_ion + 1);
        progress_tag ();
        printf (" max FZ       = %12.8f Ha/a0   (ion %d)\n", maxfz, maxfz_ion + 1);
        progress_tag ();
        printf (" max F[x,y,z] = %12.8f Ha/a0\n", max_all_f);
        progress_tag ();
        printf (" max |F|      = %12.8f Ha/a0   (ion %d)\n", maxf, maxf_ion + 1);
        if (ct.forceflag == MD_FASTRLX)
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
        progress_tag ();
        printf (" sqrt < F^2 > = %12.8f Ha/a0\n", avf);
        printf ("\n");

    }

//    if ((pct.imgpe == 0) && (verify ("pdb_atoms", NULL)))
//        write_pdb ();

}                               /* end write_force */

/******/
