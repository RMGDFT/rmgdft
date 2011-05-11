/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
/****f* QMD-MGDFT/iiforce.c *****
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
 *   void iiforce(void)
 *   Calculates ion-ion component of the forces.
 *   Uses the minimum image convention for periodic cells.
 * INPUTS
 *   nothing
 * OUTPUT
 *   The forces for each ion are stored in the main
 *   CONTROL structure ct. The values calculated here are
 *   added to the values stored in that structure.
 * PARENTS
 *   force.c
 * CHILDREN
 *   nothing
 * SOURCE
 */



#include "md.h"
#include <float.h>
#include <math.h>


void iiforce(void)
{

    int i, j;
    REAL Zi, Zj, rci, rcj, t1, t2, s1, s2, s3, n1, r;
    REAL xtal_r[3], crd_r[3];
    ION *iptr1, *iptr2;
    REAL tx, ty, tz;

    REAL time1, time2;
    time1 = my_crtc();

    n1 = TWO / sqrt(PI);

    /* Loop over ions and get the ion-ion component of the forces */
    for (i = 0; i < ct.num_ions; i++)
    {

        tx = ZERO;
        ty = ZERO;
        tz = ZERO;

        /* Get ion pointer */
        iptr1 = &ct.ions[i];


        Zi = ct.sp[iptr1->species].zvalence;
        rci = ct.sp[iptr1->species].rc;


        /* Sum contributions from the rest of the ions. */
        for (j = 0; j < ct.num_ions; j++)
        {


            if (j != i)
            {

                iptr2 = &ct.ions[j];

                Zj = ct.sp[iptr2->species].zvalence;
                rcj = ct.sp[iptr2->species].rc;


                t1 = rci * rci + rcj * rcj;
                t2 = sqrt(t1);


                /* Minimum image convention for r */
                r = minimage(iptr1, iptr2, xtal_r);

                s1 = Zi * Zj / (r * r);

                s2 = erfc(r / t2) / r;

                s3 = n1 * exp(-r * r / t1) / t2;

                to_cartesian(xtal_r, crd_r);

                iptr1->force[ct.fpt[0]][0] += crd_r[0] * s1 * (s2 + s3);
                iptr1->force[ct.fpt[0]][1] += crd_r[1] * s1 * (s2 + s3);
                iptr1->force[ct.fpt[0]][2] += crd_r[2] * s1 * (s2 + s3);

                tx += crd_r[0] * s1 * (s2 + s3);
                ty += crd_r[1] * s1 * (s2 + s3);
                tz += crd_r[2] * s1 * (s2 + s3);


            }                   /* end if */

        }                       /* end for */


    }                           /* end for */

    time2 = my_crtc();
    rmg_timings(IIFORCE_TIME, time2 - time1);

}                               /* end iiforce */


/******/
