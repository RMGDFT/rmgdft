/************************** SVN Revision Information **************************
 **    $Id: iiforce.c 3533 2016-05-02 19:10:01Z luw $    **
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



#include <float.h>
#include <math.h>
#include <complex>
#include "const.h"
#include "params.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "rmg_error.h"
#include "State.h"
#include "Kpoint.h"
#include "Functional.h"
#include "GlobalSums.h"
#include "transition.h"
#include "RmgSumAll.h"


void IIforce (double *force)
{

    int i, j;
    double Zi, Zj, rci, rcj, t1, t2, s1, s2, s3, n1, r;
    double xtal_r[3], crd_r[3];
    ION *iptr1, *iptr2;
    double x, y, z;
    double sigma;

    sigma = 0.0;
    for (i = 0; i<ct.num_species; i++) 
        sigma = std::max(sigma, ct.sp[i].rc);


    //  the ion-ion term in real space
    //  check how many unit cell to make erfc(nL/sqrt(2.0)*sigma) < tol ==1.0e-8
    //  erfc(4.06) = 9.37e-9 so we set the largest r/(sqrt(2)*sigma) = 4.06


    double rcutoff = 4.06;
    int num_cell_x = int(rcutoff * sqrt(2.0) * sigma/Rmg_L.get_xside()) + 1;
    int num_cell_y = int(rcutoff * sqrt(2.0) * sigma/Rmg_L.get_yside()) + 1;
    int num_cell_z = int(rcutoff * sqrt(2.0) * sigma/Rmg_L.get_zside()) + 1;


    n1 = 2.0 / sqrt (PI);

    /* Loop over ions and get the ion-ion component of the forces */
    for (i = pct.gridpe; i < ct.num_ions; i+=pct.grid_npes)
    {

        /* Get ion pointer */
        iptr1 = &ct.ions[i];


        Zi = ct.sp[iptr1->species].zvalence;
        rci = ct.sp[iptr1->species].rc;


        /* Sum contributions from the rest of the ions. */
        for (j = 0; j < ct.num_ions; j++)
        {

            iptr2 = &ct.ions[j];

            Zj = ct.sp[iptr2->species].zvalence;
            rcj = ct.sp[iptr2->species].rc;


            t1 = rci * rci + rcj * rcj;
            t2 = sqrt (t1);

            for(int ix = -num_cell_x; ix<= num_cell_x; ix++)
                for(int iy = -num_cell_y; iy<= num_cell_y; iy++)
                    for(int iz = -num_cell_z; iz<= num_cell_z; iz++)
                    {
                        x = iptr1->crds[0] - iptr2->crds[0] + ix * Rmg_L.a0[0] + iy * Rmg_L.a1[0] + iz * Rmg_L.a2[0];
                        y = iptr1->crds[1] - iptr2->crds[1] + ix * Rmg_L.a0[1] + iy * Rmg_L.a1[1] + iz * Rmg_L.a2[1];
                        z = iptr1->crds[2] - iptr2->crds[2] + ix * Rmg_L.a0[2] + iy * Rmg_L.a1[2] + iz * Rmg_L.a2[2];
                        r = sqrt(x*x + y*y + z*z);

                        // r= 0 means two atoms are the same one.
                        if(r > 1.0e-5) {


                            s1 = Zi * Zj / (r * r);

                            s2 = erfc (r / t2) / r;

                            s3 = n1 * exp (-r * r / t1) / t2;

                            to_cartesian (xtal_r, crd_r);

                            force[i*3 + 0] += x * s1 * (s2 + s3);
                            force[i*3 + 1] += y * s1 * (s2 + s3);
                            force[i*3 + 2] += z * s1 * (s2 + s3);


                        }                   /* end if */

                    }                       /* end for */

        }                           /* end for */

    }
}                               /* end iiforce */


