/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

/****f* QMD-MGDFT/force.c *****
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
 *   void force(REAL *rho, REAL *rhoc, REAL *vh, REAL *vxc, STATE *states)
 *   Driver routine to calculate ionic forces.
 * INPUTS
 *   rho: total charge density
 *   rhoc: compensating charge density
 *   vh:   Hartree potential
 *   vxc:  exchange-correlation potential
 *   state: points to orbital structure which include eigenvalues, 
 *          wave functions and so on. (see main.h)
 * OUTPUT
 *   Resulted forces are stored in structure ct.xxxxx 
 * PARENTS
 *   cdfastrlx.c fastrlx.c moldyn.c quench.c
 * CHILDREN
 *   iiforce.c nlforce.c lforce.c nlccforce.c
 * SOURCE
 */


#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "main.h"


/*Set this to 1 to have forces written out part by part*/
/* If you want this , you should also make sure that VERBOSE flag is enabled in
 * nlforce1.c*/
#define VERBOSE 0



void force (REAL * rho, REAL * rho_oppo, REAL * rhoc, REAL * vh, REAL * vxc, REAL * vnuc, STATE * states)
{
    int ion, st, kpt, idx, nspin = (ct.spin_flag + 1);
    REAL *vtot, *rho_tot, meanres;
    STATE *sp;
    REAL time1, time2, time3;
#if VERBOSE
    REAL *old_force;
    REAL sumx, sumy, sumz;
#endif

    int Zi;
    time3 = my_crtc ();

#if VERBOSE
    my_malloc (old_force, 3 * ct.num_ions, REAL);

    for (ion = 0; ion < ct.num_ions; ion++)
    {
        old_force[3 * ion] = ZERO;
        old_force[3 * ion + 1] = ZERO;
        old_force[3 * ion + 2] = ZERO;
    }
#endif

    my_malloc (vtot, FP0_BASIS, REAL);
    for (idx = 0; idx < FP0_BASIS; idx++)
        vtot[idx] = vxc[idx] + vh[idx] + vnuc[idx];

    /* Compute residual information */
    meanres = 0.0;
    for (kpt = 0; kpt < ct.num_kpts; kpt++)
    {
        sp = ct.kp[kpt].kstate;
        for (st = 0; st < ct.num_states; st++)
        {
            meanres += sp->res;
            sp++;
        }
    }

    if (ct.spin_flag)
	    meanres = real_sum_all (meanres, pct.spin_comm);

    ct.meanres = meanres / ((REAL) (ct.num_kpts * ct.num_states * nspin));

    /* Print out residual information */

    printf ("\n@@Force Mean Occ Subspace Res = %15.8e", ct.meanres);
    printf ("\n@@Force Max Occ Subspace Res   = %15.8e", ct.maxres);
    printf ("\n@@Force Min Occ Subspace Res   = %15.8e", ct.minres);



    /* Zero out forces */
    for (ion = 0; ion < ct.num_ions; ion++)
    {

      Zi = ct.sp[ct.ions[ion].species].zvalence;


        ct.ions[ion].force[ct.fpt[0]][0] = ct.e_field * ct.x_field_0 * Zi;
        ct.ions[ion].force[ct.fpt[0]][1] = ct.e_field * ct.y_field_0 * Zi;
        ct.ions[ion].force[ct.fpt[0]][2] = ct.e_field * ct.z_field_0 * Zi;

    }

    /* Get the ion-ion component and store. */
    iiforce ();

#if VERBOSE
    if (pct.imgpe == 0)
    {
        printf ("\n\n Ion-Ion force:");

        sumx = ZERO;
        sumy = ZERO;
        sumz = ZERO;
        for (ion = 0; ion < ct.num_ions; ion++)
        {
            printf ("\n Ion %d Force  %10.7f  %10.7f  %10.7f",
                    ion, ct.ions[ion].force[ct.fpt[0]][0] - old_force[3 * ion],
                    ct.ions[ion].force[ct.fpt[0]][1] - old_force[3 * ion + 1],
                    ct.ions[ion].force[ct.fpt[0]][2] - old_force[3 * ion + 2]);

            sumx += ct.ions[ion].force[ct.fpt[0]][0] - old_force[3 * ion];
            sumy += ct.ions[ion].force[ct.fpt[0]][1] - old_force[3 * ion + 1];
            sumz += ct.ions[ion].force[ct.fpt[0]][2] - old_force[3 * ion + 2];

            old_force[3 * ion] = ct.ions[ion].force[ct.fpt[0]][0];
            old_force[3 * ion + 1] = ct.ions[ion].force[ct.fpt[0]][1];
            old_force[3 * ion + 2] = ct.ions[ion].force[ct.fpt[0]][2];

        }
        printf ("\n II sums in x, y and z directions: %e %e %e", sumx, sumy, sumz);
    }
#endif



    /* Add in the local */
    if (ct.spin_flag)
    {
    	my_malloc (rho_tot, FP0_BASIS, REAL);
	for (idx = 0; idx < FP0_BASIS; idx++)
		rho_tot[idx] = rho[idx] + rho_oppo[idx];
	lforce(rho_tot, vh);
	my_free (rho_tot);
    }
    else
    	lforce (rho, vh);

#if VERBOSE
    if (pct.imgpe == 0)
    {
        printf ("\n\n Local force:");

        sumx = ZERO;
        sumy = ZERO;
        sumz = ZERO;
        for (ion = 0; ion < ct.num_ions; ion++)
        {
            printf ("\n Ion %d Force  %10.7f  %10.7f  %10.7f",
                    ion, ct.ions[ion].force[ct.fpt[0]][0] - old_force[3 * ion],
                    ct.ions[ion].force[ct.fpt[0]][1] - old_force[3 * ion + 1],
                    ct.ions[ion].force[ct.fpt[0]][2] - old_force[3 * ion + 2]);

            sumx += ct.ions[ion].force[ct.fpt[0]][0] - old_force[3 * ion];
            sumy += ct.ions[ion].force[ct.fpt[0]][1] - old_force[3 * ion + 1];
            sumz += ct.ions[ion].force[ct.fpt[0]][2] - old_force[3 * ion + 2];

            old_force[3 * ion] = ct.ions[ion].force[ct.fpt[0]][0];
            old_force[3 * ion + 1] = ct.ions[ion].force[ct.fpt[0]][1];
            old_force[3 * ion + 2] = ct.ions[ion].force[ct.fpt[0]][2];
        }
        printf ("\n Local sums in x, y and z directions: %e %e %e", sumx, sumy, sumz);
    }
#endif

    time1 = my_crtc ();

    /* Add in the non-local stuff */
    nlforce1 (vtot);

    time2 = my_crtc ();
    rmg_timings (NLFORCE_TIME, (time2 - time1));

#if VERBOSE
    if (pct.imgpe == 0)
    {
        printf ("\n\n Non-Local force:");

        sumx = ZERO;
        sumy = ZERO;
        sumz = ZERO;
        for (ion = 0; ion < ct.num_ions; ion++)
        {
            printf ("\n Ion %d Force  %10.7f  %10.7f  %10.7f",
                    ion, ct.ions[ion].force[ct.fpt[0]][0] - old_force[3 * ion],
                    ct.ions[ion].force[ct.fpt[0]][1] - old_force[3 * ion + 1],
                    ct.ions[ion].force[ct.fpt[0]][2] - old_force[3 * ion + 2]);

            sumx += ct.ions[ion].force[ct.fpt[0]][0] - old_force[3 * ion];
            sumy += ct.ions[ion].force[ct.fpt[0]][1] - old_force[3 * ion + 1];
            sumz += ct.ions[ion].force[ct.fpt[0]][2] - old_force[3 * ion + 2];

            old_force[3 * ion] = ct.ions[ion].force[ct.fpt[0]][0];
            old_force[3 * ion + 1] = ct.ions[ion].force[ct.fpt[0]][1];
            old_force[3 * ion + 2] = ct.ions[ion].force[ct.fpt[0]][2];
        }
        printf ("\n Non-local sums in x, y and z directions: %e %e %e", sumx, sumy, sumz);
    }
#endif



    /* The non-linear core correction part if any */
    nlccforce (rho, vxc);

#if VERBOSE
    if (pct.imgpe == 0)
    {
        printf ("\n\n Non-linear core force:");

        sumx = ZERO;
        sumy = ZERO;
        sumz = ZERO;
        for (ion = 0; ion < ct.num_ions; ion++)
        {
            printf ("\n Ion %d Force  %10.7f  %10.7f  %10.7f",
                    ion, ct.ions[ion].force[ct.fpt[0]][0] - old_force[3 * ion],
                    ct.ions[ion].force[ct.fpt[0]][1] - old_force[3 * ion + 1],
                    ct.ions[ion].force[ct.fpt[0]][2] - old_force[3 * ion + 2]);

            sumx += ct.ions[ion].force[ct.fpt[0]][0] - old_force[3 * ion];
            sumy += ct.ions[ion].force[ct.fpt[0]][1] - old_force[3 * ion + 1];
            sumz += ct.ions[ion].force[ct.fpt[0]][2] - old_force[3 * ion + 2];
        }
        printf ("\n Non-linear core force sums in x, y and z directions: %e %e %e", sumx, sumy,
                sumz);
    }
#endif

#if !GAMMA_PT
    /* Now symmetrize the forces */
    if (!(ct.kp[0].kpt[0] == 0.0 && ct.kp[0].kpt[1] == 0.0 && ct.kp[0].kpt[2] == 0.0))
        symforce ();
#endif
    my_free (vtot);

    /* Impose force constraints, if any */
    if( ct.constrainforces )
        constrain ();

#if VERBOSE
    my_free (old_force);
#endif

    rmg_timings (FORCE_TIME, (my_crtc () - time3));

}                               /* end force */


/******/
