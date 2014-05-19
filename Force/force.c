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
 *   void force(rmg_double_t *rho, rmg_double_t *rhoc, rmg_double_t *vh, rmg_double_t *vxc, STATE *states)
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
#include "common_prototypes.h"


/*Set this to 1 to have forces written out part by part*/
/* If you want this , you should also make sure that VERBOSE flag is enabled in
 * nlforce.c*/
#define VERBOSE 0



void force (rmg_double_t * rho, rmg_double_t * rho_oppo, rmg_double_t * rhoc, rmg_double_t * vh, rmg_double_t * vxc, rmg_double_t * vnuc, STATE * states)
{
    int ion, idx;
    rmg_double_t *vtott, *rho_tot;
#if VERBOSE
    rmg_double_t *old_force;
    rmg_double_t sumx, sumy, sumz;
#endif

    int Zi;

#if VERBOSE
    my_malloc (old_force, 3 * ct.num_ions, rmg_double_t);

    for (ion = 0; ion < ct.num_ions; ion++)
    {
        old_force[3 * ion] = ZERO;
        old_force[3 * ion + 1] = ZERO;
        old_force[3 * ion + 2] = ZERO;
    }
#endif

    my_malloc (vtott, get_FP0_BASIS(), rmg_double_t);
    for (idx = 0; idx < get_FP0_BASIS(); idx++)
        vtott[idx] = vxc[idx] + vh[idx] + vnuc[idx];


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
    	my_malloc (rho_tot, get_FP0_BASIS(), rmg_double_t);
	for (idx = 0; idx < get_FP0_BASIS(); idx++)
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


    /* Add in the non-local stuff */
    nlforce (vtott);



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
    my_free (vtott);

    /* Impose force constraints, if any */
    if( ct.constrainforces )
        constrain ();

#if VERBOSE
    my_free (old_force);
#endif


}                               /* end force */


/******/
