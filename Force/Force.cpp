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
 *   void force(double *rho, double *rhoc, double *vh, double *vxc, STATE *states)
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
#include "transition.h"
#include "const.h"
#include "RmgTimer.h"
#include "rmgtypedefs.h"
#include "params.h"
#include "typedefs.h"
#include "common_prototypes.h"
#include "common_prototypes1.h"
#include "rmg_error.h"
#include "Kpoint.h"


/*Set this to 1 to have forces written out part by part*/
/* If you want this , you should also make sure that VERBOSE flag is enabled in
 * nlforce.c*/
#define VERBOSE 0

template void Force<double> (double * rho, double * rho_oppo, double * rhoc, double * vh, 
        double * vxc, double * vnuc, Kpoint<double> **Kptr);
template void Force<std::complex<double> > (double * rho, double * rho_oppo, double * rhoc, double * vh, 
        double * vxc, double * vnuc, Kpoint<std::complex<double>> **Kptr);


template <typename OrbitalType> void Force (double * rho, double * rho_oppo, double * rhoc, double * vh, 
        double * vxc, double * vnuc, Kpoint<OrbitalType> **Kptr)
{
    RmgTimer RT0("Force");
    int ion, idx;
    double *vtott, *rho_tot;
#if VERBOSE
    double *old_force;
    double sumx, sumy, sumz;
#endif

    int Zi;

#if VERBOSE
    old_force = new double[ 3 * ct.num_ions];

    for (ion = 0; ion < ct.num_ions; ion++)
    {
        old_force[3 * ion] = ZERO;
        old_force[3 * ion + 1] = ZERO;
        old_force[3 * ion + 2] = ZERO;
    }
#endif

    vtott = new double[get_FP0_BASIS()];

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
    RmgTimer *RT1 = new RmgTimer("Force: ion-ion");
    iiforce ();
    delete RT1;

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
    RmgTimer *RT2 = new RmgTimer("Force: local");
    if (ct.spin_flag)
    {
        rho_tot = new double[get_FP0_BASIS()];
        for (idx = 0; idx < get_FP0_BASIS(); idx++)
            rho_tot[idx] = rho[idx] + rho_oppo[idx];
        lforce(rho_tot, vh);
        delete[] rho_tot;
    }
    else
        lforce (rho, vh);
    delete RT2;


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
    RmgTimer *RT3 = new RmgTimer("Force: non-local");
    Nlforce (vtott, Kptr);
    delete RT3;


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
    RmgTimer *RT4 = new RmgTimer("Force: 1core correction");
    nlccforce (rho, vxc);
    delete RT4;

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

    if (!ct.is_gamma) {
        RmgTimer *RT5 = new RmgTimer("Force: symmetrization");
        symforce ();
        delete RT5;
    }
    delete[] vtott;

    /* Impose force constraints, if any */
    if( ct.constrainforces )
        constrain ();

#if VERBOSE
    delete[] old_force;
#endif


}                               /* end force */


/******/
