/************************** SVN Revision Information **************************
 **    $Id: relax.c 2012 2013-05-13 17:22:57Z ebriggs $    **
******************************************************************************/

/****f* QMD-MGDFT/fastrlx.c *****
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
 *   void fastrlx(STATE *states, double *vxc, double *vh, double *vnuc,
 *                double *rho, double *rhocore, double *rhoc)
 *   drive routine for fast relax.
 * INPUTS
 *   states: all wave functions (see main.h)
 *   vxc: exchange-correlation potential
 *   vh:  Hartree potential
 *   vnuc: pseudopotential
 *   rho:  total charge density
 *   rhocore: charge density of core electrons, only useful when we 
 *            include non-linear core correction for pseudopotential.
 *   rhoc:    compensating charge density
 * OUTPUT
 *   all the above inputs are updated
 * PARENTS
 *   main.c
 * CHILDREN
 *   quench.c to_crystal.c to_cartesian.c init_nuc.c get_nlop.c scf.c
 *   get_te.c force.c rmg_fastrelax.c
 * SOURCE
 */

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/mman.h>
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
#include "transition.h"
#include "../Headers/prototypes.h"


// Instantiate gamma and non-gamma versions
template void Relax<double>(int , double *, double *, double *,
              double *, double *, double *, double *, Kpoint<double> **Kptr);
template void Relax<std::complex<double> >(int , double *, double *, double *,
              double *, double *, double *, double *, Kpoint<std::complex<double> >**Kptr);



template <typename OrbitalType> void Relax (int steps, double * vxc, double * vh, double * vnuc,
              double * rho, double * rho_oppo, double * rhocore, double * rhoc, Kpoint<OrbitalType> **Kptr)
{

    int iion;
    int CONV_FORCE, MAX_STEPS;
    int DONE, rlx_steps = 0;

    /* if ( ct.override_atoms == 1 ) quench(states, vxc, vh, vnuc, rho, rhocore, rhoc); */

	/* quench the electrons and calculate forces */
    Quench (vxc, vh, vnuc, rho, rho_oppo, rhocore, rhoc, Kptr);

    /* ---------- begin relax loop --------- */
    DONE = (ct.max_md_steps < 1 || steps < 1);

    while (!DONE)
    {

		rlx_steps++;

        if (pct.imgpe == 0)
            rmg_printf ("\nrelax: ---------- [rlx: %d/%d] ----------\n", rlx_steps, steps);

        /* not done yet ? => move atoms */
		/* move the ions */
        switch(ct.relax_method)
        {

            case FASTRELAX:
                fastrelax (&ct.iondt, ct.iondt_max, ct.iondt_inc, ct.iondt_dec, ct.relax_steps_delay, &ct.relax_steps_counter);
                break;
            case FIRE:
                fire (&ct.iondt, ct.iondt_max, ct.iondt_inc, ct.iondt_dec, ct.relax_steps_delay, &ct.relax_steps_counter);
                break;
            case QUICK_MIN:
                fastrelax (&ct.iondt, ct.iondt_max, ct.iondt_inc, ct.iondt_dec, ct.relax_steps_delay, &ct.relax_steps_counter);
                //quick_min ();
                break;
            case MD_MIN:
                fastrelax (&ct.iondt, ct.iondt_max, ct.iondt_inc, ct.iondt_dec, ct.relax_steps_delay, &ct.relax_steps_counter);
                //md_min ();
                break;
            case LBFGS:
                rmg_lbfgs();
                break;
            default:
                rmg_error_handler (__FILE__, __LINE__, "Undefined MD method");
        }

        /* ct.md_steps measures the number of updates to the atomic positions */
        ct.md_steps++;

        /* Update items that change when the ionic coordinates change */
        ReinitIonicPotentials (Kptr, vnuc, rhocore, rhoc);

        /* quench the electrons and calculate forces */
        Quench (vxc, vh, vnuc, rho, rho_oppo, rhocore, rhoc, Kptr);


        /* save data to file for future restart */
        if (ct.checkpoint)
            if ( ct.md_steps % ct.checkpoint == 0 )
                WriteRestart (ct.outfile, vh, rho, rho_oppo, vxc, Kptr);


        /* check force convergence */
        CONV_FORCE = TRUE;
        for (iion = 0; iion < ct.num_ions; iion++)
        {
            if (ct.ions[iion].movable)
            {
                double *fp;
                fp = ct.ions[iion].force[ct.fpt[0]];
                CONV_FORCE &=
                    ((fp[0] * fp[0] + fp[1] * fp[1] + fp[2] * fp[2]) < ct.thr_frc * ct.thr_frc);
            }
        }

        /* check for max relax steps */
        MAX_STEPS = (rlx_steps >= steps) || ( ct.md_steps > ct.max_md_steps);

        /* done if forces converged or reached limit of md steps */
        DONE = (CONV_FORCE || MAX_STEPS);

    }
    /* ---------- end relax loop --------- */

    if (ct.max_md_steps > 0 && steps > 0)
    {

        rmg_printf ("\n");
        //progress_tag ();

        if (CONV_FORCE)
            rmg_printf ("force convergence has been achieved. stopping ...\n");
        else
            rmg_printf ("force convergence has NOT been achieved. stopping (max number of relax steps reached) ...\n");

    }



    /*Write out final data */
    WriteRestart (ct.outfile, vh, rho, rho_oppo, vxc, Kptr);


}                               /* end fastrlx */


/******/
