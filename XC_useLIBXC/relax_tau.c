#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "main.h"





void relax_tau (int steps, STATE * states, double * vxc, double * vh, double * vnuc,
              double * rho, double * rho_oppo, double * rhocore, double * rhoc, double * tau)
{

    int iion;
    int CONV_FORCE, MAX_STEPS;
    int DONE, rlx_steps = 0;

    /* if ( ct.override_atoms == 1 ) quench(states, vxc, vh, vnuc, rho, rhocore, rhoc); */


	/* quench the electrons and calculate forces */
    quench_tau (states, vxc, vh, vnuc, rho, rho_oppo, rhocore, rhoc, tau);
    

    /* ---------- begin relax loop --------- */
    DONE = (ct.max_md_steps < 1 || steps < 1);

    while (!DONE)
    {

		rlx_steps++;

        if (pct.imgpe == 0)
            printf ("\nrelax: ---------- [rlx: %d/%d] ----------\n", rlx_steps, steps);

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
                error_handler ("Undefined MD method");
        }

        /* ct.md_steps measures the number of updates to the atomic positions */
        ct.md_steps++;

        /* Update items that change when the ionic coordinates change */
        reinit_ionic_pp (states, vnuc, rhocore, rhoc);

        /* quench the electrons and calculate forces */
        quench_tau (states, vxc, vh, vnuc, rho, rho_oppo, rhocore, rhoc, tau);


        /* save data to file for future restart */
        if (ct.checkpoint)
            if ( ct.md_steps % ct.checkpoint == 0 )
                write_restart (ct.outfile, vh, rho, rho_oppo, vxc, states);


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

        printf ("\n");
        progress_tag ();

        if (CONV_FORCE)
            printf ("force convergence has been achieved. stopping ...\n");
        else
            printf ("force convergence has NOT been achieved. stopping (max number of relax steps reached) ...\n");

    }



    /*Write out final data */
    write_restart (ct.outfile, vh, rho, rho_oppo, vxc, states);


}                               /* end fastrlx */


/******/
