/*
 *
 * Copyright 2014 The RMG Project Developers. See the COPYRIGHT file 
 * at the top-level directory of this distribution or in the current
 * directory.
 * 
 * This file is part of RMG. 
 * RMG is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * any later version.
 *
 * RMG is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
*/

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
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
#include "Atomic.h"
#include "RmgParallelFft.h"
#include "bfgs.h"


// Instantiate gamma and non-gamma versions
template void Relax<double>(int , double *, double *, double *,
              double *, double *, double *, double *, Kpoint<double> **Kptr);
template void Relax<std::complex<double> >(int , double *, double *, double *,
              double *, double *, double *, double *, Kpoint<std::complex<double> >**Kptr);



template <typename OrbitalType> void Relax (int steps, double * vxc, double * vh, double * vnuc,
              double * rho, double * rho_oppo, double * rhocore, double * rhoc, Kpoint<OrbitalType> **Kptr)
{

    int CONV_FORCE=false, MAX_STEPS;
    int DONE, rlx_steps = 0;
    static double *rhodiff;

    /* quench the electrons and calculate forces */
    if((ct.runflag != RESTART) || (ct.forceflag == MD_QUENCH))
    {
        Quench (vxc, vh, vnuc, rho, rho_oppo, rhocore, rhoc, Kptr, true);
        WriteRestart (ct.outfile, vh, rho, rho_oppo, vxc, Kptr);
    }


    /* ---------- begin relax loop --------- */
    DONE = (ct.max_md_steps < 1 || steps < 1);

    while (!DONE)
    {

		rlx_steps++;

        if (pct.imgpe == 0)
            rmg_printf ("\nrelax: ---------- [rlx: %d/%d] ----------\n", rlx_steps, steps);

        // Get atomic rho for this ionic configuration and subtract from current rho
        int FP0_BASIS = Rmg_G->get_P0_BASIS(Rmg_G->default_FG_RATIO);
        int factor = ct.noncoll_factor * ct.noncoll_factor;
        double *arho = new double[FP0_BASIS * factor];
        LcaoGetAtomicRho(arho);

        // If first step allocate rhodiff
        for(int idx = 0;idx < FP0_BASIS * factor;idx++) rho[idx] -= arho[idx];

        if(rhodiff == NULL)
        {
            rhodiff = new double[FP0_BASIS * factor];
            for(int idx = 0;idx < FP0_BASIS * factor;idx++) rhodiff[idx] = rho[idx];
        }
        else
        {
            double *trho = new double[FP0_BASIS * factor];
            for(int idx = 0;idx < FP0_BASIS * factor;idx++) trho[idx] = rho[idx];
            for(int idx = 0;idx < FP0_BASIS * factor;idx++) rho[idx] = 2.0*rho[idx] - rhodiff[idx];
            for(int idx = 0;idx < FP0_BASIS * factor;idx++) rhodiff[idx] = trho[idx];
            delete [] trho;
        }

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
                break;
            case MD_MIN:
                fastrelax (&ct.iondt, ct.iondt_max, ct.iondt_inc, ct.iondt_dec, ct.relax_steps_delay, &ct.relax_steps_counter);
                break;
            case LBFGS:
                //rmg_lbfgs();
                simple_lbfgs();
                break;
            default:
                rmg_error_handler (__FILE__, __LINE__, "Undefined MD method");
        }


        if(ct.cell_relax)
        {

            Reinit (vh, rho, rho_oppo, rhocore, rhoc, vnuc, vxc, Kptr);
        }
        // Get atomic rho for new configuration and add back to rho
        /* Update items that change when the ionic coordinates change */
        RmgTimer *RT0=new RmgTimer("1-TOTAL: run: ReinitIonicPotentials");
        ReinitIonicPotentials (Kptr, vnuc, rhocore, rhoc);
        delete RT0;

        LcaoGetAtomicRho(arho);
         for(int idx = 0;idx < FP0_BASIS * factor;idx++) rho[idx] += arho[idx];
        delete [] arho;

        // Reset mixing
        MixRho(NULL, NULL, NULL, NULL, NULL, NULL, Kptr[0]->ControlMap, true);

        /* ct.md_steps measures the number of updates to the atomic positions */
        ct.md_steps++;

        /* quench the electrons and calculate forces */
        Quench (vxc, vh, vnuc, rho, rho_oppo, rhocore, rhoc, Kptr, true);


        /* save data to file for future restart */
        WriteRestart (ct.outfile, vh, rho, rho_oppo, vxc, Kptr);

        // Extrapolate orbitals after first step
        ExtrapolateOrbitals(ct.outfile, Kptr);

        /* check force convergence */
        CONV_FORCE = TRUE;
        for (size_t ion = 0, i_end = Atoms.size(); ion < i_end; ++ion)
        {
            double *fp, fsq;
            fp = Atoms[ion].force[ct.fpt[0]];
            fsq = 0.0;
            for(int i = 0; i < 3; i++)
                if (Atoms[ion].movable[i])
                    fsq += fp[i] * fp[i];
            CONV_FORCE &= (fsq < ct.thr_frc * ct.thr_frc);
        }
    
        if(ct.cell_relax)
        {
            for(int i = 0; i < 9; i++) 
            {
                if(ct.cell_movable[i])  CONV_FORCE &= (std::abs(Rmg_L.stress_tensor[i]) * Ha_Kbar < ct.thr_stress);
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

}                               /* end fastrlx */


/******/
