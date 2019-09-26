/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
/****f* QMD-MGDFT/quench.c *****
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
 *   void quench(STATE *states, double *vxc, double *vh, double *vnuc, 
 *               double *rho, double *rhocore, double *rhoc)
 *   For a fixed atomic configuration, quench the electrons to find 
 *   the minimum total energy 
 * INPUTS
 *   states: point to orbital structure (see main.h)
 *   vxc:    exchange correlation potential
 *   vh:     Hartree potential
 *   vnuc:   Pseudopotential 
 *   rho:    total valence charge density
 *   rhocore: core chare density only for non-linear core correction
 *   rhoc:   Gaussian compensating charge density
 * OUTPUT
 *   states, vxc, vh, rho are updated
 * PARENTS
 *   cdfastrlx.c fastrlx.c md.c
 * CHILDREN
 *   scf.c force.c get_te.c subdiag.c get_ke.c
 * SOURCE
 */

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "main.h"
#include "init_var.h"
#include "prototypes_on.h"
#include "transition.h"
#include "Exxbase.h"

void quench(STATE * states, STATE * states1, double * vxc, double * vh,
            double * vnuc, double * vh_old, double *
vxc_old, double * rho, double * rho_oppo, double * rhoc, double * rhocore)
{
    int outcount = 0;
    static int CONVERGENCE = FALSE;



    for (ct.scf_steps = 0; ct.scf_steps < ct.max_scf_steps; ct.scf_steps++)
    {
        if (pct.gridpe == 0)
            printf("\n\n\n ITERATION     %d\n", ct.scf_steps);

        /* Perform a single self-consistent step */
        if (!CONVERGENCE || ct.scf_steps <= ct.freeze_rho_steps)
        {
            if(ct.LocalizedOrbitalLayout == LO_projection)
            {
                //Scf_on_proj(states, vxc, vh, vnuc, rho, rho_oppo, rhoc, 
                //        rhocore, vxc_old, vh_old, &CONVERGENCE);
                Scf_on(states, states1, vxc, vh, vnuc, rho, rho_oppo, rhoc, 
                        rhocore, vxc_old, vh_old, &CONVERGENCE);
            }
            else
            {
                Scf_on(states, states1, vxc, vh, vnuc, rho, rho_oppo, rhoc, 
                        rhocore, vxc_old, vh_old, &CONVERGENCE);
            }
        }

        if( (ct.scf_steps+1)%ct.checkpoint == 0)
            write_restart(ct.outfile, vh, vxc, vh_old, vxc_old, rho, rho_oppo, &states[0]);

        if (CONVERGENCE && ct.scf_steps > ct.freeze_rho_steps)
        {
            if (pct.gridpe == 0)
                printf ("\n\n Convergence has been achieved. stopping ...\n");
            break;
        }


    }

    if (pct.gridpe == 0)
        rmg_printf("\n final total energy = %14.8f Ha\n", ct.TOTAL);
    // Exact exchange integrals
    if(ct.exx_int_flag)
    {
        int nstates_occ = 0;
        std::vector<double> occs;
        // calculate num of occupied states
        for(int i=0;i < ct.num_states;i++) 
        {
            if(states[i].occupation[0] > 1.0e-6)
            {
                occs.push_back(states[i].occupation[0]);
                nstates_occ++;
            }
            else
            {
                //    break;
            }
        }

        int pbasis = Rmg_G->get_P0_BASIS(1);
        //  calculate the wavefuctions from localized orbitals
        double *psi = new double[2* LocalOrbital->num_tot * pbasis];
        double *Cij_global = new double[LocalOrbital->num_tot * LocalOrbital->num_tot];
        double *Cij_local = new double[LocalOrbital->num_tot * LocalOrbital->num_thispe];

        mat_dist_to_global(zz_dis, pct.desca, Cij_global);

        for(int i = 0; i < LocalOrbital->num_thispe; i++)
        {
            int i_glob = LocalOrbital->index_proj_to_global[i];
            for(int j = 0; j < LocalOrbital->num_tot; j++)
                Cij_local[j * LocalOrbital->num_thispe + i] = Cij_global[j * LocalOrbital->num_tot + i_glob];
        }


        double one(1.0), zero(0.0);
        dgemm("N", "N", &pbasis, &LocalOrbital->num_tot, &LocalOrbital->num_thispe, &one, LocalOrbital->storage_proj, &pbasis,
                Cij_local, &LocalOrbital->num_thispe, &zero, psi, &pbasis);


        Exxbase<double> Exx(*Rmg_G, Rmg_L, "tempwave", nstates_occ, occs.data(), psi, ct.exx_mode);
        if(ct.exx_mode == EXX_LOCAL_FFT)
            Exx.WriteWfsToSingleFile();

        MPI_Barrier(MPI_COMM_WORLD);
        Exx.Vexx_integrals(ct.exx_int_file);
    }


}
