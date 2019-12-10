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
#include "Exx_on.h"
#include "GpuAlloc.h"

void quench(STATE * states, STATE * states1, double * vxc, double * vh,
            double * vnuc, double * vh_old, double *
vxc_old, double * rho, double * rho_oppo, double * rhoc, double * rhocore)
{
    static int CONVERGENCE = FALSE;
    ct.exx_convergence_factor = 1.0; 
    Functional *F = new Functional ( *Rmg_G, Rmg_L, *Rmg_T, ct.is_gamma);

    double start_time = my_crtc ();
    double exx_start_time = start_time;
    double step_time;
    double exx_step_time;

    int outer_steps = 1;

    if(ct.xc_is_hybrid)
    {
        outer_steps = ct.max_exx_steps;
        Exx_onscf = new Exx_on<double>(*Rmg_G, *Rmg_halfgrid, Rmg_L, "temwave", *LocalOrbital, ct.exx_mode);
    }

    ct.FOCK = 0.0;
    ct.exx_delta = DBL_MAX;
    double f0=0.0,f1,f2=0.0, exxen = 0.0;
    for(ct.exx_steps = 0;ct.exx_steps < outer_steps;ct.exx_steps++)
    {

        exx_step_time = my_crtc ();

       // Adjust exx convergence threshold
        if(ct.xc_is_hybrid)
        {
            if(ct.exx_steps)
                ct.exx_convergence_factor = std::min(1.0e-8, std::max(1.0e-15, fabs(ct.exx_delta)) / 1000000.0) / ct.thr_energy;
            else
                ct.exx_convergence_factor = 1.0e-7 / ct.thr_energy;
            if(fabs(ct.exx_delta) < 1.0e-6) ct.exx_convergence_factor /= 10.0;
        }



        for (ct.scf_steps = 0; ct.scf_steps < ct.max_scf_steps; ct.scf_steps++)
        {
            if (pct.gridpe == 0)
                printf("\n\n\n ITERATION     %d\n", ct.scf_steps);

            step_time = my_crtc();

            bool freeze_orbital = false;

            if(ct.scf_steps >= ct.freeze_orbital_step) freeze_orbital = true;
            if(ct.exx_steps >= ct.exx_freeze_orbital_step) freeze_orbital = true;
            /* Perform a single self-consistent step */
            if (!CONVERGENCE || ct.scf_steps <= ct.freeze_rho_steps)
            {
                if(ct.LocalizedOrbitalLayout == LO_projection)
                {
                    Scf_on_proj(states, vxc, vh, vnuc, rho, rho_oppo, rhoc, 
                            rhocore, vxc_old, vh_old, &CONVERGENCE, freeze_orbital);
                    //Scf_on(states, states1, vxc, vh, vnuc, rho, rho_oppo, rhoc, 
                    //        rhocore, vxc_old, vh_old, &CONVERGENCE);
                }
                else
                {
                    Scf_on(states, states1, vxc, vh, vnuc, rho, rho_oppo, rhoc, 
                            rhocore, vxc_old, vh_old, &CONVERGENCE);
                }
            }
            step_time = my_crtc() - step_time;

            if( (ct.scf_steps+1)%ct.checkpoint == 0)
            {
                write_restart(ct.outfile, vh, vxc, vh_old, vxc_old, rho, rho_oppo, &states[0]);
            }


            if (CONVERGENCE && ct.scf_steps > ct.freeze_rho_steps)
            {
                if (pct.gridpe == 0)
                    printf ("\n\n Convergence has been achieved. stopping ...\n");
                break;
            }
            
            // Write out progress info
            double elapsed_time = my_crtc() - start_time;
            ProgressTag(step_time, elapsed_time);

        }

        if(ct.xc_is_hybrid)
        {
            
            size_t size = LocalOrbital->num_thispe * LocalOrbital->pbasis * sizeof(double);
            memcpy(Exx_onscf->PreOrbital, LocalOrbital->storage_proj, size);
            F->start_exx_rmg();

            double *rho_mat_global = (double *)GpuMallocManaged(LocalOrbital->num_tot * LocalOrbital->num_tot * sizeof(double) + 8);
            double *rho_mat_local = (double *)GpuMallocManaged(LocalOrbital->num_tot * LocalOrbital->num_thispe * sizeof(double) + 8);
            double *Sij_inverse = (double *)GpuMallocManaged(LocalOrbital->num_tot * LocalOrbital->num_tot * sizeof(double) + 8);

            mat_dist_to_global(mat_X, pct.desca, rho_mat_global);
            mat_dist_to_global(matB, pct.desca, Sij_inverse);

            
            if(ct.nspin == 1) 
            {
                double half = 0.5; 
                int ione = 1, nmat = LocalOrbital->num_tot * LocalOrbital->num_tot;
                dscal(&nmat, &half, rho_mat_global, &ione);
            }

            for(int i = 0; i < LocalOrbital->num_thispe; i++)
            {
                int i_glob = LocalOrbital->index_proj_to_global[i];
                for(int j = 0; j < LocalOrbital->num_tot; j++)
                    rho_mat_local[j * LocalOrbital->num_thispe + i] = rho_mat_global[j * LocalOrbital->num_tot + i_glob];
            }

            //f1 = Exx_onscf->Exxenergy(rho_mat_global);
            Exx_onscf->Omega(rho_mat_local, (std::abs(ct.exx_delta) > ct.vexx_fft_threshold) );
            Exx_onscf->Xij(Sij_inverse, *LocalOrbital);
            Exx_onscf->OmegaSinv(Sij_inverse, *LocalOrbital);
            f2 = Exx_onscf->Exxenergy(rho_mat_global);
            ct.exx_delta = f2 - f0;
            ct.FOCK = f2;

            f0 = f2;
            exx_step_time = my_crtc () - exx_step_time;
            // Write out progress info
            double exx_elapsed_time = my_crtc() - exx_start_time;
            ExxProgressTag(exx_step_time, exx_elapsed_time);

            GpuFreeManaged(rho_mat_local);
            GpuFreeManaged(rho_mat_global);
            GpuFreeManaged(Sij_inverse);

            if(ct.exx_steps == 0)
                UpdatePot(vxc, vh, vxc_old, vh_old, vnuc, rho, rho_oppo, rhoc, rhocore);

            if(fabs(ct.exx_delta) < ct.exx_convergence_criterion)
            {
                printf(" Finished EXX outer loop in %3d exx steps exx_delta = %8.2e, total energy = %.*f Ha\n",
                        ct.exx_steps, ct.exx_delta, 6, ct.TOTAL);
                ct.FOCK = f2;
                break;
            }
            else
            {
                printf(" Finished EXX inner loop in %3d scf steps exx_delta = %8.2e, total energy = %.*f Ha\n",
                        ct.scf_steps, ct.exx_delta, 6, ct.TOTAL);
            }


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
        double *psi = new double[2* LocalOrbital->num_tot * pbasis]();
        double *Cij_global = new double[LocalOrbital->num_tot * LocalOrbital->num_tot];
        double *Cij_local = new double[LocalOrbital->num_tot * LocalOrbital->num_thispe + 1];

        mat_dist_to_global(zz_dis, pct.desca, Cij_global);

        for(int i = 0; i < LocalOrbital->num_thispe; i++)
        {
            int i_glob = LocalOrbital->index_proj_to_global[i];
            for(int j = 0; j < LocalOrbital->num_tot; j++)
                Cij_local[j * LocalOrbital->num_thispe + i] = Cij_global[j * LocalOrbital->num_tot + i_glob];
        }


        double one(1.0), zero(0.0);
        if(LocalOrbital->num_thispe > 0)
            dgemm("N", "N", &pbasis, &LocalOrbital->num_tot, &LocalOrbital->num_thispe, &one, LocalOrbital->storage_proj, &pbasis,
                    Cij_local, &LocalOrbital->num_thispe, &zero, psi, &pbasis);


        Exxbase<double> Exx(*Rmg_G, *Rmg_halfgrid, Rmg_L, "tempwave", nstates_occ, occs.data(), psi, ct.exx_mode);
        if(ct.exx_mode == EXX_LOCAL_FFT)
            Exx.WriteWfsToSingleFile();

        MPI_Barrier(MPI_COMM_WORLD);
        Exx.Vexx_integrals(ct.exx_int_file);
    }


}
