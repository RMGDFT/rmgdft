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
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include "transition.h"
#include "const.h"
#include "State.h"
#include "Kpoint.h"
#include "TradeImages.h"
#include "RmgTimer.h"
#include "RmgThread.h"
#include "GlobalSums.h"
#include "rmgthreads.h"
#include "vhartree.h"
#include "packfuncs.h"
#include "typedefs.h"
#include "common_prototypes.h"
#include "common_prototypes1.h"
#include "rmg_error.h"
#include "Kpoint.h"
#include "Subdiag.h"
#include "Functional.h"
#include "Solvers.h"
#include "RmgParallelFft.h"




template bool Scf<double> (double *, double *, double *, double *,
          double *, double *, double *, double *, double *, double *, int ,
          int , Kpoint<double> **, std::vector<double> &);
template bool Scf<std::complex<double> > (double *, double *, double *, double *,
          double *, double *, double *, double *, double *, double *, int ,
          int , Kpoint<std::complex<double>> **, std::vector<double> &);

template <typename OrbitalType> bool Scf (double * vxc, double *vxc_in, double * vh, double *vh_in, double *vh_ext,
          double * vnuc, double * rho, double * rho_oppo, double * rhocore, double * rhoc, int spin_flag,
          int boundaryflag, Kpoint<OrbitalType> **Kptr, std::vector<double>& RMSdV)
{

    RmgTimer RT0("2-Scf steps"), *RT1;
    RmgTimer RTt("1-TOTAL: run: Scf steps");

    bool CONVERGED = false;
    double t3;
    double *vtot, *vtot_psi, *new_rho, *new_rho_oppo = NULL;
    double *vxc_psi=NULL;
    double t[3];                  /* SCF checks and average potential */
    double max_unocc_res = 0.0;

    int dimx = Rmg_G->get_PX0_GRID(Rmg_G->get_default_FG_RATIO());
    int dimy = Rmg_G->get_PY0_GRID(Rmg_G->get_default_FG_RATIO());
    int dimz = Rmg_G->get_PZ0_GRID(Rmg_G->get_default_FG_RATIO());
    int FP0_BASIS = dimx * dimy * dimz;
    int P0_BASIS;

    double *rho_tf = NULL;
    double *rho_save = NULL; 
    

    P0_BASIS =  Rmg_G->get_P0_BASIS(1);
    FP0_BASIS = Rmg_G->get_P0_BASIS(Rmg_G->default_FG_RATIO);

    /* allocate memory for eigenvalue send array and receive array */
    new_rho = new double[ct.nspin * FP0_BASIS];
    if(ct.nspin == 2) new_rho_oppo = &new_rho[FP0_BASIS];
    vtot = new double[FP0_BASIS];
    vtot_psi = new double[P0_BASIS];

    /* save old vhxc + vnuc */
    for (int idx = 0; idx < FP0_BASIS; idx++) {
        vtot[idx] = vxc[idx] + vh[idx] + vnuc[idx];
    }

    if(ct.noncoll) vxc_psi = new double[4*P0_BASIS];


    /* Evaluate XC energy and potential */
    RT1 = new RmgTimer("2-Scf steps: exchange/correlation");
    Functional *F = new Functional ( *Rmg_G, Rmg_L, *Rmg_T, ct.is_gamma);
    F->v_xc(rho, rhocore, ct.XC, ct.vtxc, vxc, ct.nspin );
    if(ct.filter_dpot && (Rmg_G->default_FG_RATIO > 1)) 
    {
        for(int is = 0; is < ct.nspin; is++)
            FftFilter(&vxc[is*FP0_BASIS], *fine_pwaves, sqrt(ct.filter_factor) / (double)ct.FG_RATIO, LOW_PASS);
    }
    //if(pct.gridpe==0)printf("\nXC = %f  %f\n", ct.XC, ct.vtxc);
    delete F;
    delete RT1;

    double rms_target = std::min(std::max(ct.rms/ct.hartree_rms_ratio, 1.0e-12), 1.0e-6);

    /*Simplified solvent model, experimental */
    if (ct.num_tfions > 0)
    {

        rho_tf = new double[FP0_BASIS];
        rho_save = new double[FP0_BASIS];

        //get_dipole (rho);
        rmg_printf("\nCalling get_tf_rho\n");
        get_tf_rho(rho_tf);
        rmg_printf("get_tf_rho Done\n");

        for (int idx = 0; idx < FP0_BASIS; idx++)
        {
            /*Save original rho, copy it into rho_save*/
            rho_save[idx] = rho[idx];

            /*Add rho_tf to rho*/
            rho[idx] += rho_tf[idx];

            /*Check that rho_tf chargeis indeed 0, can be removed if it works well*/
            t[0] += rho_tf[idx];
        }                           /* idx */
    }


    RT1 = new RmgTimer("2-Scf steps: Hartree");
    double hartree_residual = VhDriver(rho, rhoc, vh, vh_ext, rms_target);
    delete(RT1);

    /*Simplified solvent model, experimental */
    if (ct.num_tfions > 0)
    {
        for (int idx = 0; idx < FP0_BASIS; idx++)
        { 
            rho[idx] = rho_save[idx];
        }
    }


    // Save input hartree potential
    for (int idx = 0; idx < FP0_BASIS; idx++) vh_in[idx] = vh[idx];

    /* check convergence */
    t[0] = t[1] = t[2] = 0.0;

    for (int idx = 0; idx < FP0_BASIS; idx++)
    {
        t3 = -vtot[idx];
        vtot[idx] = vxc[idx] + vh[idx] + vnuc[idx];
        t3 += vtot[idx];
        t[0] += rho[idx] * t3;
        t[1] += t3 * t3;
        t[2] += vh[idx];
    }                           /* idx */

    MPI_Allreduce(MPI_IN_PLACE, t, 3, MPI_DOUBLE, MPI_SUM, pct.grid_comm);
    MPI_Allreduce(MPI_IN_PLACE, t, 3, MPI_DOUBLE, MPI_SUM, pct.spin_comm);
    t[0] *= get_vel_f();

    /* get the averaged value over each spin and each fine grid */
    t[1] = sqrt (t[1] / ((double) (ct.nspin * ct.psi_fnbasis)));  
    t[2] /= ((double) (ct.nspin * ct.psi_fnbasis));   

    ct.rms = t[1];

    if (ct.scf_steps)
    {
        //rmg_printf("scf check: <rho dv>   = %8.2e\n", t[0]);
        RMSdV.emplace_back(t[1]);
        if(ct.poisson_solver == MULTIGRID_SOLVER) 
            rmg_printf("hartree residual      = %8.2e\n", hartree_residual);
        rmg_printf("average potential <V> = %8.2e\n", t[2]);
    }

    if(!Verify ("freeze_occupied", true, Kptr[0]->ControlMap)) {
  	if (ct.scf_steps && t[1] < ct.thr_rms) CONVERGED = true;
    }

    // Transfer vtot from the fine grid to the wavefunction grid
    GetVtotPsi (vtot_psi, vtot, Rmg_G->default_FG_RATIO);
    if(ct.noncoll)
        for(int is = 0; is < ct.nspin; is++)
        {
            GetVtotPsi (&vxc_psi[is*P0_BASIS], &vxc[is*FP0_BASIS], Rmg_G->default_FG_RATIO);
        }

    /*Generate the Dnm_I */
    get_ddd (vtot, vxc, true);

    if(ct.ldaU_mode != LDA_PLUS_U_NONE)
    {
        RmgTimer("3-MgridSubspace: ldaUop x psi");
        int pstride = Kptr[0]->ldaU->ldaU_m;
        int occ_size = ct.nspin * Atoms.size() * pstride * pstride;
        std::complex<double> *ns_occ_g = new std::complex<double>[occ_size]();

        for (int kpt =0; kpt < ct.num_kpts_pe; kpt++)
        {
            LdaplusUxpsi(Kptr[kpt], 0, Kptr[kpt]->nstates, Kptr[kpt]->orbitalsint_local);
            Kptr[kpt]->ldaU->calc_ns_occ(Kptr[kpt]->orbitalsint_local, 0, Kptr[kpt]->nstates);
            for(int idx = 0; idx < occ_size; idx++)
                ns_occ_g[idx] += Kptr[kpt]->ldaU->ns_occ.data()[idx];
        }

        MPI_Allreduce(MPI_IN_PLACE, (double *)ns_occ_g, occ_size * 2, MPI_DOUBLE, MPI_SUM, pct.kpsub_comm);

        if(Rmg_Symm) Rmg_Symm->symm_nsocc(ns_occ_g, pstride);

        for (int kpt =0; kpt < ct.num_kpts_pe; kpt++)
        {
            for(int idx = 0; idx < occ_size; idx++)
                Kptr[kpt]->ldaU->ns_occ.data()[idx] = ns_occ_g[idx];
        }
        if(ct.verbose) Kptr[0]->ldaU->write_ldaU();


        delete ns_occ_g;
    }


    // Loop over k-points
    for(int kpt = 0;kpt < ct.num_kpts_pe;kpt++) {

        if (Verify ("kohn_sham_solver","multigrid", Kptr[0]->ControlMap) || 
                ((ct.scf_steps < 4) && (ct.md_steps == 0) && (ct.runflag != RESTART )) ||
                (ct.xc_is_hybrid && Functional::is_exx_active())) {
            RmgTimer *RT1 = new RmgTimer("2-Scf steps: MgridSubspace");
            Kptr[kpt]->MgridSubspace(vtot_psi, vxc_psi);
            delete RT1;
        }
        else if(Verify ("kohn_sham_solver","davidson", Kptr[0]->ControlMap)) {
            int notconv;
            RmgTimer *RT1 = new RmgTimer("2-Scf steps: Davidson");
            Kptr[kpt]->Davidson(vtot_psi, vxc_psi, notconv);
            delete RT1;
        }

        // Needed to ensure consistency with some types of kpoint parrelization
        MPI_Barrier(pct.grid_comm);
    } // end loop over kpoints

    if (ct.nspin == 2)
        GetOppositeEigvals (Kptr);

    for(int kpt = 0;kpt < ct.num_kpts_pe;kpt++) {
        int factor = 1;
        if(ct.nspin == 2) factor = 2;
        ct.kp[kpt+pct.kstart].eigs.resize(ct.num_states * factor);
        for(int st = 0; st < ct.num_states; st++)
            ct.kp[kpt+pct.kstart].eigs[st] = Kptr[kpt]->Kstates[st].eig[0];

        if(ct.nspin == 2)
        {
            for(int st = 0; st < ct.num_states; st++)
                ct.kp[kpt+pct.kstart].eigs[st + ct.num_states] = Kptr[kpt]->Kstates[st].eig[1];
        }
    }



    /* Take care of occupation filling */
    ct.efermi = Fill (Kptr, ct.occ_width, ct.nel, ct.occ_mix, ct.num_states, ct.occ_flag, ct.mp_order);


    if (ct.occ_flag == 1 )
    {
        rmg_printf ("\n");
        //progress_tag ();
        rmg_printf ("FERMI ENERGY = %15.8f eV\n", ct.efermi * Ha_eV);
    }

    // Calculate total energy 
    // Eigenvalues are based on in potentials and density
    GetTe (rho, rho_oppo, rhocore, rhoc, vh, vxc, Kptr, !ct.scf_steps);

    /* Generate new density */
    RT1 = new RmgTimer("2-Scf steps: GetNewRho");
    GetNewRho(Kptr, new_rho);
    if (ct.nspin == 2)
        get_rho_oppo (new_rho,  new_rho_oppo);
    delete(RT1);

    // Get Hartree potential for the output density
    int ratio = Rmg_G->default_FG_RATIO;
    double vel = Rmg_L.get_omega() / ((double)(Rmg_G->get_NX_GRID(ratio) * Rmg_G->get_NY_GRID(ratio) * Rmg_G->get_NZ_GRID(ratio)));

    /*Simplified solvent model, experimental */
    if (ct.num_tfions > 0)
    {

        for (int idx = 0; idx < FP0_BASIS; idx++)
        {
            /*Save original rho, copy it into rho_save*/
            rho_save[idx] = new_rho[idx];

            /*Add rho_tf to rho*/
            new_rho[idx] += rho_tf[idx];

        }                           /* idx */

    }

    double *vh_out = new double[FP0_BASIS];
    RT1 = new RmgTimer("2-Scf steps: Hartree");
    VhDriver(new_rho, rhoc, vh_out, vh_ext, rms_target);
    delete RT1;

    /*Simplified solvent model, experimental */
    if (ct.num_tfions > 0)
    {
        for (int idx = 0; idx < FP0_BASIS; idx++)
        { 
            new_rho[idx] = rho_save[idx];
        }

        /*GlobalSums (t, 1, pct.img_comm);
          t[0] *= get_vel_f();
          rmg_printf("\n vh_out sum (using new_rho) is %.8e\n", t[0]);*/
    }


    // Compute convergence measure (2nd order variational term) and average by nspin
    double sum = 0.0;
    for(int i = 0;i < FP0_BASIS;i++) sum += (vh_out[i] - vh[i]) * (new_rho[i] - rho[i]);
    sum = 0.5 * vel * sum;

    MPI_Allreduce(MPI_IN_PLACE, &sum, 1, MPI_DOUBLE, MPI_SUM, pct.grid_comm);
    MPI_Allreduce(MPI_IN_PLACE, &sum, 1, MPI_DOUBLE, MPI_SUM, pct.spin_comm);
    ct.scf_accuracy = sum;

    // Compute variational energy correction term if any
    sum = EnergyCorrection(Kptr, rho, new_rho, vh, vh_out);
    ct.scf_correction = sum;
    // Check if this convergence threshold has been reached
    if(!Verify ("freeze_occupied", true, Kptr[0]->ControlMap)) {
        if (ct.scf_steps && fabs(ct.scf_accuracy) < ct.adaptive_thr_energy) CONVERGED = true;
        if (ct.scf_steps && fabs(ct.scf_accuracy < 1.0e-15)) CONVERGED = true;
    }

    if(CONVERGED || (ct.scf_steps == (ct.max_scf_steps-1)))
    {
        // If the multigrid solver is selected the total energy calculation from
        // the loop above is not variational so we do a single davidson step
        // with a multiplier of 1 to get a variational value for the energy.
        if (Verify ("kohn_sham_solver","multigrid", Kptr[0]->ControlMap) && !ct.noncoll)
        {
            for (int idx = 0; idx < FP0_BASIS; idx++)
            {
                vtot[idx] = vxc[idx] + vh[idx] + vnuc[idx];
            }
            // Transfer vtot from the fine grid to the wavefunction grid
            GetVtotPsi (vtot_psi, vtot, Rmg_G->default_FG_RATIO);

            /*Generate the Dnm_I */
            get_ddd (vtot, vxc, true);

            // Set davidson parameters
            ct.davidx = 1;
            ct.david_max_steps = 0;
            ct.scf_correction = 0.0;
            int notconv;

            for(int kpt = 0;kpt < ct.num_kpts_pe;kpt++)
            {
                Kptr[kpt]->Davidson(vtot_psi, vxc_psi, notconv);
                MPI_Barrier(pct.grid_comm);
            }
            if (ct.nspin == 2) GetOppositeEigvals (Kptr);
            GetTe (rho, rho_oppo, rhocore, rhoc, vh, vxc, Kptr, !ct.scf_steps);
        }

        // Evaluate XC energy and potential from the output density
        // for the force correction
        RT1 = new RmgTimer("2-Scf steps: exchange/correlation");
        for(int i = 0;i < FP0_BASIS;i++) vxc_in[i] = vxc[i];
        Functional *F = new Functional ( *Rmg_G, Rmg_L, *Rmg_T, ct.is_gamma);
        F->v_xc(new_rho, rhocore, ct.XC, ct.vtxc, vxc, ct.nspin );
        delete F;
        delete RT1;

    }


    /*Takes care of mixing and checks whether the charge density is negative*/
    RT1 = new RmgTimer("2-Scf steps: MixRho");
    MixRho(new_rho, rho, rhocore, vh, vh_out, rhoc, Kptr[0]->ControlMap, false);
    delete RT1;

    if (ct.nspin == 2)
        get_rho_oppo (rho,  rho_oppo);



    /* Make sure we see output, e.g. so we can shut down errant runs */
    fflush( ct.logfile );
#if !(__CYGWIN__ || __MINGW64__ || _WIN32 || _WIN64)
    fsync( fileno(ct.logfile) );
#endif

    /* release memory */
    delete [] new_rho;
    delete [] vtot;
    delete [] vtot_psi;
    if(ct.noncoll) delete [] vxc_psi;

    if (ct.num_tfions > 0)
    {
        delete [] rho_tf;
        delete [] rho_save;
    }


    if(Verify ("freeze_occupied", true, Kptr[0]->ControlMap)) {

        if(ct.scf_steps && (max_unocc_res < ct.gw_threshold)) {
            rmg_printf("\nGW: convergence criteria of %10.5e has been met.\n", ct.gw_threshold);
            rmg_printf("GW:  Highest occupied orbital index              = %d\n", Kptr[0]->highest_occupied);
            //            rmg_printf("GW:  Highest unoccupied orbital meeting criteria = %d\n", Kptr[0]->max_unocc_res_index);

            CONVERGED = true;
        }

    }

    for(int i = 0;i < FP0_BASIS;i++) vh[i] = vh_out[i];
    delete [] vh_out; 
    return CONVERGED;
}                               /* end scf */


/******/
