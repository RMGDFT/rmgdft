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
#include "vhartree.h"
#include "common_prototypes.h"
#include "common_prototypes1.h"
#include "Functional.h"
#include "GridObject.h"



template bool Scf<double> (spinobj<double> &, fgobj<double> &, double *,
                           int, int, Kpoint<double> **, std::vector<double> &);
template bool Scf<std::complex<double> > (spinobj<double> &, fgobj<double> &, double *,
                           int, int, Kpoint<std::complex<double>> **, std::vector<double> &);

template <typename OrbitalType> bool Scf (
          spinobj<double> &vxc_in, fgobj<double> &vh_in, double *vh_ext,
          int spin_flag, int boundaryflag, Kpoint<OrbitalType> **Kptr, std::vector<double>& RMSdV)
{
    if(ct.verbose) {
        printf("PE: %d  start scf \n", pct.gridpe);
        fflush(NULL);
    }


    spinobj<double> &rho = *(Kptr[0]->rho);
    spinobj<double> &vxc = *(Kptr[0]->vxc);
    fgobj<double> &rhoc = *(Kptr[0]->rhoc);
    fgobj<double> &rhocore = *(Kptr[0]->rhocore);
    fgobj<double> &vnuc = *(Kptr[0]->vnuc);
    fgobj<double> &vh = *(Kptr[0]->vh);
    //fgobj<double> &vh_ext = *(Kptr[0]->vh_ext);

    RmgTimer RT0("2-Scf steps"), *RT1;
    RmgTimer RTt("1-TOTAL: run: Scf steps");

    bool CONVERGED = false;
    double t3;
    double *vxc_psi=NULL;
    double t[3];                  /* SCF checks and average potential */
    double max_unocc_res = 0.0;

    spinobj<double> new_rho;
    fgobj<double> vtot;
    fgobj<double> vh_dipole_corr;
    wfobj<double> vtot_psi;
    int FP0_BASIS = vtot.size();
    int P0_BASIS = vtot_psi.size();

    if(ct.dipole_corr[0]+ct.dipole_corr[0]+ct.dipole_corr[0] >0)
    {
        double dipole[3];
        get_dipole(rho.data(), dipole);
        DipoleCorrection(dipole,  vh_dipole_corr.data());
    }
    else
    {
        for (int idx = 0; idx < vtot.pbasis; idx++) {
            vh_dipole_corr[idx] = 0.0;
        }
    }
    /* save old vhxc + vnuc */
    for (int idx = 0; idx < vtot.pbasis; idx++) {
        vtot[idx] = vxc[idx] + vh[idx] + vnuc[idx] + vh_dipole_corr[idx];
    }

    if(ct.noncoll) vxc_psi = new double[4*P0_BASIS];

    /* Evaluate XC energy and potential */
    RT1 = new RmgTimer("2-Scf steps: exchange/correlation");
    Functional *F = new Functional ( *Rmg_G, Rmg_L, *Rmg_T, ct.is_gamma);
    F->v_xc(rho.data(), rhocore.data(), ct.XC, ct.vtxc, vxc.data(), ct.nspin );
    //if(pct.gridpe==0)printf("\nXC = %f  %f\n", ct.XC, ct.vtxc);
    delete F;
    delete RT1;

    double rms_target = std::min(std::max(ct.rms/ct.hartree_rms_ratio, 1.0e-12), 1.0e-6);

    RT1 = new RmgTimer("2-Scf steps: Hartree");
    double hartree_residual = VhDriver(rho.data(), rhoc.data(), vh.data(), vh_ext, rms_target);
    delete(RT1);

    /* check convergence */
    t[0] = t[1] = t[2] = 0.0;

    for (int idx = 0; idx < vtot.pbasis; idx++)
    {
        t3 = -vtot[idx];
        vtot[idx] = vxc[idx] + vh[idx] + vnuc[idx] + vh_dipole_corr[idx];
        t3 += vtot[idx];
        t[0] += rho[idx] * t3;
        t[1] += t3 * t3;
        t[2] += (vh[idx] - vh_in[idx])*(vh[idx] - vh_in[idx]);
    }                           /* idx */

    MPI_Allreduce(MPI_IN_PLACE, t, 3, MPI_DOUBLE, MPI_SUM, pct.grid_comm);
    MPI_Allreduce(MPI_IN_PLACE, t, 3, MPI_DOUBLE, MPI_SUM, pct.spin_comm);
    MPI_Allreduce(MPI_IN_PLACE, t, 3, MPI_DOUBLE, MPI_MAX, pct.img_comm);
    t[0] *= rho.vel;
    t[2] *= vh.vel;

    /* get the averaged value over each spin and each fine grid */
    if(ct.AFM)
    {
        t[1] = sqrt (t[1] / ((double) (ct.psi_fnbasis)));  
        t[2] = sqrt (t[2] / ((double) (ct.psi_fnbasis)));
    }
    else
    {
        t[1] = sqrt (t[1] / ((double) (ct.nspin * ct.psi_fnbasis)));  
        t[2] = sqrt (t[2] / ((double) (ct.nspin * ct.psi_fnbasis)));
    }
    ct.rms = t[1];
    ct.rms_vh = t[2];

    // Save input hartree potential
    vh_in = vh;

    if(std::isnan(ct.rms))
        rmg_error_handler(__FILE__, __LINE__, "NaN encountered in computational stream. Terminating.\n");

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
    GetVtotPsi (vtot_psi.data(), vtot.data(), Rmg_G->default_FG_RATIO);
    if(ct.noncoll)
        for(int is = 0; is < ct.nspin; is++)
        {
            double *vptr = vxc.data();
            GetVtotPsi (&vxc_psi[is*P0_BASIS], &vptr[is*FP0_BASIS], Rmg_G->default_FG_RATIO);
        }

    /*Generate the Dnm_I */
    get_ddd (vtot.data(), vxc.data(), true);


    // Loop over k-points
    for(int kpt = 0;kpt < ct.num_kpts_pe;kpt++) {

        if (Verify ("kohn_sham_solver","multigrid", Kptr[0]->ControlMap) || 
                ((ct.scf_steps < ct.davidson_premg) && (ct.md_steps == 0) && (ct.runflag != RESTART )) ||
                (ct.xc_is_hybrid && Functional::is_exx_active())) {
            RmgTimer *RT1 = new RmgTimer("2-Scf steps: MgridSubspace");
            Kptr[kpt]->MgridSubspace(vtot_psi.data(), vxc_psi);
            delete RT1;
        }
        else if(Verify ("kohn_sham_solver","davidson", Kptr[0]->ControlMap)) {
            int notconv;
            RmgTimer *RT1 = new RmgTimer("2-Scf steps: Davidson");
            Kptr[kpt]->Davidson(vtot_psi.data(), vxc_psi, notconv);
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


    //    if(ct.ldaU_mode != LDA_PLUS_U_NONE && (ct.ns_occ_rms >1.0e-15 || ct.scf_steps ==0) && ct.scf_steps  < ct.freeze_ldaU_steps )
    if(ct.ldaU_mode != LDA_PLUS_U_NONE && (ct.ns_occ_rms >1.0e-15 || ct.scf_steps ==0))
    {
        RmgTimer("3-MgridSubspace: ldaUop x psi");
        int pstride = Kptr[0]->ldaU->ldaU_m;
        int num_ldaU_ions = Kptr[0]->ldaU->num_ldaU_ions;
        int occ_size = ct.nspin * num_ldaU_ions * pstride * pstride;

        doubleC_4d_array old_ns_occ, new_ns_occ;
        old_ns_occ.resize(boost::extents[ct.nspin][num_ldaU_ions][Kptr[0]->ldaU->ldaU_m][Kptr[0]->ldaU->ldaU_m]);
        new_ns_occ.resize(boost::extents[ct.nspin][num_ldaU_ions][Kptr[0]->ldaU->ldaU_m][Kptr[0]->ldaU->ldaU_m]);

        old_ns_occ =  Kptr[0]->ldaU->ns_occ;

        for (int kpt =0; kpt < ct.num_kpts_pe; kpt++)
        {
            LdaplusUxpsi(Kptr[kpt], 0, Kptr[kpt]->nstates, Kptr[kpt]->orbitalsint_local);
            Kptr[kpt]->ldaU->calc_ns_occ(Kptr[kpt]->orbitalsint_local, 0, Kptr[kpt]->nstates);
            for(int idx = 0; idx< occ_size; idx++)
            {
                new_ns_occ.data()[idx] += Kptr[kpt]->ldaU->ns_occ.data()[idx];
            }
        }

        MPI_Allreduce(MPI_IN_PLACE, (double *)new_ns_occ.data(), occ_size * 2, MPI_DOUBLE, MPI_SUM, pct.kpsub_comm);

        Rmg_Symm->symm_nsocc(new_ns_occ.data(), pstride, Kptr[0]->ldaU->map_to_ldaU_ion, Kptr[0]->ldaU->ldaU_ion_index);
        if(ct.AFM)
        {
            Rmg_Symm->nsocc_AFM(new_ns_occ, Kptr[0]->ldaU->ldaU_m, Kptr[0]->ldaU->map_to_ldaU_ion, Kptr[0]->ldaU->ldaU_ion_index);
        }


        // for spin-polarized case, only the first spin on the pe are symmetrized so communicate to opposite spin
        if(ct.nspin == 2 && !ct.AFM)
        {
            MPI_Status status;
            int len = 2 * Kptr[0]->ldaU->num_ldaU_ions * pstride * pstride;
            double *sendbuf = (double *)new_ns_occ.data();
            double *recvbuf = sendbuf + len;
            MPI_Sendrecv(sendbuf, len, MPI_DOUBLE, (pct.spinpe+1)%2, pct.gridpe,
                    recvbuf, len, MPI_DOUBLE, (pct.spinpe+1)%2, pct.gridpe, pct.spin_comm, &status);

        }


        ct.ns_occ_rms = 0.0;
        for(int idx = 0; idx < occ_size; idx++)
        {
            ct.ns_occ_rms += std::norm (old_ns_occ.data()[idx] - new_ns_occ.data()[idx]);
        }
        ct.ns_occ_rms = sqrt(ct.ns_occ_rms / occ_size);

        if( ct.scf_steps  < ct.freeze_ldaU_steps )
        {
            MixLdaU(occ_size *2, (double *)new_ns_occ.data(), (double *)old_ns_occ.data(), Kptr[0]->ControlMap, false);
        }

        for (int kpt =0; kpt < ct.num_kpts_pe; kpt++)
        {
            Kptr[kpt]->ldaU->ns_occ = old_ns_occ;
            Kptr[kpt]->ldaU->calc_energy();
        }
        if(ct.verbose) Kptr[0]->ldaU->write_ldaU();

    }

    /* Take care of occupation filling */
    if(ct.verbose) {
        printf("PE: %d  start find Fermi \n", pct.gridpe);
        fflush(NULL);
    }
    ct.efermi = Fill (Kptr, ct.occ_width, ct.nel, ct.occ_mix, ct.num_states, ct.occ_flag, ct.mp_order);
    if(ct.verbose) {
        printf("PE: %d  done find Fermi \n", pct.gridpe);
        fflush(NULL);
    }

    rmg_printf ("\n");
    //progress_tag ();
    rmg_printf ("FERMI ENERGY = %15.8f eV\n", ct.efermi * Ha_eV);

    // Calculate total energy 
    // Eigenvalues are based on in potentials and density
    if(ct.verbose) {
        printf("PE: %d  start GetTe \n", pct.gridpe);
        fflush(NULL);
    }
    GetTe (rho, rhocore, rhoc, vh, vxc, Kptr, !ct.scf_steps);
    if(ct.verbose) {
        printf("PE: %d  donet GetTe \n", pct.gridpe);
        fflush(NULL);
    }

    /* Generate new density */
    RT1 = new RmgTimer("2-Scf steps: GetNewRho");
    if(ct.verbose) {
        printf("PE: %d  start GetNewRho \n", pct.gridpe);
        fflush(NULL);
    }
    GetNewRho(Kptr, new_rho.data());
    if(ct.verbose) {
        printf("PE: %d  done GetNewRho \n", pct.gridpe);
        fflush(NULL);
    }
    new_rho.get_oppo();
    delete(RT1);

    // Get Hartree potential for the output density
    fgobj<double> vh_out;
    RT1 = new RmgTimer("2-Scf steps: Hartree");
    VhDriver(new_rho.data(), rhoc.data(), vh_out.data(), vh_ext, rms_target);
    delete RT1;

    // Compute convergence measure (2nd order variational term) and average by nspin
    double sum = 0.0;
    for(int i = 0;i < rho.pbasis;i++) sum += (vh_out[i] - vh[i]) * (new_rho[i] - rho[i]);
    if(ct.AFM) 
    {
        for(int i = 0;i < rho.pbasis;i++) sum += (vh_out[i] - vh[i]) * (new_rho.dw[i] - rho.dw[i]);
    }
    sum = 0.5 * rho.vel * sum;

    MPI_Allreduce(MPI_IN_PLACE, &sum, 1, MPI_DOUBLE, MPI_SUM, pct.grid_comm);
    MPI_Allreduce(MPI_IN_PLACE, &sum, 1, MPI_DOUBLE, MPI_SUM, pct.spin_comm);
    MPI_Allreduce(MPI_IN_PLACE, &sum, 1, MPI_DOUBLE, MPI_MAX, pct.img_comm);
    ct.scf_accuracy = sum;

    // Compute variational energy correction term if any
    sum = EnergyCorrection(Kptr, rho.data(), new_rho.data(), vh.data(), vh_out.data());
    ct.scf_correction = sum;

    // Check if this convergence threshold has been reached
    if(!Verify ("freeze_occupied", true, Kptr[0]->ControlMap)) {
        if (ct.scf_steps && fabs(ct.scf_accuracy) < ct.adaptive_thr_energy) CONVERGED = true;
    }

    if( (CONVERGED || (ct.scf_steps == (ct.max_scf_steps-1))) )
    {
        // If the multigrid solver is selected the total energy calculation from
        // the loop above is not variational but the following block of code
        // will give us a variational energy.
        if (Verify ("kohn_sham_solver","multigrid", Kptr[0]->ControlMap) && !ct.noncoll 
                && ct.potential_acceleration_constant_step > 1.0e-5)
        {
            ct.scf_correction = 0.0;
            for (int idx = 0; idx < vtot.pbasis; idx++)
            {
                vtot[idx] = vxc[idx] + vh[idx] + vnuc[idx];
            }
            // Transfer vtot from the fine grid to the wavefunction grid
            GetVtotPsi (vtot_psi.data(), vtot.data(), Rmg_G->default_FG_RATIO);

            /*Generate the Dnm_I */
            get_ddd (vtot.data(), vxc.data(), true);

            // We save ct.potential_acceleration_constant_step and set it to zero
            // when we do the diagonalization since otherwise the energy is not variational
            double potential_acceleration_step = ct.potential_acceleration_constant_step;
            ct.potential_acceleration_constant_step = 0.0;
            for(int kpt = 0;kpt < ct.num_kpts_pe;kpt++)
            {
                Kptr[kpt]->Subdiag(vtot_psi.data(), vxc_psi, ct.subdiag_driver);
                Kptr[kpt]->BetaProjector->project(Kptr[kpt], Kptr[kpt]->newsint_local, 0, 
                        Kptr[kpt]->nstates*ct.noncoll_factor, 
                        Kptr[kpt]->nl_weight);
                MPI_Barrier(pct.grid_comm);
            }

            if(ct.ldaU_mode != LDA_PLUS_U_NONE)
            {
                int pstride = Kptr[0]->ldaU->ldaU_m;
                int num_ldaU_ions = Kptr[0]->ldaU->num_ldaU_ions;
                int occ_size = ct.nspin * num_ldaU_ions * pstride * pstride;

                doubleC_4d_array old_ns_occ, new_ns_occ;
                new_ns_occ.resize(boost::extents[ct.nspin][num_ldaU_ions][Kptr[0]->ldaU->ldaU_m][Kptr[0]->ldaU->ldaU_m]);

                for (int kpt =0; kpt < ct.num_kpts_pe; kpt++)
                {
                    LdaplusUxpsi(Kptr[kpt], 0, Kptr[kpt]->nstates, Kptr[kpt]->orbitalsint_local);
                    Kptr[kpt]->ldaU->calc_ns_occ(Kptr[kpt]->orbitalsint_local, 0, Kptr[kpt]->nstates);
                    for(int idx = 0; idx< occ_size; idx++)
                    {
                        new_ns_occ.data()[idx] += Kptr[kpt]->ldaU->ns_occ.data()[idx];
                    }
                }

                MPI_Allreduce(MPI_IN_PLACE, (double *)new_ns_occ.data(), occ_size * 2, MPI_DOUBLE, MPI_SUM, pct.kpsub_comm);

                Rmg_Symm->symm_nsocc(new_ns_occ.data(), pstride, Kptr[0]->ldaU->map_to_ldaU_ion, Kptr[0]->ldaU->ldaU_ion_index);
                if(ct.AFM)
                {
                    Rmg_Symm->nsocc_AFM(new_ns_occ, Kptr[0]->ldaU->ldaU_m, Kptr[0]->ldaU->map_to_ldaU_ion, Kptr[0]->ldaU->ldaU_ion_index);
                }


                // for spin-polarized case, only the first spin on the pe are symmetrized so communicate to opposite spin
                if(ct.nspin == 2 && !ct.AFM)
                {
                    MPI_Status status;
                    int len = 2 * Kptr[0]->ldaU->num_ldaU_ions * pstride * pstride;
                    double *sendbuf = (double *)new_ns_occ.data();
                    double *recvbuf = sendbuf + len;
                    MPI_Sendrecv(sendbuf, len, MPI_DOUBLE, (pct.spinpe+1)%2, pct.gridpe,
                            recvbuf, len, MPI_DOUBLE, (pct.spinpe+1)%2, pct.gridpe, pct.spin_comm, &status);

                }


                for (int kpt =0; kpt < ct.num_kpts_pe; kpt++)
                {
                    Kptr[kpt]->ldaU->ns_occ = new_ns_occ;
                    Kptr[kpt]->ldaU->calc_energy();
                }
                if(ct.verbose) Kptr[0]->ldaU->write_ldaU();
            }


            if (ct.nspin == 2) GetOppositeEigvals (Kptr);
            GetTe (rho, rhocore, rhoc, vh, vxc, Kptr, !ct.scf_steps);
            ct.potential_acceleration_constant_step = potential_acceleration_step;
        }

        // Evaluate XC energy and potential from the output density
        // for the force correction
        RT1 = new RmgTimer("2-Scf steps: exchange/correlation");
        vxc_in = vxc;
        Functional *F = new Functional ( *Rmg_G, Rmg_L, *Rmg_T, ct.is_gamma);
        F->v_xc(new_rho.data(), rhocore.data(), ct.XC, ct.vtxc, vxc.data(), ct.nspin );
        delete F;
        delete RT1;

    }


    /*Takes care of mixing and checks whether the charge density is negative*/
    RT1 = new RmgTimer("2-Scf steps: MixRho");
    MixRho(new_rho.data(), rho.data(), rhocore.data(), vh.data(), vh_out.data(), rhoc.data(), Kptr[0]->ControlMap, false);
    delete RT1;

    rho.get_oppo();


    /* Make sure we see output, e.g. so we can shut down errant runs */
    fflush( ct.logfile );
#if !(__CYGWIN__ || __MINGW64__ || _WIN32 || _WIN64)
    fsync( fileno(ct.logfile) );
#endif

    /* release memory */
    if(ct.noncoll) delete [] vxc_psi;

    if(Verify ("freeze_occupied", true, Kptr[0]->ControlMap)) {

        if(ct.scf_steps && (max_unocc_res < ct.gw_threshold)) {
            rmg_printf("\nGW: convergence criteria of %10.5e has been met.\n", ct.gw_threshold);
            rmg_printf("GW:  Highest occupied orbital index              = %d\n", Kptr[0]->highest_occupied);
            //            rmg_printf("GW:  Highest unoccupied orbital meeting criteria = %d\n", Kptr[0]->max_unocc_res_index);

            CONVERGED = true;
        }

    }

    vh = vh_out;
    if(ct.verbose) {
        printf("PE: %d  done scf \n", pct.gridpe);
        fflush(NULL);
    }
    return CONVERGED;
}                               /* end scf */

