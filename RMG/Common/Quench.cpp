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


#include <sys/stat.h>
#include <sys/types.h>
#include <fcntl.h>
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
#include "RmgException.h"
#include "Kpoint.h"
#include "transition.h"
#include "Plots.h"
#include "Functional.h"
#include "Exxbase.h"
#include "../Headers/macros.h"
#include "Stress.h"
#include "Voronoi.h"
#include "GpuAlloc.h"
#include "Wannier.h"
#include "RmgSumAll.h"
#include "FDOpt.h"
#include "BerryPhase.h"




// Local function prototypes
void PlotConvergence(std::vector<double> &RMSdV, bool CONVERGED);
void ChargeAnalysis(spinobj<double> &rho, std::unordered_map<std::string, InputKey *>& ControlMap, Voronoi &);


// Instantiate gamma and non-gamma versions
template bool Quench<double> (Kpoint<double> **Kptr, bool);
template bool Quench<std::complex<double> > (Kpoint<std::complex<double>> **Kptr, bool);

template <typename OrbitalType> bool Quench (Kpoint<OrbitalType> **Kptr, bool compute_forces)
{
    if(ct.verbose) {
        printf("PE: %d  start Quench \n", pct.gridpe);
        fflush(NULL);
    }

    /* reset total number of scf steps */
    ct.scf_converged = false;
    ct.scf_is_converging = true;

    spinobj<double> &rho = *(Kptr[0]->rho);
    spinobj<double> &vxc = *(Kptr[0]->vxc);
    fgobj<double> &rhoc = *(Kptr[0]->rhoc);
    fgobj<double> &rhocore = *(Kptr[0]->rhocore);
    fgobj<double> &vnuc = *(Kptr[0]->vnuc);
    fgobj<double> &vh = *(Kptr[0]->vh);

    ct.adaptive_thr_energy = ct.thr_energy;
    static std::vector<double> RMSdV;
    static std::vector<double> etot;
    static std::vector<double> dStress;

    Functional *F = new Functional ( *Rmg_G, Rmg_L, *Rmg_T, ct.is_gamma);
    std::string tempwave = std::string(ct.outfile) + "_serial";
    int FP0_BASIS =  Rmg_G->get_P0_BASIS(Rmg_G->get_default_FG_RATIO());
    fgobj<double> vh_in;
    spinobj<double> vxc_in;

    // Now we save the initial wavefunctions in case we need to restart
    std::string initfile(ct.outfile);
    initfile = initfile + std::string("_init");
    if (!(!strcmp ("/dev/null", ct.outfile)) || (!strcmp ("/dev/null/", ct.outfile)) )
        WriteRestart ((char *)initfile.c_str(), vh.data(), rho.up.data(), rho.dw.data(), vxc.data(), Kptr);


    /* ---------- begin scf loop ---------- */

    double start_time = my_crtc ();
    double exx_start_time = start_time;
    double step_time;
    double exx_step_time=0.0, exx_elapsed_time=0.0;
    double f0=0.0,f1,f2=0.0,exxen=0.0;

    RmgTimer *RT = new RmgTimer("Init Voronoi");
    Voronoi *Voronoi_charge = new Voronoi();
    delete RT;
    int outer_steps = 1;

    Exxbase<OrbitalType> *Exx_scf = NULL;

    ct.FOCK = 0.0;
    if(ct.xc_is_hybrid)
    {
        outer_steps = ct.max_exx_steps;

        std::vector<double> occs;
        occs.resize(Kptr[0]->nstates);
        for(int i=0;i < Kptr[0]->nstates;i++) occs[i] = Kptr[0]->Kstates[i].occupation[0];

        Exx_scf = new Exxbase<OrbitalType>(*Kptr[0]->G, *Rmg_halfgrid, *Kptr[0]->L, tempwave, Kptr[0]->nstates, occs.data(),
                Kptr[0]->orbital_storage, ct.exx_mode);

        occs.clear();
        if(ct.runflag == RESTART)
        {
            F->start_exx_rmg();
            Exx_scf->Vexx(Kptr[0]->vexx, false);
            ct.FOCK = Exx_scf->Exxenergy(Kptr[0]->vexx);
            exxen = ct.FOCK;
        }
    }

    ct.exx_delta = DBL_MAX;
    ct.vexx_rms = DBL_MAX;
    for(ct.exx_steps = 0;ct.exx_steps < outer_steps;ct.exx_steps++)
    { 
        exx_step_time = my_crtc ();
        RMSdV.clear();
        dStress.clear();

        ct.total_scf_steps = 0;
        // Adjust exx convergence threshold
        if(ct.xc_is_hybrid)
        {
            Exx_scf->vexx_RMS[ct.exx_steps] = 0.0;
            if(ct.exx_steps)
                ct.adaptive_thr_energy = std::max(1.0e-15,std::min(1.0e-8, ct.vexx_rms / 1000.0));
            else
                ct.adaptive_thr_energy = 1.0e-8;
        }

        for (ct.scf_steps = 0, ct.scf_converged = false;
                ct.total_scf_steps < ct.max_scf_steps && !ct.scf_converged; ct.scf_steps++, ct.total_scf_steps++)
        {

            /* perform a single self-consistent step */
            step_time = my_crtc ();
            if(ct.forceflag == NSCF)
            {
                ct.scf_converged = Nscf (vxc.data(), vxc_in.data(), vh.data(), vh_in.data(), ct.vh_ext, vnuc.data(), rho.up.data(), rho.dw.data(), rhocore.data(), rhoc.data(), ct.spin_flag, ct.boundaryflag, Kptr, RMSdV);
            }
            else
            {
                double tmix = AutoMix();
                if(tmix != ct.mix && ct.mix > 0.02)
                {
                    bool reset = true;
                    if(ct.scf_accuracy > 1.0e-3)
                    {
                        std::string initfile(ct.outfile);
                        initfile = initfile + std::string("_init");
                        ReadData ((char *)initfile.c_str(), vh.data(), rho.data(), vxc.data(), Kptr);
                        rho.get_oppo();
                    }
                    ct.mix /= 2.0;
                    MixRho(NULL, NULL, NULL, NULL, NULL, NULL, Kptr[0]->ControlMap, reset);
                    RMSdV.clear();
                    dStress.clear();
                    ct.scf_steps = 0;
                    ct.exx_steps = 0;
                }
                if(ct.mix < 0.02)
                {
                    rmg_printf ("\nCharge density mixing has fallen below  %10.4f.\n", ct.mix);
                    rmg_printf ("and usually means something is wrong with the job.\n");
                }

                // Should add a check here for a minimum density mixing to trigger an error
                ct.scf_converged = Scf (vxc_in, vh_in, ct.vh_ext, ct.spin_flag, ct.boundaryflag, Kptr, RMSdV);
            }
            step_time = my_crtc () - step_time;

            // Save data to file for future restart at checkpoint interval if this is a quench run.
            // For Relaxation and molecular dynamics we save at the end of each ionic step.
            if (ct.checkpoint && (ct.scf_steps % ct.checkpoint == 0) )
                WriteRestart (ct.outfile, vh.data(), rho.up.data(), rho.dw.data(), vxc.data(), Kptr);

            /* output the eigenvalues with occupations */
            if (ct.write_eigvals_period && (ct.scf_steps % ct.write_eigvals_period == 0))
            {
                OutputEigenvalues (Kptr, 0, ct.scf_steps);
                rmg_printf ("\nTotal charge in supercell = %16.8f\n", ct.tcharge);
            }

            /*Perform charge analysis if requested*/
            ChargeAnalysis(rho, Kptr[0]->ControlMap, *Voronoi_charge);

#if PLPLOT_LIBS
            // Generate convergence plots
            PlotConvergence(RMSdV, ct.scf_converged);
#endif

            // Write out progress info
            double elapsed_time = my_crtc() - start_time;
            ProgressTag(step_time, elapsed_time);

        }
        /* ---------- end scf loop ---------- */

        if(ct.BerryPhase)
        {
            Rmg_BP->CalcBP(Kptr);
        }

        etot.push_back(ct.TOTAL);

        // If a hybrid calculation compute vexx
        if(ct.xc_is_hybrid)
        {
            F->start_exx_rmg();
            // f1 is fock energy calculated using vexx from previous orbitals and the current orbitals
            f1 = Exx_scf->Exxenergy(Kptr[0]->vexx);
            Exx_scf->Vexx(Kptr[0]->vexx, (fabs(ct.exx_delta) > ct.vexx_fft_threshold));
            // f2 is fock energy calculated using vexx from current orbitals and the current orbitals
            f2 = Exx_scf->Exxenergy(Kptr[0]->vexx);
            ct.exx_delta = f1 - 0.5*(f2 + f0);
            f0 = f2;
            ct.FOCK = exxen + f2 - f1;
            exxen = f2;

            exx_step_time = my_crtc () - exx_step_time;
            exx_elapsed_time = my_crtc() - exx_start_time;

            double deltaE = 1.0;
            if(ct.exx_steps > 0) deltaE = std::abs(etot[ct.exx_steps] - etot[ct.exx_steps-1]);
            if(Exx_scf->vexx_RMS[ct.exx_steps] < ct.exx_convergence_criterion || deltaE < 1.0e-7)
            { 
                rmg_printf(" Finished EXX outer loop in %3d exx steps, elapsed time = %6.2f, vexx_rms = %8.2e, total energy = %.*f Ha\n",
                        ct.exx_steps, exx_elapsed_time, Exx_scf->vexx_RMS[ct.exx_steps], 6, ct.TOTAL);
                ct.FOCK = f2;
                break;
            }
            else
            {
                rmg_printf(" Finished EXX inner loop in %3d scf steps, exx step time = %6.2f, vexx_rms = %8.2e, total energy = %.*f Ha\n",
                        ct.scf_steps, exx_step_time, Exx_scf->vexx_RMS[ct.exx_steps], 6, ct.TOTAL);
            }
        }

    }
    /* ---------- end exx loop ---------- */


    if(ct.scf_converged)
    {
        if(ct.rms < ct.thr_rms)
        {
            rmg_printf("\n Convergence criterion reached: potential RMS (%.2e) is lower than threshold (%.2e)\n", ct.rms, ct.thr_rms);
            if (pct.imgpe == 0 && pct.images == 1)
                fprintf(stdout,"\n Convergence criterion reached: potential RMS (%.2e) is lower than threshold (%.2e)\n", ct.rms, ct.thr_rms);
        }
        else if(fabs(ct.scf_accuracy) < ct.thr_energy)
        {
            rmg_printf("\n Convergence criterion reached: Energy variation (%.2e) is lower than threshold (%.2e)\n", fabs(ct.scf_accuracy), ct.thr_energy);
            if (pct.imgpe == 0 && pct.images ==1)
                fprintf(stdout, "\n Convergence criterion reached: Energy variation (%.2e) is lower than threshold (%.2e)", fabs(ct.scf_accuracy), ct.thr_energy);
        }

        rmg_printf ("\n");
        rmg_printf ("potential convergence has been achieved. stopping ...\n");

        /*Write PDOS if converged*/
        //	if (ct.pdos_flag)
        //	    get_pdos (Kptr[0]->kstates, ct.Emin, ct.Emax, ct.E_POINTS);

    }
    else
    {
        rmg_printf("\n Convergence criterion not met but max_scf_steps %d was reached.\n", ct.max_scf_steps);
        if (pct.imgpe == 0 && pct.images ==1)
            fprintf(stdout, "\n Convergence criterion not met but max_scf_steps %d was reached.\n", ct.max_scf_steps);
    }


    rmg_printf ("\n");
    progress_tag ();

    if(0)
    {
        if(ct.is_gamma)
        {
            FDOpt *Fopt = new FDOpt();

            double *one_orbital = (double *)Kptr[0]->Kstates[0].psi;
            Fopt->Analyze_fft(one_orbital);
            delete Fopt;
        }

    }



    double *v_psi, *vxc_psi;
    int pbasis = Kptr[0]->pbasis;
    v_psi = new double[pbasis];
    vxc_psi = new double[pbasis]();
    int nstates = Kptr[0]->nstates;
    OrbitalType *Hcore = (OrbitalType *)RmgMallocHost(ct.num_kpts_pe * nstates * nstates * sizeof(OrbitalType));
    OrbitalType *Hcore_kin = (OrbitalType *)RmgMallocHost(ct.num_kpts_pe * nstates * nstates * sizeof(OrbitalType));
    OrbitalType *Hcore_localpp = (OrbitalType *)RmgMallocHost(ct.num_kpts_pe * nstates * nstates * sizeof(OrbitalType));

    bool compute_direct = (ct.write_qmcpack_restart ||
            ct.write_qmcpack_restart_localized) && ct.norm_conserving_pp;

    double efactor = ct.energy_output_conversion[ct.energy_output_units];
    const char *eunits = ct.energy_output_string[ct.energy_output_units].c_str();
    rmg_printf ("\nfinal total energy from eig sum = %16.8f %s\n", efactor*ct.TOTAL, eunits);

    // test conditions
    check_tests();

    if(compute_direct)
    {
        double kin_energy=0.0, pseudo_energy= 0.0, total_e = 0.0, E_localpp = 0.0, E_nonlocalpp;
        compute_vxc(rho, rhocore, ct.XC, ct.vtxc, vxc, ct.nspin );

        GetNewRho(Kptr, rho.data());
        rho.get_oppo();

        VhDriver(rho.data(), rhoc.data(), vh.data(), ct.vh_ext, 1.0e-12);

        GetTe (rho, rhocore, rhoc, vh, vxc, Kptr, true);

        GetVtotPsi (v_psi, vnuc.data(), Rmg_G->default_FG_RATIO);
        for(int kpt = 0; kpt < ct.num_kpts_pe; kpt++) {
            Kptr[kpt]->ComputeHcore(v_psi, vxc_psi, &Hcore[kpt*nstates * nstates], &Hcore_kin[kpt*nstates * nstates], &Hcore_localpp[kpt * nstates * nstates]);

            for (int st = 0; st < nstates; st++)
            {
                double occ = Kptr[kpt]->Kstates[st].occupation[0] * Kptr[kpt]->kp.kweight;
                kin_energy += std::real(Hcore_kin[kpt * nstates * nstates + st * nstates + st]) * occ;
                pseudo_energy += std::real(Hcore[kpt * nstates * nstates + st * nstates + st]) * occ;
                E_localpp += std::real(Hcore_localpp[kpt * nstates * nstates + st * nstates + st]) * occ;
            }
        }
        kin_energy = RmgSumAll(kin_energy, pct.kpsub_comm);
        kin_energy = RmgSumAll(kin_energy, pct.spin_comm);
        pseudo_energy = RmgSumAll(pseudo_energy, pct.kpsub_comm);
        pseudo_energy = RmgSumAll(pseudo_energy, pct.spin_comm);
        E_localpp = RmgSumAll(E_localpp, pct.kpsub_comm);
        E_localpp = RmgSumAll(E_localpp, pct.spin_comm);

        pseudo_energy -= kin_energy;
        E_nonlocalpp = pseudo_energy - E_localpp;


        //Kptr[0]->ldaU->Ecorrect
        total_e = kin_energy + pseudo_energy + ct.ES + ct.XC + ct.II + ct.Evdw + ct.ldaU_E + ct.FOCK;
        /* Print contributions to total energies into output file */
        std::vector<double> occs;
        occs.resize(Kptr[0]->nstates);
        for(int i=0;i < Kptr[0]->nstates;i++) occs[i] = Kptr[0]->Kstates[i].occupation[0];
        Exxbase<OrbitalType> *Exx = NULL;
        Exx = new Exxbase<OrbitalType>(*Kptr[0]->G, *Rmg_halfgrid, *Kptr[0]->L, tempwave, Kptr[0]->nstates, occs.data(),
                Kptr[0]->orbital_storage, ct.exx_mode);
        if(ct.exx_mode == EXX_LOCAL_FFT)
            Exx->WriteWfsToSingleFile();
        if(ct.qmc_nband == 0) ct.qmc_nband = ct.num_states;
        if(ct.qmc_nband > ct.num_states)
            throw RmgFatalException() << "qmc_nband " << ct.qmc_nband << " is larger than ct.num_states " << ct.num_states << "\n";

        Exx->SetHcore(Hcore, Hcore_kin, nstates);
        if(ct.exx_int_flag)
        {
            Exx->Vexx_integrals(ct.exx_int_file);
        }
        ct.exx_fraction = 1.0;
        ct.exx_steps = 0;
        Exx->Vexx(Kptr[0]->vexx, false);
        double E_exchange = Exx->Exxenergy(Kptr[0]->vexx);

        rmg_printf ("\n@@ TOTAL ENEGY Components \n");
        rmg_printf ("@@ ION_ION            = %15.6f %s\n", efactor*ct.II, eunits);
        rmg_printf ("@@ ELECTROSTATIC      = %15.6f %s\n", efactor*ct.ES, eunits);
        rmg_printf ("@@ EXC                = %15.6f %s\n", efactor*ct.XC, eunits);
        rmg_printf ("@@ Kinetic            = %15.6f %s\n", efactor*kin_energy, eunits);
        rmg_printf ("@@ E_localpp          = %15.6f %s\n", efactor*E_localpp, eunits);
        rmg_printf ("@@ E_nonlocalpp       = %15.6f %s\n", efactor*E_nonlocalpp, eunits);
        if(ct.vdw_corr)
            rmg_printf ("@@ vdw correction     = %15.6f %s\n", efactor*ct.Evdw, eunits);
        if((ct.ldaU_mode != LDA_PLUS_U_NONE) && (ct.num_ldaU_ions > 0))
            rmg_printf ("@@ LdaU correction    = %15.6f %s\n", efactor*ct.ldaU_E, eunits);
        rmg_printf ("final total energy from direct =  %16.8f %s\n", efactor*total_e, eunits);
        ct.TOTAL_DIRECT = total_e;

        double Madelung = MadelungConstant();
        rmg_printf ("MadelungConstant      =  %16.8f \n", Madelung);

        double vme = 0.5 * ct.nel * Madelung;
        total_e = kin_energy + pseudo_energy + ct.ES + E_exchange + ct.II + ct.Evdw + ct.ldaU_E;
        total_e += vme;
        rmg_printf ("\n Hartree Fock total energy \n");
        rmg_printf ("@@ ION_ION            = %15.6f %s\n", efactor*ct.II, eunits);
        rmg_printf ("@@ ELECTROSTATIC      = %15.6f %s\n", efactor*ct.ES, eunits);
        rmg_printf ("@@ Exchange           = %15.6f %s\n", efactor*E_exchange, eunits);
        rmg_printf ("@@ Kinetic            = %15.6f %s\n", efactor*kin_energy, eunits);
        rmg_printf ("@@ E_localpp          = %15.6f %s\n", efactor*E_localpp, eunits);
        rmg_printf ("@@ E_nonlocalpp       = %15.6f %s\n", efactor*E_nonlocalpp, eunits);
        rmg_printf ("@@ Madelung           = %15.6f %s\n", efactor*vme, eunits);
        if(ct.vdw_corr)
            rmg_printf ("@@ vdw correction     = %15.6f %s\n", efactor*ct.Evdw, eunits);
        if((ct.ldaU_mode != LDA_PLUS_U_NONE) && (ct.num_ldaU_ions > 0))
            rmg_printf ("@@ LdaU correction    = %15.6f %s\n", efactor*ct.ldaU_E, eunits);
        if(ct.BerryPhase) 
        {
            total_e += Rmg_BP->enthalpy_elec;
            rmg_printf ("@@ Electric Enthalpy  = %15.6f %s\n", efactor*Rmg_BP->enthalpy_elec, eunits);
        }

        rmg_printf ("total energy Hartree Fock(+Madelung)      =  %16.8f %s\n", efactor*total_e, eunits);
        rmg_printf ("\n WARNING: Madelung term should not be included, add it to compare with qmcpack\n" );


        fflush(NULL);
        delete Exx;

    }


    RmgFreeHost(Hcore);
    RmgFreeHost(Hcore_kin);
    RmgFreeHost(Hcore_localpp);
    delete [] v_psi;
    delete [] vxc_psi;


    /* output final eigenvalues with occupations */
    OutputEigenvalues (Kptr, 0, ct.scf_steps);
    OutputDos(Kptr);
    rmg_printf ("\nTotal charge in supercell = %16.8f\n", ct.tcharge);

    // Output RMSdV for convergence analysis
    if(pct.imgpe == 0) {
        // Write out convergence file
        std::string ConvergenceFile(ct.basename);
        ConvergenceFile = ConvergenceFile + ".rmsdv.xmgr";
        mode_t mode = O_CREAT |O_TRUNC |O_RDWR;
        //  if(ct.md_steps > 0) mode = O_RDWR | O_APPEND;

        int fhand = open(ConvergenceFile.c_str(), mode, S_IREAD | S_IWRITE);
        if (fhand < 0)
            throw RmgFatalException() <<  "Unable to write file in " << __FILE__ << " at line " << __LINE__ << "\n";
        char tbuf[100];

        int idx = 0;
        snprintf(tbuf, sizeof(tbuf), "@    title \"RMG convergence plot\"\n@    xaxis  label \"SCF steps\"@    yaxis  label \"log10(RMSdV)\"\n");

        for(auto it = RMSdV.begin();it != RMSdV.end();it++) {
            snprintf(tbuf, sizeof(tbuf), "%d %12.6f\n", idx, log10(*it));
            write(fhand, tbuf, strlen(tbuf));
            idx++;
        }

        snprintf(tbuf, sizeof(tbuf), "&\n");

    }

    if(ct.stress)
    {
        bool local_only = false;
        fgobj<double> vtot;
        for(int idx = 0; idx < FP0_BASIS; idx++) vtot[idx] = vh[idx] + vnuc[idx] + vxc[idx];
        Stress<OrbitalType> Stress_cal(Kptr, Rmg_L, *Rmg_G, *fine_pwaves, Atoms, Species, 
                ct.XC, vxc, rho, rhocore, vtot, ct.stress_tensor, local_only);
    }


    /*When running MD, force pointers need to be rotated before calculating new forces */
    if(compute_forces)
    {
        ct.fpt[0] = 0;
        ct.fpt[1] = 1;
        ct.fpt[2] = 2;
        ct.fpt[3] = 3;
        ct.sqrt_interpolation = false;

        spinobj<double> trho;
        trho = rho;
        GetNewRho(Kptr, trho.data());

        for (size_t ion = 0, i_end = Atoms.size(); ion < i_end; ++ion)
        {
            Atoms[ion].RotateForces();
        }
        Force (trho.up.data(), trho.dw.data(), rhoc.data(), vh.data(), vh_in.data(), vxc.data(), vxc_in.data(), vnuc.data(), Kptr);
        fgobj<double> vtot;
        for(int idx = 0; idx < vtot.pbasis; idx++)
        {
            vtot[idx] = vh[idx] + vxc[idx] + vnuc[idx];
        }

        get_ddd (vtot.data(), vxc.data(), true);

    }


    rmg_printf (" volume and energy per atom = %18.8f  %18.8f eV\n", Rmg_L.get_omega()*a0_A*a0_A*a0_A/Atoms.size(),ct.TOTAL * Ha_eV/Atoms.size());


    /*Calculate and write dipole moment if requested*/
    if (ct.dipole_moment)
    {
        double dipole[3];
        get_dipole(rho.data(), dipole);

        double DEBYE_CONVERSION = 2.54174618772314;

        /*Now we need to convert to debye units */
        if (pct.imgpe==0)
        {
            rmg_printf("\n\n Dipole moment [Debye]: (%16.8e,%16.8e, %16.8e)", 
                    DEBYE_CONVERSION *dipole[0], 
                    DEBYE_CONVERSION *dipole[1], 
                    DEBYE_CONVERSION *dipole[2]);
        }

    }

    WriteRestart (ct.outfile, vh.data(), rho.up.data(), rho.dw.data(), vxc.data(), Kptr);
    if(ct.wannier90)
    {
        int scdm = ct.wannier90_scdm;
        double scdm_mu = ct.wannier90_scdm_mu;
        double scdm_sigma = ct.wannier90_scdm_sigma;
        int n_wannier = ct.num_wanniers;

        Wannier<OrbitalType> Wan(*Kptr[0]->G, *Kptr[0]->L, "WfsForWannier90/wfs", Kptr[0]->nstates, 
                n_wannier, scdm, scdm_mu, scdm_sigma, Kptr[0]->orbital_storage, Kptr);
        for(int kpt = 0; kpt < ct.num_kpts_pe; kpt++) {
            int kpt_global = kpt + pct.kstart;
            Write_Wfs_forWannier(kpt_global, Kptr[kpt], Wan.exclude_bands, "WfsForWannier90/wfs");
        }
        Wan.AmnMmn("WfsForWannier90/wfs");

    }

    /* output the forces */
    if (pct.imgpe == 0)
        write_force ();

    if(ct.verbose) {
        printf("PE: %d  done Quench \n", pct.gridpe);
        fflush(NULL);
    }

    return ct.scf_converged;

}                               /* end quench */


#if PLPLOT_LIBS
void PlotConvergence(std::vector<double> &RMSdV, bool ct.scf_converged)
{
    if(pct.imgpe == 0) {
        std::vector<double> x;
        std::string ConvergencePlot(ct.basename);
        ConvergencePlot = ConvergencePlot + ".rmsdv.png";
        std::string ConvergenceTitle("RMG convergence:");
        if(ct.scf_converged) {
            ConvergenceTitle = ConvergenceTitle + " quench completed.";
        }
        else {
            ConvergenceTitle = ConvergenceTitle + " quenching electrons.";
        }
        LinePlotLog10y(ConvergencePlot.c_str(), "SCF Steps", "log10(RMS[dV])", ConvergenceTitle.c_str(), x, RMSdV);
    }
}
#endif


void ChargeAnalysis(spinobj<double> &rho, std::unordered_map<std::string, InputKey *>& ControlMap, Voronoi &Voronoi_charge)
{
    /*Perform charge analysis if requested*/
    if (ct.charge_analysis_period)
    {
        if (ct.scf_steps % ct.charge_analysis_period == 0)
        {
            if (Verify("charge_analysis","Voronoi", ControlMap))
            {
                double timex = my_crtc ();
                double *localrho = new double[Atoms.size()];
                if(ct.nspin == 1)
                {
                    Voronoi_charge.LocalCharge(rho.data(), localrho);
                    for(size_t ion = 0; ion < Atoms.size(); ion++)
                        Atoms[ion].partial_charge = Voronoi_charge.localrho_atomic[ion] - localrho[ion];
                }

                if(ct.nspin == 2)
                {
                    double *localrho_up = new double[Atoms.size()];
                    double *localrho_dn = new double[Atoms.size()];
                    Voronoi_charge.LocalCharge(rho.up.data(), localrho_up);
                    Voronoi_charge.LocalCharge(rho.dw.data(), localrho_dn);

                    for (size_t ion = 0, i_end = Atoms.size(); ion < i_end; ++ion)
                    {
                        Atoms[ion].mag = localrho_up[ion]-localrho_dn[ion];
                    }
                    for(size_t ion = 0; ion < Atoms.size(); ion++)
                        Atoms[ion].partial_charge = Voronoi_charge.localrho_atomic[ion] - localrho_up[ion] - localrho_dn[ion];
                    delete [] localrho_up;
                    delete [] localrho_dn;
                }

                if(ct.noncoll)
                {
                    Voronoi_charge.LocalCharge(rho.data(), localrho);
                    for(size_t ion = 0; ion < Atoms.size(); ion++)
                        Atoms[ion].partial_charge = Voronoi_charge.localrho_atomic[ion] - localrho[ion];
                    double *localrho_x = new double[Atoms.size()];
                    double *localrho_y = new double[Atoms.size()];
                    double *localrho_z = new double[Atoms.size()];

                    Voronoi_charge.LocalCharge(rho.cx.data(), localrho_x);
                    Voronoi_charge.LocalCharge(rho.cy.data(), localrho_y);
                    Voronoi_charge.LocalCharge(rho.cz.data(), localrho_z);
                    rmg_printf("\n@ION  Ion  Species      Magnetization_xyz(Voronoi)\n");
                    for (size_t ion = 0, i_end = Atoms.size(); ion < i_end; ++ion)
                    {
                        ION &Atom = Atoms[ion];
                        SPECIES &AtomType = Species[Atom.species];

                        rmg_printf ("@ION  %3lu  %4s        %7.3f  %7.3f  %7.3f\n",
                                ion + 1, AtomType.atomic_symbol,localrho_x[ion], localrho_y[ion], localrho_z[ion]);
                    }

                    delete [] localrho_x;
                    delete [] localrho_y;
                    delete [] localrho_z;
                }
                delete [] localrho;

                rmg_printf("\n Vdd took %f seconds\n", my_crtc () - timex);
            }
            WriteChargeAnalysis();
        }
    }
}

