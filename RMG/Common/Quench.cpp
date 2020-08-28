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


// Local function prototypes
void PlotConvergence(std::vector<double> &RMSdV, bool CONVERGED);
void ChargeAnalysis(double *rho, std::unordered_map<std::string, InputKey *>& ControlMap, Voronoi &);


// Instantiate gamma and non-gamma versions
template bool Quench<double> (double *, double *, double *, double *, double *, double *, double *, Kpoint<double> **Kptr, bool);
template bool Quench<std::complex<double> > (double *, double *, double *, double *, double *, double *, double *, Kpoint<std::complex<double>> **Kptr, bool);

template <typename OrbitalType> bool Quench (double * vxc, double * vh, double * vnuc, double * rho,
             double * rho_oppo, double * rhocore, double * rhoc, Kpoint<OrbitalType> **Kptr, bool compute_forces)
{

    bool CONVERGED = false;
    ct.adaptive_thr_energy = ct.thr_energy;
    static std::vector<double> RMSdV;
    Functional *F = new Functional ( *Rmg_G, Rmg_L, *Rmg_T, ct.is_gamma);

    int FP0_BASIS =  Rmg_G->get_P0_BASIS(Rmg_G->get_default_FG_RATIO());
    double *vh_in = new double[FP0_BASIS];
    double *vxc_in = new double[FP0_BASIS * ct.nspin];


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

        Exx_scf = new Exxbase<OrbitalType>(*Kptr[0]->G, *Rmg_halfgrid, *Kptr[0]->L, "tempwave", Kptr[0]->nstates, occs.data(),
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

        // Adjust exx convergence threshold
        if(ct.xc_is_hybrid)
        {
            Exx_scf->vexx_RMS[ct.exx_steps] = 0.0;
            if(ct.exx_steps)
                ct.adaptive_thr_energy = std::min(1.0e-9, ct.vexx_rms / 100000.0);
            else
                ct.adaptive_thr_energy = 1.0e-9;
        }

        for (ct.scf_steps = 0, CONVERGED = false;
                ct.scf_steps < ct.max_scf_steps && !CONVERGED; ct.scf_steps++, ct.total_scf_steps++)
        {

            /* perform a single self-consistent step */
            step_time = my_crtc ();
            CONVERGED = Scf (vxc, vxc_in, vh, vh_in, ct.vh_ext, vnuc, rho, rho_oppo,
                    rhocore, rhoc, ct.spin_flag, ct.boundaryflag, Kptr, RMSdV);
            step_time = my_crtc () - step_time;

            // Save data to file for future restart at checkpoint interval if this is a quench run.
            // For Relaxation and molecular dynamics we save at the end of each ionic step.
            if (ct.checkpoint && (ct.scf_steps % ct.checkpoint == 0) && (ct.forceflag == MD_QUENCH))
                WriteRestart (ct.outfile, vh, rho, rho_oppo, vxc, Kptr);

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
            PlotConvergence(RMSdV, CONVERGED);
#endif

            // Write out progress info
            double elapsed_time = my_crtc() - start_time;
            ProgressTag(step_time, elapsed_time);

        }
        /* ---------- end scf loop ---------- */

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

            if(Exx_scf->vexx_RMS[ct.exx_steps] < ct.exx_convergence_criterion)
            { 
                printf(" Finished EXX outer loop in %3d exx steps, elapsed time = %6.2f, vexx_rms = %8.2e, total energy = %.*f Ha\n",
                        ct.exx_steps, exx_elapsed_time, Exx_scf->vexx_RMS[ct.exx_steps], 6, ct.TOTAL);
                ct.FOCK = f2;
                break;
            }
            else
            {
                printf(" Finished EXX inner loop in %3d scf steps, exx step time = %6.2f, vexx_rms = %8.2e, total energy = %.*f Ha\n",
                        ct.scf_steps, exx_step_time, Exx_scf->vexx_RMS[ct.exx_steps], 6, ct.TOTAL);
            }
        }

    }
    /* ---------- end exx loop ---------- */


    if(CONVERGED)
    {
        if(ct.rms < ct.thr_rms)
        {
            rmg_printf("\n Convergence criterion reached: potential RMS (%.2e) is lower than threshold (%.2e)\n", ct.rms, ct.thr_rms);
            if (pct.imgpe == 0 && pct.images == 1)
                fprintf(stdout,"\n Convergence criterion reached: potential RMS (%.2e) is lower than threshold (%.2e)", ct.rms, ct.thr_rms);
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
    rmg_printf ("final total energy = %16.8f Ha\n", ct.TOTAL);


    // Exact exchange integrals
    if(ct.exx_int_flag)
    {
        if(!ct.is_gamma) 
        {
            throw RmgFatalException() <<  "Exx integrals not programed for non-gamma kpoint" << __FILE__ << " at line " << __LINE__ << "\n";
        }
        std::vector<double> occs;
        occs.resize(Kptr[0]->nstates);
        for(int i=0;i < Kptr[0]->nstates;i++) occs[i] = Kptr[0]->Kstates[i].occupation[0];
        Exxbase<OrbitalType> Exx(*Kptr[0]->G, *Rmg_halfgrid, *Kptr[0]->L, "tempwave", Kptr[0]->nstates, occs.data(), 
                Kptr[0]->orbital_storage, ct.exx_mode);
        if(ct.exx_mode == EXX_LOCAL_FFT)
            Exx.WriteWfsToSingleFile();
        
        double *v_psi, *vxc_psi;
        int pbasis = Kptr[0]->pbasis;
        v_psi = new double[pbasis];
        vxc_psi = new double[pbasis]();
        int nstates = Kptr[0]->nstates;
        OrbitalType *Hcore = (OrbitalType *)GpuMallocManaged(nstates * nstates * sizeof(OrbitalType));

        bool is_xc_hybrid = ct.xc_is_hybrid;

        ct.xc_is_hybrid = false;

        GetVtotPsi (v_psi, vnuc, Rmg_G->default_FG_RATIO);

        ct.xc_is_hybrid = false;

        Kptr[0]->ComputeHcore(v_psi, vxc_psi, Hcore);
        Exx.SetHcore(Hcore, nstates);
        Exx.Vexx_integrals(ct.exx_int_file);
        ct.xc_is_hybrid = is_xc_hybrid;
        delete [] v_psi;
        delete [] vxc_psi;
    }


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
        double *vtot = new double[FP0_BASIS];
        for(int idx = 0; idx < FP0_BASIS; idx++) vtot[idx] = vh[idx] + vnuc[idx] + vxc[idx];
        Stress<OrbitalType> Stress_cal(Kptr, Rmg_L, *Rmg_G, *fine_pwaves, Atoms, Species, 
                ct.XC, vxc, rho, rhocore, vtot);
        delete [] vtot;
    }

    /*When running MD, force pointers need to be rotated before calculating new forces */
    if(compute_forces)
    {
        ct.fpt[0] = 0;
        ct.fpt[1] = 1;
        ct.fpt[2] = 2;
        ct.fpt[3] = 3;

        for (size_t ion = 0, i_end = Atoms.size(); ion < i_end; ++ion)
        {
            Atoms[ion].RotateForces();
        }
        Force (rho, rho_oppo, rhoc, vh, vh_in, vxc, vxc_in, vnuc, Kptr);
    }


    rmg_printf (" volume and energy per atom = %18.8f  %18.8f eV\n", Rmg_L.get_omega()*a0_A*a0_A*a0_A/Atoms.size(),ct.TOTAL * Ha_eV/Atoms.size());

    if (Verify("charge_analysis","Voronoi", Kptr[0]->ControlMap))
    {
        double timex = my_crtc ();
        double *localrho = new double[Atoms.size()];
        Voronoi_charge->LocalCharge(rho, localrho);
        for(size_t ion = 0; ion < Atoms.size(); ion++)
            Atoms[ion].partial_charge = Voronoi_charge->localrho_atomic[ion] - localrho[ion];
        if(ct.noncoll)
        {
            Voronoi_charge->LocalCharge(&rho[FP0_BASIS], localrho);
            for(size_t ion = 0; ion < Atoms.size(); ion++)
                printf("\n Voronoi rho_x %ld %f", ion, localrho[ion]);
            Voronoi_charge->LocalCharge(&rho[2*FP0_BASIS], localrho);
            for(size_t ion = 0; ion < Atoms.size(); ion++)
                printf("\n Voronoi rho_y %ld %f", ion, localrho[ion]);
            Voronoi_charge->LocalCharge(&rho[3*FP0_BASIS], localrho);
            for(size_t ion = 0; ion < Atoms.size(); ion++)
                printf("\n Voronoi rho_z %ld %f", ion, localrho[ion]);

        }
        delete [] localrho;

        rmg_printf("\n Vdd took %f seconds\n", my_crtc () - timex);
    }

    /*Calculate and write dipole moment if requested*/
    if (ct.dipole_moment)
    {
        double dipole[3];
        get_dipole(rho, dipole);

        double DEBYE_CONVERSION = 2.54174618772314;

        /*Now we need to convert to debye units */
        if (pct.imgpe==0)
        {
            printf("\n\n Dipole moment [Debye]: (%.3f,%.3f, %.3f)", 
                    DEBYE_CONVERSION *dipole[0], 
                    DEBYE_CONVERSION *dipole[1], 
                    DEBYE_CONVERSION *dipole[2]);
        }

    }

    WriteRestart (ct.outfile, vh, rho, rho_oppo, vxc, Kptr);
    if(ct.wannier90)
    {
        int scdm = ct.wannier90_scdm;
        double scdm_mu = ct.wannier90_scdm_mu;
        double scdm_sigma = ct.wannier90_scdm_sigma;
        int n_wannier = ct.num_wanniers;
 
        Wannier<OrbitalType> Wan(*Kptr[0]->G, *Kptr[0]->L, "tempwave", Kptr[0]->nstates, 
                n_wannier, scdm, scdm_mu, scdm_sigma, Kptr[0]->orbital_storage, Kptr);

    }

    /* output the forces */
    if (pct.imgpe == 0)
        write_force ();

    delete [] vxc_in;
    delete [] vh_in;
    return CONVERGED;


}                               /* end quench */


#if PLPLOT_LIBS
void PlotConvergence(std::vector<double> &RMSdV, bool CONVERGED)
{
    if(pct.imgpe == 0) {
        std::vector<double> x;
        std::string ConvergencePlot(ct.basename);
        ConvergencePlot = ConvergencePlot + ".rmsdv.png";
        std::string ConvergenceTitle("RMG convergence:");
        if(CONVERGED) {
            ConvergenceTitle = ConvergenceTitle + " quench completed.";
        }
        else {
            ConvergenceTitle = ConvergenceTitle + " quenching electrons.";
        }
        LinePlotLog10y(ConvergencePlot.c_str(), "SCF Steps", "log10(RMS[dV])", ConvergenceTitle.c_str(), x, RMSdV);
    }
}
#endif


void ChargeAnalysis(double *rho, std::unordered_map<std::string, InputKey *>& ControlMap, Voronoi &Vdd)
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
                Vdd.LocalCharge(rho, localrho);
                for(size_t ion = 0; ion < Atoms.size(); ion++)
                    Atoms[ion].partial_charge = Vdd.localrho_atomic[ion] - localrho[ion];
                delete [] localrho;
                WriteChargeAnalysis();
                rmg_printf("\n Vdd took %f seconds\n", my_crtc () - timex);
            }
        }
    }
}

