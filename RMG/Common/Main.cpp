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
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <fcntl.h> 
#include <unistd.h>
#include <unordered_map>
#include <csignal>

#include "const.h"
#include "RmgTimer.h"
#include "RmgException.h"
#include "rmgtypedefs.h"
#include "params.h"
#include "typedefs.h"
#include "rmg_error.h"
#include "transition.h"
#include "Kpoint.h"
#include "InputKey.h"
#include "blas.h"
#include "RmgThread.h"
#include "rmgthreads.h"


#include "../Headers/common_prototypes.h"
#include "../Headers/common_prototypes1.h"
#include "prototypes_tddft.h"
#include "Exxbase.h"


void initialize (int argc, char **argv);

template <typename OrbitalType> void run (Kpoint<OrbitalType> **Kptr);

void report (void);

void finish (void);

std::vector<ION> Atoms;
std::vector<SPECIES> Species;


/* Electronic charge density or charge density of own spin in polarized case */
double *rho;

/*  Electronic charge density of pposite spin density*/
double *rho_oppo;  


/* Core Charge density */
double *rhocore;


/* Compensating charge density */
double *rhoc;


/* Hartree potential */
double *vh;

/* Nuclear local potential */
double *vnuc;

/* Exchange-correlation potential */
double *vxc;

// Pointer to Kpoint class arrays for gamma and non-gamma
Kpoint<double> **Kptr_g;
Kpoint<std::complex<double> > **Kptr_c;

double *tau;
/* Main control structure which is declared extern in main.h so any module */
/* may access it.					                 */
CONTROL ct;

/* PE control structure which is also declared extern in main.h */
PE_CONTROL pct;

std::unordered_map<std::string, InputKey *> ControlMap;

std::atomic<bool> shutdown_request(false);

extern "C" void term_handler(int signal)
{
    shutdown_request.store(true); 
}

void CheckShutdown(void)
{
    if(shutdown_request.load())
    {
        DeleteNvmeArrays();
        for (int kpt = 0; kpt < ct.num_kpts_pe; kpt++)
        {

            int kpt1 = kpt + pct.kstart;
            if(ct.is_gamma)
            {
                Kptr_g[kpt1]->DeleteNvmeArrays();
            }
            else
            {
                Kptr_c[kpt1]->DeleteNvmeArrays();
            }
        }
        MPI_Abort( MPI_COMM_WORLD, 0 );
        kill(getpid(), SIGKILL);
    }
}

int main (int argc, char **argv)
{

    // Set branch type and save argc and argv in control structure
    ct.rmg_branch = RMG_BASE;
    ct.save_args(argc, argv);

    // Signal handlers to cleanup in case user terminates
    std::signal(SIGTERM, term_handler);
    std::signal(SIGINT, term_handler);

    RmgTimer *RT = new RmgTimer("1-TOTAL");
    char *tptr;

    /* Define a default output stream, gets redefined to log file later */

    ct.logfile = stdout;

// for RMG, the projectors |beta> are multiplied by exp(-ik.r) for non-gamma point
// for ON and NEGF, the phase exp(-ik.r) is in the matrix separtion.
    ct.proj_nophase = 0; 


    // Get RMG_MPI_THREAD_LEVEL environment variable
    ct.mpi_threadlevel = MPI_THREAD_SERIALIZED;
    if(NULL != (tptr = getenv("RMG_MPI_THREAD_LEVEL"))) {
        ct.mpi_threadlevel = atoi(tptr);
    }

    try {

        RmgTimer *RT1 =  new RmgTimer("1-TOTAL: Init");
        initialize (argc, argv);
        delete(RT1);



        RmgTimer *RT2 = new RmgTimer("1-TOTAL: run");
        if(ct.is_gamma)
            run<double> ((Kpoint<double> **)Kptr_g);
        else
            run<std::complex<double> >((Kpoint<std::complex<double>> **)Kptr_c);
        delete(RT2);

    }

    // Catch exceptions issued by us.
    catch(RmgFatalException const &e) {
        std::cout << e.rwhat() << std::endl;
        finish ();
        exit(0);
    }

    // By std
    catch (std::exception &e) {
        std::cout << "Caught a std exception: " << e.what () << std::endl;
        finish ();
        exit(0);
    } 

    // Catchall
    catch (...) {
        std::cout << "Caught an unknown exception of some type." << std::endl;
        finish ();
        exit(0);
    } 

    delete(RT);   // Destructor has to run before report
    report ();

    finish ();

    // Shutdown threads gracefully otherwise Cray perftools has issues
    RmgTerminateThreads();

}


void initialize(int argc, char **argv) 
{

    int FP0_BASIS;

    /* start the benchmark clock */
    ct.time0 = my_crtc ();
    RmgTimer *RT0 = new RmgTimer("2-Init");
    RmgTimer *RT = new RmgTimer("2-Init: KpointClass");

    /* Initialize all I/O including MPI group comms */
    /* Also reads control and pseudopotential files*/
    InitIo (argc, argv, ControlMap);

    FP0_BASIS = Rmg_G->get_P0_BASIS(Rmg_G->default_FG_RATIO);

    int num_images = pct.images;
    num_images = 1;
    lbfgs_init(Atoms.size(), num_images);

    int spinfac = 1;
    if(ct.spin_flag) spinfac = 2;
    rho = new double[spinfac * FP0_BASIS]();
    rhocore = new double[FP0_BASIS];
    rhoc = new double[FP0_BASIS];
    vh = new double[FP0_BASIS];
    vnuc = new double[FP0_BASIS];
    vxc = new double[spinfac * FP0_BASIS];
    if (ct.xctype == MGGA_TB09) 
    	tau = new double[FP0_BASIS];

    /* for spin polarized calculation set pointer to memory for density of the opposite spin */
    rho_oppo = rho + FP0_BASIS;

    /* Check if Bweight is needed */
    ct.need_Bweight = true;
//    if((ct.discretization == CENTRAL_DISCRETIZATION) && ct.norm_conserving_pp) ct.need_Bweight = false;
    if(ct.discretization == CENTRAL_DISCRETIZATION) ct.need_Bweight = false;

    /* Initialize some k-point stuff */
    Kptr_g = new Kpoint<double> * [ct.num_kpts_pe];
    Kptr_c = new Kpoint<std::complex<double> > * [ct.num_kpts_pe];

    ct.is_gamma = true;
    for (int kpt = 0; kpt < ct.num_kpts; kpt++) {
        double v1, v2, v3;
        v1 = twoPI * ct.kp[kpt].kpt[0] / Rmg_L.get_xside();
        v2 = twoPI * ct.kp[kpt].kpt[1] / Rmg_L.get_yside();
        v3 = twoPI * ct.kp[kpt].kpt[2] / Rmg_L.get_zside();

        ct.kp[kpt].kvec[0] = v1;
        ct.kp[kpt].kvec[1] = v2;
        ct.kp[kpt].kvec[2] = v3;
        ct.kp[kpt].kmag = v1 * v1 + v2 * v2 + v3 * v3;

        if(ct.kp[kpt].kmag != 0.0) ct.is_gamma = false;
    }

    if(ct.is_gamma) 
    {
        ct.is_use_symmetry = 0;
    }

    for (int kpt = 0; kpt < ct.num_kpts_pe; kpt++)
    {

        int kpt1 = kpt + pct.kstart;
        if(ct.is_gamma) {

            // Gamma point
            Kptr_g[kpt] = new Kpoint<double> (ct.kp[kpt1], kpt, pct.grid_comm, Rmg_G, Rmg_T, &Rmg_L, ControlMap);

        }
        else {

            // General case
            Kptr_c[kpt] = new Kpoint<std::complex<double>> (ct.kp[kpt1], kpt, pct.grid_comm, Rmg_G, Rmg_T, &Rmg_L, ControlMap);

        }
        ct.kp[kpt].kidx = kpt;
    }


    MPI_Barrier (pct.img_comm);

    /* Record the time it took from the start of run until we hit init */
    delete(RT);

    /* Perform any necessary initializations */
    if(ct.is_gamma) {
        Init (vh, rho, rho_oppo, rhocore, rhoc, vnuc, vxc, Kptr_g);
    }
    else {
        Init (vh, rho, rho_oppo, rhocore, rhoc, vnuc, vxc, Kptr_c);
    }


    /* Flush the results immediately */
    fflush (NULL);


    /* Wait until everybody gets here */
    /* MPI_Barrier(MPI_COMM_WORLD); */
    MPI_Barrier(pct.img_comm);

    delete(RT0);

}

template <typename OrbitalType> void run (Kpoint<OrbitalType> **Kptr)
{


    /* Dispatch to the correct driver routine */
    switch (ct.forceflag)
    {

        case MD_QUENCH:            /* Quench the electrons */
            if (ct.xctype == MGGA_TB09)
                //relax_tau (0, states, vxc, vh, vnuc, rho, rho_oppo, rhocore, rhoc, tau);
                ;
            else 
            Relax (0, vxc, vh, vnuc, rho, rho_oppo, rhocore, rhoc, Kptr);
            break;

        case MD_FASTRLX:           /* Fast relax */
            Relax (ct.max_md_steps, vxc, vh, vnuc, rho, rho_oppo, rhocore, rhoc, Kptr);
            break;

        case NEB_RELAX:           /* nudged elastic band relax */
            NEB_relax (ct.max_neb_steps, vxc, vh, vnuc, rho, rho_oppo, rhocore, rhoc, Kptr);
            break;

        case MD_CVE:               /* molecular dynamics */
        case MD_CVT:
        case MD_CPT:
            ct.fpt[0] = 0;  // Eventually fix all references to fpt in the code and this will not be needed
            ct.fpt[1] = 1;
            ct.fpt[2] = 2;
            ct.fpt[3] = 3;
            MolecularDynamics (Kptr, vxc, vh, vnuc, rho, rho_oppo, rhoc, rhocore);
            break;

        case BAND_STRUCTURE:
            BandStructure (Kptr, vxc, vh, vnuc);
            if(ct.rmg2bgw) WriteBGW_Rhog(rho, rho_oppo);
            OutputBandPlot(Kptr);
            return;

        case TDDFT:
            if(!ct.restart_tddft) Relax (0, vxc, vh, vnuc, rho, rho_oppo, rhocore, rhoc, Kptr);
            RmgTddft (vxc, vh, vnuc, rho, rho_oppo, rhocore, rhoc, Kptr);
            break;
        
        case Exx_only:
            {
                std::vector<double> occs;
                occs.resize(Kptr[0]->nstates);
                for(int i=0;i < Kptr[0]->nstates;i++) occs[i] = Kptr[0]->Kstates[i].occupation[0];
                Exxbase<OrbitalType> Exx(*Kptr[0]->G, *Kptr[0]->L, "tempwave", Kptr[0]->nstates, occs.data(), 
                        Kptr[0]->orbital_storage, ct.exx_mode);
                Exx.Vexx_integrals(ct.exx_int_file);
                break;
            }

        default:
            rmg_error_handler (__FILE__, __LINE__, "Undefined MD method");


    }


    if(Verify ("output_rho_xsf", true, Kptr[0]->ControlMap))
        Output_rho_xsf(rho, pct.grid_comm);

}                               /* end run */

void report ()
{

    /* write planar averages of quantities */
    if (ct.zaverage == 1)
    {
        /* output the average potential */
        write_avgv (vh, vnuc);
        write_avgd (rho);
    }
    else if (ct.zaverage == 2)
    {
        //write_zstates (states);
        ;
    }


    /* If milliken population info is requested then compute and output it */
    /*if (ct.domilliken)
      mulliken (states);*/


    if (ct.xctype == MGGA_TB09) 
        delete [] tau;

    /* Write timing information */
    if(pct.imgpe == 0) fclose(ct.logfile);
    int override_rank = 0;
    if(pct.imgpe==0) MPI_Comm_rank (pct.img_comm, &override_rank);
    int num_owned_ions;
    if(ct.is_gamma)
    {
        num_owned_ions = Kptr_g[0]->BetaProjector->get_num_owned_ions();
    }
    else
    {
        num_owned_ions = Kptr_c[0]->BetaProjector->get_num_owned_ions();
    }
    RmgPrintTimings(pct.img_comm, ct.logname, ct.scf_steps, num_owned_ions * ct.num_kpts_pe, override_rank);


}                               /* end report */


void finish ()
{

    DeleteNvmeArrays();
    for (int kpt = 0; kpt < ct.num_kpts_pe; kpt++)
    {

        int kpt1 = kpt + pct.kstart;
        if(ct.is_gamma)
        {
            Kptr_g[kpt1]->DeleteNvmeArrays();
        }
        else
        {
            Kptr_c[kpt1]->DeleteNvmeArrays();
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    /*Exit MPI */
    MPI_Finalize ();

#if GPU_ENABLED
    //cublasDestroy(ct.cublas_handle);
    //cudaDeviceReset();
#endif

}                               /* end finish */


/******/


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

// Instantiate gamma and non-gamma versions
template bool Quench<double> (double *, double *, double *, double *, double *, double *, double *, Kpoint<double> **Kptr, bool);
template bool Quench<std::complex<double> > (double *, double *, double *, double *, double *, double *, double *, Kpoint<std::complex<double>> **Kptr, bool);

template <typename OrbitalType> bool Quench (double * vxc, double * vh, double * vnuc, double * rho,
        double * rho_oppo, double * rhocore, double * rhoc, Kpoint<OrbitalType> **Kptr, bool compute_forces)
{

    bool CONVERGED;
    static std::vector<double> RMSdV;

    int FP0_BASIS =  Rmg_G->get_P0_BASIS(Rmg_G->get_default_FG_RATIO());
    double *vh_in = new double[FP0_BASIS];
    double *vxc_in = new double[FP0_BASIS];

    /* ---------- begin scf loop ---------- */

    double start_time = my_crtc ();
    double step_time;
    double elapsed_time;

    for (ct.scf_steps = 0, CONVERGED = false;
            ct.scf_steps < ct.max_scf_steps && !CONVERGED; ct.scf_steps++, ct.total_scf_steps++)
    {



        /* perform a single self-consistent step */
        step_time = my_crtc ();
        CONVERGED = Scf (vxc, vxc_in, vh, vh_in, ct.vh_ext, vnuc, rho, rho_oppo, rhocore, rhoc, ct.spin_flag, ct.boundaryflag, Kptr, RMSdV);
        step_time = my_crtc () - step_time;

        // Save data to file for future restart at checkpoint interval if this is a quench run.
        // For Relaxation and molecular dynamics we save at the end of each ionic step.
        if (ct.checkpoint)
            if ((ct.scf_steps % ct.checkpoint == 0) && (ct.forceflag == MD_QUENCH))
                WriteRestart (ct.outfile, vh, rho, rho_oppo, vxc, Kptr);

        /* output the eigenvalues with occupations */
        if (ct.write_eigvals_period)
        {
            if (ct.scf_steps % ct.write_eigvals_period == 0)
            {
                if (pct.imgpe == 0)
                {
                    OutputEigenvalues (Kptr, 0, ct.scf_steps);
                    rmg_printf ("\nTotal charge in supercell = %16.8f\n", ct.tcharge);
                }
            }
        }


        /*Perform charge analysis if requested*/
        if (ct.charge_analysis_period)
        {
            if (ct.scf_steps % ct.charge_analysis_period == 0)
            {
                if (Verify("charge_analysis","Voronoi", Kptr[0]->ControlMap))
                {
                    double timex = my_crtc ();
                    Vdd(rho);
                    WriteChargeAnalysis();
                    rmg_printf("\n Vdd took %f seconds\n", my_crtc () - timex);
                }
            }
        }

#if PLPLOT_LIBS
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
#endif

        elapsed_time = my_crtc() - start_time;
        if (pct.imgpe == 0) {
            rmg_printf (" quench: [md: %3d/%-d  scf: %3d/%-d  step time: %6.2f  scf time: %8.2f secs  RMS[dV]: %8.2e ]\n\n\n",
                    ct.md_steps, ct.max_md_steps, ct.scf_steps, ct.max_scf_steps, step_time, elapsed_time, ct.rms);

            /*Also print to stdout*/
            if(pct.images == 1)
                fprintf (stdout,"\n quench: [md: %3d/%-d  scf: %3d/%-d  step time: %6.2f  scf time: %8.2f secs  RMS[dV]: %8.2e ]",
                        ct.md_steps, ct.max_md_steps, ct.scf_steps, ct.max_scf_steps, step_time, elapsed_time, ct.rms);
        }

    }
    /* ---------- end scf loop ---------- */



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
        //progress_tag ();
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
    // Experimental for now. Exchange correlation type must be manually set to
    // gaupbe in the input file (and gaupbe is the only divergence type supported).
    if(ct.exx_int_flag)
    {
    }


    /* output final eigenvalues with occupations */
    OutputEigenvalues (Kptr, 0, ct.scf_steps);
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
        Vdd(rho);
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

    }                               /* end getpoi_bc.c */





    /* output the forces */
    if (pct.imgpe == 0)
        write_force ();

    delete [] vxc_in;
    delete [] vh_in;
    return CONVERGED;


}                               /* end quench */

