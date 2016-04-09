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
#include "../Headers/prototypes.h"
#include "RmgParallelFft.h"




static int firststep = true;

template bool Scf<double> (double *, double *, double *,
          double *, double *, double *, double *, double *, int ,
          int , int , int , Kpoint<double> **, std::vector<double> &);
template bool Scf<std::complex<double> > (double *, double *, double *,
          double *, double *, double *, double *, double *, int ,
          int , int , int , Kpoint<std::complex<double>> **, std::vector<double> &);

template <typename OrbitalType> bool Scf (double * vxc, double * vh, double *vh_ext,
          double * vnuc, double * rho, double * rho_oppo, double * rhocore, double * rhoc, int spin_flag,
          int hartree_min_sweeps, int hartree_max_sweeps , int boundaryflag, Kpoint<OrbitalType> **Kptr, std::vector<double>& RMSdV)
{

    RmgTimer RT0("Scf steps"), *RT1;
    int nspin = (spin_flag + 1);
    bool CONVERGED = false;
    double t3;
    double *vtot, *vtot_psi, *new_rho;
    double time1;
    double t[3];                  /* SCF checks and average potential */
    double max_unocc_res = 0.0;


    int P0_BASIS, FP0_BASIS;

    /* to hold the send data and receive data of eigenvalues */
    double *rho_tot=NULL;   
    
    time1 = my_crtc ();

    P0_BASIS =  Rmg_G->get_P0_BASIS(1);
    FP0_BASIS = Rmg_G->get_P0_BASIS(Rmg_G->default_FG_RATIO);

    /* allocate memory for eigenvalue send array and receive array */
    if (spin_flag)
    {
    	rho_tot = new double[FP0_BASIS];
    }

    new_rho = new double[FP0_BASIS];
    vtot = new double[FP0_BASIS];
    vtot_psi = new double[P0_BASIS];

    /* save old vhxc + vnuc */
    for (int idx = 0; idx < FP0_BASIS; idx++) {
        vtot[idx] = vxc[idx] + vh[idx] + vnuc[idx];
    }

    /* Calculate total energy and the new exchange-correlation potential */
    GetTe (rho, rho_oppo, rhocore, rhoc, vh, vxc, Kptr, !ct.scf_steps);


    double rms_target = std::max(ct.rms/ct.hartree_rms_ratio, 1.0e-12);
    if (spin_flag)        
    {
	/*calculate the total charge density in order to calculate hartree potential*/
	for (int idx = 0; idx < FP0_BASIS; idx++) {
            rho_tot[idx] = rho[idx] + rho_oppo[idx];
        }
	
	/* Generate hartree potential */
        RT1 = new RmgTimer("Scf steps: Hartree");
        int dimx = Rmg_G->get_PX0_GRID(Rmg_G->get_default_FG_RATIO());
        int dimy = Rmg_G->get_PY0_GRID(Rmg_G->get_default_FG_RATIO());
        int dimz = Rmg_G->get_PZ0_GRID(Rmg_G->get_default_FG_RATIO());
        int pbasis = dimx * dimy * dimz;
        double *rho_neutral = new double[pbasis];

        /* Subtract off compensating charges from rho */
        for (int idx = 0; idx < pbasis; idx++) {
            rho_neutral[idx] = rho[idx] - rhoc[idx];
        }

        double residual = CPP_get_vh (Rmg_G, &Rmg_L, Rmg_T, rho_neutral, vh_ext, hartree_min_sweeps, 
                    hartree_max_sweeps, ct.poi_parm.levels, ct.poi_parm.gl_pre, 
                    ct.poi_parm.gl_pst, ct.poi_parm.mucycles, rms_target,
                    ct.poi_parm.gl_step, ct.poi_parm.sb_step, boundaryflag, Rmg_G->get_default_FG_RATIO(), false);
        rmg_printf("Hartree residual = %14.6e\n", residual);
 
        /* Pack the portion of the hartree potential used by the wavefunctions
         * back into the wavefunction hartree array. */
        CPP_pack_dtos (Rmg_G, vh, vh_ext, dimx, dimy, dimz, boundaryflag);
        delete [] rho_neutral;

        delete(RT1);

     }  	
    else
    {
    	/* Generate hartree potential */
        RT1 = new RmgTimer("Scf steps: Hartree");
    	get_vh (rho, rhoc, vh, hartree_min_sweeps, hartree_max_sweeps, ct.poi_parm.levels, rms_target, boundaryflag);
        delete(RT1);
    }


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

    GlobalSums (t, 3, pct.img_comm);
    t[0] *= get_vel_f();
    
    /* get the averaged value over each spin and each fine grid */
    t[1] = sqrt (t[1] / ((double) (nspin * ct.psi_fnbasis)));  
    t[2] /= ((double) (nspin * ct.psi_fnbasis));   
    
    ct.rms = t[1];

    if (!firststep)
    {
        rmg_printf ("\n");
        //progress_tag ();
        rmg_printf ("SCF CHECKS: <rho dv>  = %15.8e\n", t[0]);
        //progress_tag ();
        rmg_printf ("SCF CHECKS: RMS[dv]   = %15.8e\n", t[1]);
        RMSdV.emplace_back(t[1]);
        //progress_tag ();
        rmg_printf ("AVERAGE POTENTIAL <V> = %15.8e\n", t[2]);
    }

    if(!Verify ("freeze_occupied", true, Kptr[0]->ControlMap)) {

        if (!firststep && t[1] < ct.thr_rms) CONVERGED = true;

    }

    // Transfer vtot from the fine grid to the wavefunction grid
    GetVtotPsi (vtot_psi, vtot, Rmg_G->default_FG_RATIO);

    /*Generate the Dnm_I */
    get_ddd (vtot);

    // Loop over k-points
    for(int kpt = 0;kpt < ct.num_kpts;kpt++) {

        if (Verify ("kohn_sham_solver","multigrid", Kptr[0]->ControlMap) || (ct.scf_steps < 4)) {
            MgridSubspace(Kptr[kpt], vtot_psi);
        }
        else if(Verify ("kohn_sham_solver","davidson", Kptr[0]->ControlMap)) {
            int notconv;
            Davidson(Kptr[kpt], vtot_psi, notconv);
        }

    } // end loop over kpoints


    if (spin_flag)
        GetOppositeEigvals (Kptr);


    /* Take care of occupation filling */
    ct.efermi = Fill (Kptr, ct.occ_width, ct.nel, ct.occ_mix, ct.num_states, ct.occ_flag);


    if (ct.occ_flag == 1 )
    {
        rmg_printf ("\n");
        //progress_tag ();
        rmg_printf ("FERMI ENERGY = %15.8f eV\n", ct.efermi * Ha_eV);
    }


    if (firststep)
        firststep = false;

    /* Generate new density */
    RT1 = new RmgTimer("Scf steps: Get rho");
    GetNewRho(Kptr, new_rho);

    /*Takes care of mixing and checks whether the charge density is negative*/
    MixRho(new_rho, rho, rhocore, Kptr[0]->ControlMap);
    

    if (spin_flag)
	get_rho_oppo (rho,  rho_oppo);
    
    delete(RT1);


    /* Make sure we see output, e.g. so we can shut down errant runs */
    fflush( ct.logfile );
#if !(__CYGWIN__ || __MINGW64__ || _WIN32 || _WIN64)
    fsync( fileno(ct.logfile) );
#endif

    /* release memory */
    delete [] new_rho;
    delete [] vtot;
    delete [] vtot_psi;

    if (spin_flag)
    {
    	delete [] rho_tot;
    }

    rmg_printf("\n SCF STEP TIME = %10.2f\n",my_crtc () - time1);

    if(Verify ("freeze_occupied", true, Kptr[0]->ControlMap)) {

        if(!firststep && (max_unocc_res < ct.gw_threshold)) {
            rmg_printf("\nGW: convergence criteria of %10.5e has been met.\n", ct.gw_threshold);
            rmg_printf("GW:  Highest occupied orbital index              = %d\n", Kptr[0]->highest_occupied);
//            rmg_printf("GW:  Highest unoccupied orbital meeting criteria = %d\n", Kptr[0]->max_unocc_res_index);

            CONVERGED = true;
        }

    }

    return CONVERGED;
}                               /* end scf */


/******/
