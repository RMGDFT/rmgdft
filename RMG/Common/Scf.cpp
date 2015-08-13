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
#include "BaseThread.h"
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
#include "../Headers/prototypes.h"




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

    RmgTimer RT0("Scf steps");
    int st1, diag_this_step;
    int nspin = (spin_flag + 1);
    bool CONVERGED = false;
    double t3;
    double *vtot, *vtot_psi, *new_rho;
    double time1;
    double t[3];                  /* SCF checks and average potential */
    double mean_occ_res = DBL_MAX;
    double mean_unocc_res = DBL_MAX;
    double max_occ_res = 0.0;
    double max_unocc_res = 0.0;
    double min_occ_res = DBL_MAX;
    double min_unocc_res = DBL_MAX;

    int ist, istop, P0_BASIS, FP0_BASIS;
    BaseThread *T = BaseThread::getBaseThread(0);

    /* to hold the send data and receive data of eigenvalues */
    double *rho_tot;   
    
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

    /* Generate exchange-correlation potential */
    RmgTimer *RT1 = new RmgTimer("Scf steps: Get vxc");
    get_vxc (rho, rho_oppo, rhocore, vxc);
    delete(RT1);

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
                    ct.poi_parm.gl_pst, ct.poi_parm.mucycles, ct.rms/ct.hartree_rms_ratio,
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
        int dimx = Rmg_G->get_PX0_GRID(Rmg_G->get_default_FG_RATIO());
        int dimy = Rmg_G->get_PY0_GRID(Rmg_G->get_default_FG_RATIO());
        int dimz = Rmg_G->get_PZ0_GRID(Rmg_G->get_default_FG_RATIO());
        int pbasis = dimx * dimy * dimz;
        double *rho_tf = new double[pbasis];
        double *rho_tmp = new double[pbasis];

#if 1
	get_dipole (rho);
        rmg_printf("Calling get_tf_rho\n");
	get_tf_rho(rho_tf);
        rmg_printf("get_tf_rho Done\n");
    
	/* check TF charge */
	t[0] = 0.0;

	for (int idx = 0; idx < FP0_BASIS; idx++)
	{
	    rho_tmp[idx] = rho[idx] + rho_tf[idx];
	    t[0] += rho_tf[idx];
	}                           /* idx */
	
	GlobalSums (t, 1, pct.img_comm);
	t[0] *= get_vel_f();
        rmg_printf("tf_rho sum is %.8e\n", t[0]);
#endif
    	get_vh (rho_tmp, rhoc, vh, hartree_min_sweeps, 100, ct.poi_parm.levels, 10e-9, boundaryflag);
        delete(RT1);
        delete [] rho_tf;
        delete [] rho_tmp;
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

    get_vtot_psi (vtot_psi, vtot, get_FG_RATIO());

    /*Generate the Dnm_I */
    get_ddd (vtot);

    // Loop over k-points
    for(int kpt = 0;kpt < ct.num_kpts;kpt++) {

        mean_occ_res = 0.0;
        mean_unocc_res = 0.0;

        for(int vcycle = 0;vcycle < ct.eig_parm.mucycles;vcycle++) {

            // Update betaxpsi        
            RT1 = new RmgTimer("Scf steps: Beta x psi");
            Betaxpsi (Kptr[kpt]);
            delete(RT1);

            AppNls(Kptr[kpt], Kptr[kpt]->oldsint_local);
//            Kptr[kpt]->mix_betaxpsi(0);

            /* Update the wavefunctions */
            RT1 = new RmgTimer("Scf steps: Mg_eig");
            istop = Kptr[kpt]->nstates / T->get_threads_per_node();
            istop = istop * T->get_threads_per_node();

            for(st1=0;st1 < istop;st1+=T->get_threads_per_node()) {
              SCF_THREAD_CONTROL thread_control[MAX_RMG_THREADS];
              for(ist = 0;ist < T->get_threads_per_node();ist++) {
                  thread_control[ist].job = HYBRID_EIG;
                  thread_control[ist].vtot = vtot_psi;
                  thread_control[ist].sp = &Kptr[kpt]->Kstates[st1 + ist];
                  thread_control[ist].p3 = (void *)Kptr[kpt];
                  T->set_pptr(ist, &thread_control[ist]);
              }

              // Thread tasks are set up so run them
              T->run_thread_tasks(T->get_threads_per_node());

            }

            // Process any remaining states in serial fashion
            for(st1 = istop;st1 < Kptr[kpt]->nstates;st1++) {
                if(ct.is_gamma) {
                    MgEigState<double,float> ((Kpoint<double> *)Kptr[kpt], (State<double> *)&Kptr[kpt]->Kstates[st1], vtot_psi);
                }
                else {
                    MgEigState<std::complex<double>, std::complex<float> > ((Kpoint<std::complex<double>> *)Kptr[kpt], (State<std::complex<double> > *)&Kptr[kpt]->Kstates[st1], vtot_psi);
                }
            }
            delete(RT1);

        }

        if(Verify ("freeze_occupied", true, Kptr[kpt]->ControlMap)) {

            // Orbital residual measures (used for some types of calculations
            Kptr[kpt]->max_unocc_res_index = (int)(ct.gw_residual_fraction * (double)Kptr[kpt]->nstates);
            Kptr[kpt]->mean_occ_res = 0.0;
            Kptr[kpt]->min_occ_res = DBL_MAX;
            Kptr[kpt]->max_occ_res = 0.0;
            Kptr[kpt]->mean_unocc_res = 0.0;
            Kptr[kpt]->min_unocc_res = DBL_MAX;
            Kptr[kpt]->max_unocc_res = 0.0;
            Kptr[kpt]->highest_occupied = 0;
            for(int istate = 0;istate < Kptr[kpt]->nstates;istate++) {
                if(Kptr[kpt]->Kstates[istate].occupation[0] > 0.0) {
                    Kptr[kpt]->mean_occ_res += Kptr[kpt]->Kstates[istate].res;
                    mean_occ_res += Kptr[kpt]->Kstates[istate].res;
                    if(Kptr[kpt]->Kstates[istate].res >  Kptr[kpt]->max_occ_res)  Kptr[kpt]->max_occ_res = Kptr[kpt]->Kstates[istate].res;
                    if(Kptr[kpt]->Kstates[istate].res <  Kptr[kpt]->min_occ_res)  Kptr[kpt]->min_occ_res = Kptr[kpt]->Kstates[istate].res;
                    if(Kptr[kpt]->Kstates[istate].res >  max_occ_res)  max_occ_res = Kptr[kpt]->Kstates[istate].res;
                    if(Kptr[kpt]->Kstates[istate].res <  min_occ_res)  min_occ_res = Kptr[kpt]->Kstates[istate].res;
                    Kptr[kpt]->highest_occupied = istate;
                }
                else {
                    if(istate <= Kptr[kpt]->max_unocc_res_index) {
                        Kptr[kpt]->mean_unocc_res += Kptr[kpt]->Kstates[istate].res;
                        mean_unocc_res += Kptr[kpt]->Kstates[istate].res;
                        if(Kptr[kpt]->Kstates[istate].res >  Kptr[kpt]->max_unocc_res)  Kptr[kpt]->max_unocc_res = Kptr[kpt]->Kstates[istate].res;
                        if(Kptr[kpt]->Kstates[istate].res <  Kptr[kpt]->min_unocc_res)  Kptr[kpt]->min_unocc_res = Kptr[kpt]->Kstates[istate].res;
                        if(Kptr[kpt]->Kstates[istate].res >  max_unocc_res)  max_unocc_res = Kptr[kpt]->Kstates[istate].res;
                        if(Kptr[kpt]->Kstates[istate].res <  min_unocc_res)  min_unocc_res = Kptr[kpt]->Kstates[istate].res;
                    }
                }
            }
            Kptr[kpt]->mean_occ_res = Kptr[kpt]->mean_occ_res / (double)(Kptr[kpt]->highest_occupied + 1);
            Kptr[kpt]->mean_unocc_res = Kptr[kpt]->mean_unocc_res / (double)(Kptr[kpt]->max_unocc_res_index -(Kptr[kpt]->highest_occupied + 1));
            mean_occ_res = mean_occ_res / (double)(ct.num_kpts*(Kptr[kpt]->highest_occupied + 1));
            mean_unocc_res = mean_unocc_res / (double)(ct.num_kpts*Kptr[kpt]->max_unocc_res_index -(Kptr[kpt]->highest_occupied + 1));

            rmg_printf("Mean/Min/Max unoccupied wavefunction residual for kpoint %d  =  %10.5e  %10.5e  %10.5e\n", kpt, Kptr[kpt]->mean_unocc_res, Kptr[kpt]->min_unocc_res, Kptr[kpt]->max_unocc_res);

        }


        /*wavefunctions have changed, projectors have to be recalculated */
        RT1 = new RmgTimer("Scf steps: Beta x psi");
        Betaxpsi (Kptr[kpt]);
        delete(RT1);


        /* Now we orthognalize and optionally do subspace diagonalization
         * In the gamma point case, orthogonalization is not required when doing subspace diagonalization
         * For non-gamma point we have to do first orthogonalization and then, optionally subspace diagonalization
         * the reason is for non-gamma subdiag is not coded to solve generalized eigenvalue problem, it can
         * only solve the regular eigenvalue problem and that requires that wavefunctions are orthogonal to start with.*/

        diag_this_step = (ct.diag && ct.scf_steps % ct.diag == 0 && ct.scf_steps < ct.end_diag);

        /* do diagonalizations if requested, if not orthogonalize */
        if (diag_this_step) {

            Subdiag (Kptr[kpt], vh, vnuc, vxc, ct.subdiag_driver);
            Betaxpsi (Kptr[kpt]);
            Kptr[kpt]->mix_betaxpsi(0);
            // Projectors are rotated along with orbitals in Subdiag so no need to recalculate
            // after diagonalizing.

        }
        else {

            RT1 = new RmgTimer("Scf steps: Orthogonalization");
            Kptr[kpt]->orthogonalize(Kptr[kpt]->orbital_storage);
            delete(RT1);

            // wavefunctions have changed, projectors have to be recalculated */
            RT1 = new RmgTimer("Scf steps: Beta x psi");
            Betaxpsi (Kptr[kpt]);
            delete(RT1);
            Kptr[kpt]->mix_betaxpsi(1);

        }
            

        /* If sorting is requested then sort the states. */
        if (ct.sortflag) {
            Kptr[kpt]->sort_orbitals();
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
    MixRho(new_rho, rho, rhocore, FP0_BASIS, get_FPX0_GRID(), get_FPY0_GRID(), get_FPZ0_GRID(), Kptr[0]->ControlMap);
    

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
            rmg_printf("GW:  Highest unoccupied orbital meeting criteria = %d\n", Kptr[0]->max_unocc_res_index);

            CONVERGED = true;
        }

    }

    return CONVERGED;
}                               /* end scf */


/******/
