/************************** SVN Revision Information **************************
 **    $Id: scf.c 2253 2014-04-02 16:16:15Z ebriggs $    **
******************************************************************************/

/****f* QMD-MGDFT/scf.c *****
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
 *   void scf(STATE *states, double *vxc, double *vh, double *vnuc,
 *            double *rho, double *rhocore, double *rhoc, int *CONVERGENCE)
 *   Performs a single self consistent step over a full set of orbitals.
 *   This includes a loop over k-points.
 * INPUTS
 *   states: points to orbital structure
 *   vxc: exchange correlation potential
 *   vh:  Hartree potential
 *   vnuc: pseudopotential
 *   rho: total valence charge density
 *   rhocore:  core charge density
 *   rhoc: Gaussian compensating charge density
 * OUTPUT
 *   states, vxc, vh, rho are updated
 *   CONVERGENCE: 1 converged, 0 not
 * PARENTS
 *   cdfastrlx.c fastrlx.c moldyn.c quench.c
 * CHILDREN
 *   get_vxc.c get_vh.c mg_eig_state.c ortho_full.c fill.c get_rho.c
 * SOURCE
 */


#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "const.h"
#include "common_prototypes.h"
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

// Transitional headers
#include "transition.h"
#include "typedefs.h"


static int firststep = true;

template <typename ScfType> bool scf (Kpoint<ScfType> * kpoint, double * vxc, double * vh, double *vh_ext,
          double * vnuc, double * rho, double * rho_oppo, double * rhocore, double * rhoc, int spin_flag,
          int hartree_min_sweeps, int hartree_max_sweeps , int boundaryflag)
{

    RmgTimer RT0("Scf steps");
    int st1, idx, diag_this_step;
    int nspin = (spin_flag + 1);
    bool CONVERGED = false;
    double t3;
    double *vtot, *vtot_psi, *new_rho;
    double time1, time2, time3;
    double t[3];                  /* SCF checks and average potential */
    int ist, istop, P0_BASIS, FP0_BASIS;
    BaseThread *T = BaseThread::getBaseThread(0);


    /* to hold the send data and receive data of eigenvalues */
    double *rho_tot;   
    
    time3 = my_crtc ();

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
    for (idx = 0; idx < FP0_BASIS; idx++)
        vtot[idx] = vxc[idx] + vh[idx] + vnuc[idx];


    time1 = my_crtc (); 
 
    /* Generate exchange-correlation potential */
    RmgTimer *RT1 = new RmgTimer("Scf steps: Get vxc");
    get_vxc (rho, rho_oppo, rhocore, vxc);
    delete(RT1);

    if (spin_flag)        
    {
	/*calculate the total charge density in order to calculate hartree potential*/
	for (idx = 0; idx < FP0_BASIS; idx++)
		rho_tot[idx] = rho[idx] + rho_oppo[idx];
	
	/* Generate hartree potential */
        RT1 = new RmgTimer("Scf steps: Hartree");
//        get_vh (rho_tot, rhoc, vh, ct.hartree_min_sweeps, ct.hartree_max_sweeps, ct.poi_parm.levels, ct.rms/ct.hartree_rms_ratio, ct.boundaryflag);
        int dimx = Rmg_G->get_PX0_GRID(Rmg_G->get_default_FG_RATIO());
        int dimy = Rmg_G->get_PY0_GRID(Rmg_G->get_default_FG_RATIO());
        int dimz = Rmg_G->get_PZ0_GRID(Rmg_G->get_default_FG_RATIO());
        int idx, pbasis = dimx * dimy * dimz;
        double *rho_neutral = new double[pbasis];

        /* Subtract off compensating charges from rho */
        for (idx = 0; idx < pbasis; idx++)
            rho_neutral[idx] = rho[idx] - rhoc[idx];

        double residual = CPP_get_vh (Rmg_G, &Rmg_L, Rmg_T, rho_neutral, vh_ext, hartree_min_sweeps, hartree_max_sweeps, ct.poi_parm.levels, ct.poi_parm.gl_pre,
                    ct.poi_parm.gl_pst, ct.poi_parm.mucycles, ct.rms/ct.hartree_rms_ratio,
                    ct.poi_parm.gl_step, ct.poi_parm.sb_step, boundaryflag, Rmg_G->get_default_FG_RATIO(), false);
        //cout << "Hartree residual = " << residual << endl;
 
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
    	get_vh (rho, rhoc, vh, hartree_min_sweeps, hartree_max_sweeps, ct.poi_parm.levels, ct.rms/ct.hartree_rms_ratio, boundaryflag);
        delete(RT1);
    }



    /* check convergence */
    t[0] = t[1] = t[2] = 0.0;


    for (idx = 0; idx < FP0_BASIS; idx++)
    {
        t3 = -vtot[idx];
        vtot[idx] = vxc[idx] + vh[idx] + vnuc[idx];
        t3 += vtot[idx];
        t[0] += rho[idx] * t3;
        t[1] += t3 * t3;
        t[2] += vh[idx];
    }                           /* idx */


    
    
    idx = 3;
    GlobalSums (t, idx, pct.img_comm);
    t[0] *= get_vel_f();
    
    /* get the averaged value over each spin and each fine grid */
    t[1] = sqrt (t[1] / ((double) (nspin * ct.psi_fnbasis)));  
    t[2] /= ((double) (nspin * ct.psi_fnbasis));   
    
    ct.rms = t[1];

    if (!firststep)
    {
        printf ("\n");
        //progress_tag ();
        printf ("SCF CHECKS: <rho dv>  = %15.8e\n", t[0]);
        //progress_tag ();
        printf ("SCF CHECKS: RMS[dv]   = %15.8e\n", t[1]);
        //progress_tag ();
        printf ("AVERAGE POTENTIAL <V> = %15.8e\n", t[2]);
    }

    if (!firststep && t[1] < ct.thr_rms)
    {
	    CONVERGED = true;
    }

    get_vtot_psi (vtot_psi, vtot, get_FG_RATIO());

    /*Generate the Dnm_I */
    get_ddd (vtot);

    for(int vcycle = 0;vcycle < ct.eig_parm.mucycles;vcycle++) {

        RT1 = new RmgTimer("Scf steps: Beta x psi");
        betaxpsi (states);
        delete(RT1);

#if BATCH_NLS
        app_nls_batch (states, pct.nv, pct.ns, pct.Bns, pct.oldsintR_local);
#endif

        /* Update the wavefunctions */
        RT1 = new RmgTimer("Scf steps: Mg_eig");
        istop = kpoint->nstates / T->get_threads_per_node();
        istop = istop * T->get_threads_per_node();
        for(st1=0;st1 < istop;st1+=T->get_threads_per_node()) {
          SCF_THREAD_CONTROL thread_control[MAX_RMG_THREADS];
          for(ist = 0;ist < T->get_threads_per_node();ist++) {
              thread_control[ist].job = HYBRID_EIG;
              thread_control[ist].vtot = vtot_psi;
              thread_control[ist].sp = &states[st1 + ist];
              T->set_pptr(ist, &thread_control[ist]);
//              set_pptr(ist, &thread_control[ist]);
          }

          // Thread tasks are set up so run them
          T->run_thread_tasks(T->get_threads_per_node());

        }

        // Process any remaining states in serial fashion
        for(st1 = istop;st1 < kpoint->nstates;st1++) {
//            MgEigState<ScfType> (Rmg_G, Rmg_T, 
            mg_eig_state_driver (&states[st1], 0, vtot_psi, ct.mg_eig_precision);
        }
        delete(RT1);

    }


    time2 = my_crtc ();

    /*wavefunctions have changed, projectors have to be recalculated */
    RT1 = new RmgTimer("Scf steps: Beta x psi");
    betaxpsi (states);
    delete(RT1);



    /* Now we orthognalize and optionally do subspace diagonalization
     * In the gamma point case, orthogonalization is not required when doing subspace diagonalization
     * For non-gamma point we have to do first orthogonalization and then, optionally subspace diagonalization
     * the reason is for non-gamma subdiag is not coded to solve generalized eigenvalue problem, it can
     * only solve the regular eigenvalue problem and that requires that wavefunctions are orthogonal to start with.*/
    diag_this_step = (ct.diag && ct.scf_steps % ct.diag == 0 && ct.scf_steps < ct.end_diag);
#if GAMMA_PT
    /* do diagonalizations if requested, if not orthogonalize */
    if (diag_this_step) {
        RT1 = new RmgTimer("Scf steps: Diagonalization");
        subdiag_gamma (kpoint->Kstates, vh, vnuc, vxc);
        delete(RT1);
    }
    else {
        RT1 = new RmgTimer("Scf steps: Orthogonalization");
        ortho (states, 0);
        delete(RT1);
    }
#else

    RT1 = new RmgTimer("Scf steps: Orthogonalization");
//    ortho (kpoint->Kstates, kpt);
    ortho (kpoint->Kstates, 0);
    delete(RT1);
    
    if (diag_this_step)
    {
	/*Projectores need to be updated prior to subspace diagonalization*/
        RT1 = new RmgTimer("Scf steps: Beta x psi");
        betaxpsi (states);
        delete(RT1);
        
        RT1 = new RmgTimer("Scf steps: Diagonalization");
        subdiag_nongamma (kpoint->kstates, vh, vnuc, vxc);
        delete(RT1);

    }
#endif
    
    
    /*wavefunctions have changed, projectors have to be recalculated */
    time1 = my_crtc ();
    RT1 = new RmgTimer("Scf steps: Beta x psi");
    betaxpsi (states);
    delete(RT1);
    
    /*Get oldsintR*/
    if (diag_this_step)
	mix_betaxpsi(0);
    else 
	mix_betaxpsi(1);
    


    if (spin_flag)
	get_opposite_eigvals (states);

	

#if 0
    /* Take care of occupation filling */
    if  (!firststep)
	ct.efermi = fill (states, ct.occ_width, ct.nel, ct.occ_mix, kpoint->nstates, ct.occ_flag);





    if (ct.occ_flag == 1 && !firststep)
    {
        printf ("\n");
        //progress_tag ();
        printf ("FERMI ENERGY = %15.8f eV\n", ct.efermi * Ha_eV);
    }
#endif

    if (firststep)
        firststep = false;

    /* Generate new density */
    RT1 = new RmgTimer("Scf steps: Get rho");
    get_new_rho (states, new_rho);

    /*Takes care of mixing and checks whether the charge density is negative*/
    mix_rho(new_rho, rho, rhocore, FP0_BASIS, get_FPX0_GRID(), get_FPY0_GRID(), get_FPZ0_GRID());

    if (spin_flag)
	get_rho_oppo (rho,  rho_oppo);
    
    delete(RT1);


    /* If sorting is requested then sort the states. */
    if (ct.sortflag)
        kpoint->sort_orbitals();


    /* Make sure we see output, e.g. so we can shut down errant runs */
    fflush( ct.logfile );
	fsync( fileno(ct.logfile) );

    /* release memory */
    delete [] new_rho;
    delete [] vtot;
    delete [] vtot_psi;

    /* free the memory */
    if (spin_flag)
    {
    	delete [] rho_tot;
    }

    printf("\n SCF STEP TIME = %10.2f\n",my_crtc () - time3);

    return CONVERGED;
}                               /* end scf */


/******/
