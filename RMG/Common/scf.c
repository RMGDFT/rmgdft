/************************** SVN Revision Information **************************
 **    $Id$    **
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
 *   void scf(STATE *states, REAL *vxc, REAL *vh, REAL *vnuc,
 *            REAL *rho, REAL *rhocore, REAL *rhoc, int *CONVERGENCE)
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
#include "main.h"

#if HYBRID_MODEL
#include "hybrid.h"
#include <pthread.h>
#endif


static int firststep = TRUE;

bool scf (STATE * states, REAL * vxc, REAL * vh, REAL * vnuc,
          REAL * rho, REAL * rho_oppo, REAL * rhocore, REAL * rhoc )
{

    int kpt, st1, idx, ik, st, diag_this_step, nspin = (ct.spin_flag + 1), istop, vcycle;
    bool CONVERGED = false;
    REAL t3;
    REAL *vtot, *vtot_psi, *new_rho;
    REAL time1, time2, time3;
    REAL t[3];                  /* SCF checks and average potential */
    int ist;

    MPI_Status status, stat[2]; 
    MPI_Request req[2];   
    /* to hold the send data and receive data of eigenvalues */
    REAL *eigval_sd, *eigval_rv, *rho_tot;   
    
    time3 = my_crtc ();

    /* allocate memory for eigenvalue send array and receive array */
    if (ct.spin_flag)
    {
    	my_malloc (rho_tot, pct.FP0_BASIS, REAL);
    }

    my_malloc (new_rho, pct.FP0_BASIS, REAL);
    my_malloc (vtot, pct.FP0_BASIS, REAL);
    my_malloc (vtot_psi,pct.P0_BASIS, REAL);

    /* save old vhxc + vnuc */
    for (idx = 0; idx < pct.FP0_BASIS; idx++)
        vtot[idx] = vxc[idx] + vh[idx] + vnuc[idx];


    time1 = my_crtc (); 
 
    /* Generate exchange-correlation potential */
    get_vxc (rho, rho_oppo, rhocore, vxc);
    rmg_timings (SCF_XC_TIME, (my_crtc () - time1));

    if (ct.spin_flag)        
    {
	/*calculate the total charge density in order to calculate hartree potential*/
	for (idx = 0; idx < pct.FP0_BASIS; idx++)
		rho_tot[idx] = rho[idx] + rho_oppo[idx];
	
	/* Generate hartree potential */
        get_vh (rho_tot, rhoc, vh, ct.hartree_min_sweeps, ct.hartree_max_sweeps, ct.poi_parm.levels, ct.rms/ct.hartree_rms_ratio);
     }  	
    else
    {
    	/* Generate hartree potential */
    	get_vh (rho, rhoc, vh, ct.hartree_min_sweeps, ct.hartree_max_sweeps, ct.poi_parm.levels, ct.rms/ct.hartree_rms_ratio);
    }



    /* check convergence */
    t[0] = t[1] = t[2] = 0.0;


    for (idx = 0; idx < pct.FP0_BASIS; idx++)
    {
        t3 = -vtot[idx];
        vtot[idx] = vxc[idx] + vh[idx] + vnuc[idx];
        t3 += vtot[idx];
        t[0] += rho[idx] * t3;
        t[1] += t3 * t3;
        t[2] += vh[idx];
    }                           /* idx */


    
    
    idx = 3;
    global_sums (t, &idx, pct.img_comm);
    t[0] *= ct.vel_f;
    
    /* get the averaged value over each spin and each fine grid */
    t[1] = sqrt (t[1] / ((REAL) (nspin * ct.psi_fnbasis)));  
    t[2] /= ((REAL) (nspin * ct.psi_fnbasis));   
    
    ct.rms = t[1];

    if (!firststep)
    {
        printf ("\n");
        progress_tag ();
        printf ("SCF CHECKS: <rho dv>  = %15.8e\n", t[0]);
        progress_tag ();
        printf ("SCF CHECKS: RMS[dv]   = %15.8e\n", t[1]);
        progress_tag ();
        printf ("AVERAGE POTENTIAL <V> = %15.8e\n", t[2]);
    }

    if (!firststep && t[1] < ct.thr_rms)
    {
	    CONVERGED = true;
    }

    get_vtot_psi (vtot_psi, vtot, FG_NX);

    /*Generate the Dnm_I */
    get_ddd (vtot);
#if MPI

    time1 = my_crtc ();
#if HYBRID_MODEL
    for(vcycle = 0;vcycle < ct.eig_parm.mucycles;vcycle++) {
        betaxpsi (states);
#if BATCH_NLS
        app_nls_batch (states, pct.nv, pct.ns, pct.oldsintR_local);
#endif

        enter_threaded_region();
        scf_barrier_init(ct.THREADS_PER_NODE);
        /* Update the wavefunctions */
        istop = ct.num_kpts * ct.num_states / ct.THREADS_PER_NODE;
        istop = istop * ct.THREADS_PER_NODE;
        for(st1=0;st1 < istop;st1+=ct.THREADS_PER_NODE) {
          for(ist = 0;ist < ct.THREADS_PER_NODE;ist++) {
              thread_control[ist].job = HYBRID_EIG;
              thread_control[ist].vtot = vtot_psi;
              thread_control[ist].sp = &states[st1 + ist];
          }

          // Thread tasks are set up so wake them
          wake_threads(ct.THREADS_PER_NODE);

          // Then wait for them to finish this task
          wait_for_threads(ct.THREADS_PER_NODE);
        }
        scf_barrier_destroy();
        leave_threaded_region();

        // Process any remaining states in serial fashion
        for(st1 = istop;st1 < ct.num_kpts * ct.num_states;st1++) {
            mg_eig_state_driver (&states[st1], 0, vtot_psi, ct.mg_eig_precision);
        }

    }

#else
    /* Update the wavefunctions */
    for(vcycle = 0;vcycle < ct.eig_parm.mucycles;vcycle++) {
        betaxpsi (states);
#if BATCH_NLS
        app_nls_batch (states, pct.nv, pct.ns, pct.oldsintR_local);
#endif
        for (st1 = 0; st1 < ct.num_kpts * ct.num_states; st1++) {
            mg_eig_state_driver (&states[st1], 0, vtot_psi, ct.mg_eig_precision);
        }
    }
#endif

    time2 = my_crtc ();
    rmg_timings (EIG_TIME, (time2 - time1));

    /*wavefunctions have changed, projectors have to be recalculated */
    time1 = my_crtc ();
    betaxpsi (states);
    rmg_timings (SCF_BETAXPSI, (my_crtc () - time1));






    /* Now we orthognalize and optionally do subspace diagonalization
     * In the gamma point case, orthogonalization is not required when doing subspace diagonalization
     * For non-gamma point we have to do first orthogonalization and then, optionally subspace diagonalization
     * the reason is for non-gamma subdiag is not coded to solve generalized eigenvalue problem, it can
     * only solve the regular eigenvalue problem and that requires that wavefunctions are orthogonal to start with.*/
    diag_this_step = (ct.diag && ct.scf_steps % ct.diag == 0 && ct.scf_steps < ct.end_diag);
#if GAMMA_PT
    /* do diagonalizations if requested, if not orthogonalize */
    if (diag_this_step)
        subdiag_gamma (ct.kp[0].kstate, vh, vnuc, vxc);
    else
        ortho (states, 0);
#else
    for (kpt =0; kpt < ct.num_kpts; kpt++)
        ortho (&states[kpt *ct.num_states], kpt);
    
    if (diag_this_step)
    {
	/*Projectores need to be updated prior to subspace diagonalization*/
        time1 = my_crtc ();
	
        betaxpsi (states);
        
        rmg_timings (SCF_BETAXPSI, (my_crtc () - time1));
        
        for (ik = 0; ik < ct.num_kpts; ik++)
            subdiag_nongamma (ct.kp[ik].kstate, vh, vnuc, vxc);
    }
#endif
    
    
    /*wavefunctions have changed, projectors have to be recalculated */
    time1 = my_crtc ();
    betaxpsi (states);
    
    /*Get oldsintR*/
    if (diag_this_step)
	mix_betaxpsi(0);
    else 
	mix_betaxpsi(1);
    
    rmg_timings (SCF_BETAXPSI, (my_crtc () - time1));


    if (ct.spin_flag)
	get_opposite_eigvals (states);

	


    /* Take care of occupation filling */
    if  (!firststep)
	ct.efermi = fill (states, ct.occ_width, ct.nel, ct.occ_mix, ct.num_states, ct.occ_flag);

#endif



    if (ct.occ_flag == 1 && !firststep)
    {
        printf ("\n");
        progress_tag ();
        printf ("FERMI ENERGY = %15.8f eV\n", ct.efermi * Ha_eV);
    }

    if (firststep)
        firststep = FALSE;

    /* Generate new density */
    time1 = my_crtc ();
    get_new_rho (states, new_rho);

    /*Takes care of mixing and checks whether the charge density is negative*/
    mix_rho(new_rho, rho, rhocore, pct.FP0_BASIS, pct.FPX0_GRID, pct.FPY0_GRID, pct.FPZ0_GRID);

    if (ct.spin_flag)
	get_rho_oppo (rho,  rho_oppo);
    
    time2 = my_crtc ();
    rmg_timings (RHO_TIME, (time2 - time1));


    /* If sorting is requested then sort the states. */
    for (kpt = 0; kpt < ct.num_kpts; kpt++)
        if (ct.sortflag)
            sortpsi (ct.kp[kpt].kstate);


    /* Make sure we see output, e.g. so we can shut down errant runs */
    fflush( ct.logfile );
	fsync( fileno(ct.logfile) );

    /* release memory */
    my_free (new_rho);
    my_free (vtot);
    my_free (vtot_psi);

    /* free the memory */
    if (ct.spin_flag)
    {
    	my_free (rho_tot);
    }

    rmg_timings (SCF_TIME, (my_crtc () - time3));
    printf("\n SCF STEP TIME = %10.2f\n",my_crtc () - time3);

    return CONVERGED;
}                               /* end scf */


/******/
