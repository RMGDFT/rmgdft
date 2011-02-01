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


int static firststep = TRUE;

void scf (STATE * states, REAL * vxc, REAL * vh, REAL * vnuc,
          REAL * rho, REAL * rhocore, REAL * rhoc, int *CONVERGENCE)
{

    int kpt, st1, idx, ik;
    REAL t3;
    REAL *vtot, *vtot_psi, *new_rho;
    REAL time1, time2, time3;
    REAL t[3];                  /* SCF checks and average potential */


    time3 = my_crtc ();

    my_malloc (new_rho, FP0_BASIS, REAL);
    my_malloc (vtot, FP0_BASIS, REAL);
    my_malloc (vtot_psi, P0_BASIS, REAL);

    /* save old vhxc + vnuc */
    for (idx = 0; idx < FP0_BASIS; idx++)
        vtot[idx] = vxc[idx] + vh[idx] + vnuc[idx];


    time1 = my_crtc ();
    /* Generate exchange-correlation potential */
    get_vxc (rho, rhocore, vxc);
    rmg_timings (SCF_XC_TIME, (my_crtc () - time1), 0);


    /* Generate hartree potential */
    get_vh (rho, rhoc, vh, 15, ct.poi_parm.levels);

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
    global_sums (t, &idx);
    t[0] *= ct.vel_f;
    t[1] = sqrt (t[1] / ((REAL) (ct.psi_fnbasis)));
    t[2] /= ((REAL) (ct.psi_fnbasis));

    if (pct.thispe == 0 && !firststep)
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
        *CONVERGENCE = TRUE;

    get_vtot_psi (vtot_psi, vtot);

    /*Generate the Dnm_I */
    get_ddd (vtot);

#if MPI

    time1 = my_crtc ();
    /* Update the wavefunctions */
    for (st1 = 0; st1 < ct.num_kpts * ct.num_states; st1++)
        mg_eig_state (&states[st1], 0, vtot_psi);

    time2 = my_crtc ();
    rmg_timings (EIG_TIME, (time2 - time1), 0);

    /*wavefunctions have changed, projectors have to be recalculated */
    betaxpsi (states);






    /* Now we orthognalize and optionally do subspace diagonalization
     * In the gamma point case, orthogonalization is not required when doing subspace diagonalization
     * For non-gamma point we have to do first orthogonalization and then, optionally subspace diagonalization
     * the reason is for non-gamma subdiag is not coded to solve generalized eigenvalue problem, it can
     * only solve the regular eigenvalue problem and that requires that wavefunctions are orthogonal to start with.*/
#if GAMMA_PT
    /* do diagonalizations if requested, if not orthogonalize */
    if (ct.diag && ct.scf_steps % ct.diag == 0 && ct.scf_steps < ct.end_diag)
        subdiag_gamma (ct.kp[0].kstate, vh, vnuc, vxc);
    else
        ortho_full (states);
#else
    ortho_full (states);
    if (ct.diag && ct.scf_steps > 0 && ct.scf_steps % ct.diag == 0 && ct.scf_steps < ct.end_diag)
        for (ik = 0; ik < ct.num_kpts; ik++)
            subdiag_nongamma (ct.kp[ik].kstate, vh, vnuc, vxc);
#endif


    /* Take care of occupation filling */
    if (!firststep)
        ct.efermi = fill (states, ct.occ_width, ct.nel, ct.occ_mix, ct.num_states, ct.occ_flag);


#endif

    if (pct.thispe == 0 && ct.occ_flag == 1 && !firststep)
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
    mix_rho(new_rho, rho, rhocore, FP0_BASIS, FPX0_GRID, FPY0_GRID, FPZ0_GRID);

    time2 = my_crtc ();
    rmg_timings (RHO_TIME, (time2 - time1), 0);


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

    rmg_timings (SCF_TIME, (my_crtc () - time3), 0);

}                               /* end scf */


/******/
