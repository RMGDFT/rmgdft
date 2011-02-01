/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

/****f* QMD-MGDFT/quench.c *****
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
 *   void quench(STATE *states, REAL *vxc, REAL *vh, REAL *vnuc, 
 *               REAL *rho, REAL *rhocore, REAL *rhoc)
 *   For a fixed atomic configuration, quench the electrons to find 
 *   the minimum total energy 
 * INPUTS
 *   states: point to orbital structure (see main.h)
 *   vxc:    exchange correlation potential
 *   vh:     Hartree potential
 *   vnuc:   Pseudopotential 
 *   rho:    total valence charge density
 *   rhocore: core chare density only for non-linear core correction
 *   rhoc:   Gaussian compensating charge density
 * OUTPUT
 *   states, vxc, vh, rho are updated
 * PARENTS
 *   cdfastrlx.c fastrlx.c main.c
 * CHILDREN
 *   scf.c force.c get_te.c subdiag.c get_ke.c
 * SOURCE
 */



#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "main.h"

void quench (STATE * states, REAL * vxc, REAL * vh, REAL * vnuc, REAL * rho,
             REAL * rhocore, REAL * rhoc)
{

    static int CONVERGENCE;
    int numacc = 1, ic;
    /*int ist, ik;
       REAL KE; */

    /* ---------- begin scf loop ---------- */


    for (ct.scf_steps = 0, CONVERGENCE = FALSE;
         ct.scf_steps < ct.max_scf_steps && !CONVERGENCE; ct.scf_steps++, ct.total_scf_steps++)
    {

        if (pct.thispe == 0)
            printf ("\n\nquench: ------ [md: %d/%d  scf: %d/%d] ------\n",
                    ct.md_steps, ct.max_md_steps, ct.scf_steps, ct.max_scf_steps);


        /* perform a single self-consistent step */
        scf (states, vxc, vh, vnuc, rho, rhocore, rhoc, &CONVERGENCE);



        /* ??? */
        if (ct.scf_steps == 0)
            CONVERGENCE = FALSE;        /* I guess just in case (?) */


        /* check if we need to output intermediate results */
        if (ct.outcount == 0 || (ct.scf_steps % ct.outcount) == 0)
        {
            /* get the total energy */
            get_te (rho, rhocore, rhoc, vh, vxc, states);
        }


        /* output the eigenvalues with occupations */
        if (ct.write_eigvals_period)
        {
            if (ct.scf_steps % ct.write_eigvals_period == 0)
            {
                if (pct.thispe == 0)
                {
                    output_eigenvalues (states, 0, ct.scf_steps);
                    printf ("\nTotal charge in supercell = %16.8f\n", ct.tcharge);
                }
            }
        }

        /* do diagonalizations if requested */
        /*This is now done in scf.c */
#if 0
        if (ct.diag &&
            ct.scf_steps > 0 &&
            ct.scf_steps % ct.diag == 0 && ct.scf_steps < ct.end_diag && !CONVERGENCE)
            for (ik = 0; ik < ct.num_kpts; ik++)
#if GAMMA_PT
                subdiag_gamma (ct.kp[ik].kstate, vh, vnuc, vxc);
#else
                subdiag_nongamma (ct.kp[ik].kstate, vh, vnuc, vxc);
#endif
#endif



    }

    /* ---------- end scf loop ---------- */

    if (pct.thispe == 0)
    {
        if (CONVERGENCE)
        {
            printf ("\n");
            progress_tag ();
            printf ("potential convergence has been achieved. stopping ...\n");
        }

        printf ("\n");
        progress_tag ();
        printf ("final total energy = %14.7f Ha\n", ct.TOTAL);
    }



    /* output final eigenvalues with occupations */
    if (pct.thispe == 0)
    {
        output_eigenvalues (states, 0, ct.scf_steps);
        printf ("\nTotal charge in supercell = %16.8f\n", ct.tcharge);
    }

    wvfn_residual (states);




    /*When running MD, force pointers need to be rotated before calculating new forces */
    if ((ct.forceflag == MD_CVE) || (ct.forceflag == MD_CVT) || (ct.forceflag == MD_CPT))
    {

        /* rotate the force pointers */
        switch (ct.mdorder)
        {
        case ORDER_2:
            numacc = 1;
            break;
        case ORDER_3:
            numacc = 2;
            break;
        case ORDER_5:
            numacc = 4;
            break;
        }
        for (ic = (numacc - 1); ic > 0; ic--)
        {
            ct.fpt[ic] = ct.fpt[ic - 1];
        }
        ct.fpt[0] = ct.fpt[numacc - 1] + 1;
        if (ct.fpt[0] > (numacc - 1) || numacc == 1)
            ct.fpt[0] = 0;

    }


    /* compute the forces */
    /* Do not calculate forces for quenching when we are not converged */
    if ((CONVERGENCE) || (ct.forceflag != MD_QUENCH))
        force (rho, rhoc, vh, vxc, vnuc, states);

    /* output the forces */
    if (pct.thispe == 0)
        write_force ();



#if 0
    /* compute kinetic energy of all states */
    if (pct.thispe == 0)
        printf ("\n");
    for (ist = 0; ist < ct.num_states; ist++)
    {
        KE = get_ke (&states[ist], 0);
        if (pct.thispe == 0)
        {
            progress_tag ();
            printf ("kinetic energy for state %3d = %14.6f\n", ist, KE);
        }
    }
#endif



}                               /* end quench */




/******/
