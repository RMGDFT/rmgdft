/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

/****f* QMD-MGDFT/get_te.c *****
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
 *   void get_te(REAL *rho, REAL *rhocore, REAL *rhoc, REAL *vh, REAL *vxc,
 *               STATE *states)
 *   Gets total energy of the system. Stores result in control structure.
 * INPUTS
 *   rho:  total charge density
 *   rhocore: charge density of core electrons, only useful when we 
 *            include non-linear core correction for pseudopotential.
 *   rhoc:    compensating charge density
 *   vh:  Hartree potential
 *   vxc: exchange-correlation potential
 *   states: point to orbital structure
 * OUTPUT
 *   total energy is printed out
 * PARENTS
 *   cdfastrlx.c fastrlx.c moldyn.c quench.c
 * CHILDREN
 *   exclda_pz81.c xcgga.c
 * SOURCE
 */



#include <math.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include "main.h"


void get_te (REAL * rho, REAL * rhocore, REAL * rhoc, REAL * vh, REAL * vxc, STATE * states)
{
    int state, kpt, idx, i, j, three = 3;
    REAL r, esum[3], t1, eigsum, xcstate, xtal_r[3];
    REAL vel;
    REAL *exc, *nrho, *nrho_buff;
    ION *iptr1, *iptr2;
    REAL time1, time2;

    time1 = my_crtc ();

    vel = ct.vel_f;

    /* Grab some memory */
    my_malloc (exc, 2 * FP0_BASIS, REAL);
    nrho = exc + FP0_BASIS;


    /* Loop over states and get sum of the eigenvalues */
    eigsum = 0.0;
    for (kpt = 0; kpt < ct.num_kpts; kpt++)
    {
        t1 = 0.0;
        for (state = 0; state < ct.num_states; state++)
        {

            t1 += (states[state + kpt * ct.num_states].occupation *
                   states[state + kpt * ct.num_states].eig);

        }
        eigsum += t1 * ct.kp[kpt].kweight;
    }


    /* Evaluate electrostatic energy correction terms */
    esum[0] = 0.0;
    for (idx = 0; idx < FP0_BASIS; idx++)
        esum[0] += (rho[idx] + rhoc[idx]) * vh[idx];


    time2 = my_crtc ();
    for (idx = 0; idx < FP0_BASIS; idx++)
    {
        nrho[idx] = rhocore[idx] + rho[idx];
    }

    
 
    /* Evaluate XC energy correction terms */
    switch (ct.xctype)
    {

    case LDA_PZ81:             /* LDA Perdew Zunger 81 */
        /* exclda_pz81 (nrho, exc); */

	/* incoporate both the Perdew Zunger 1981 and Ortiz Ballone 1994, default is PZ 1981 */
        xclda (nrho, vxc, exc);
        break;

    case GGA_BLYP:             /* GGA X-Becke C-Lee Yang Parr */
        xcgga (nrho, vxc, exc, ct.xctype);
        break;

    case GGA_XB_CP:            /* GGA X-Becke C-Perdew */
        xcgga (nrho, vxc, exc, ct.xctype);
        break;

    case GGA_XP_CP:            /* GGA X-Perdew C-Perdew */
        xcgga (nrho, vxc, exc, ct.xctype);
        break;

    case GGA_PBE:              /* GGA Perdew, Burke, Ernzerhof */
        xcgga (nrho, vxc, exc, ct.xctype);
	break;

    default:
        error_handler ("Unknown exchange-correlation functional");

    }                           /* end switch */


    esum[1] = 0.0;
    esum[2] = 0.0;

    for (idx = 0; idx < FP0_BASIS; idx++)
    {
        esum[1] += (rhocore[idx] + rho[idx]) * exc[idx];
        esum[2] += rho[idx] * vxc[idx];
    }

    rmg_timings (GET_TE_XC_TIME, (my_crtc () - time2), 0);


    /*Sum emergies over all processors */
    global_sums (esum, &three);


    /*Electrostatic E */
    ct.ES = 0.5 * vel * esum[0];

    /* XC E */
    ct.XC = vel * esum[1];


    /*XC potential energy */
    xcstate = vel * esum[2];


    time2 = my_crtc ();
    /* Evaluate total ion-ion energy */
    ct.II = 0.0;
    for (i = 0; i < ct.num_ions; i++)
        ct.II -= (ct.sp[ct.ions[i].species].zvalence *
                  ct.sp[ct.ions[i].species].zvalence /
                  ct.sp[ct.ions[i].species].rc) / sqrt (2.0 * PI);


    for (i = 0; i < ct.num_ions; i++)
    {

        iptr1 = &ct.ions[i];
        for (j = i + 1; j < ct.num_ions; j++)
        {

            iptr2 = &ct.ions[j];

            r = minimage (iptr1, iptr2, xtal_r);

            t1 = sqrt (ct.sp[iptr1->species].rc * ct.sp[iptr1->species].rc +
                       ct.sp[iptr2->species].rc * ct.sp[iptr2->species].rc);

            ct.II += (ct.sp[iptr1->species].zvalence *
                      ct.sp[iptr2->species].zvalence * erfc (r / t1) / r);
        }
    }

    rmg_timings (GET_TE_II_TIME, (my_crtc () - time2), 0);


    /* Sum them all up */
    ct.TOTAL = eigsum - ct.ES - xcstate + ct.XC + ct.II;
    if (pct.thispe == 0)
    {
        printf ("\n\n");

        progress_tag ();
        printf ("@@ EIGENVALUE SUM     = %16.9f Ha\n", eigsum);
        progress_tag ();
        printf ("@@ ION_ION            = %16.9f Ha\n", ct.II);
        progress_tag ();
        printf ("@@ ELECTROSTATIC      = %16.9f Ha\n", -ct.ES);
        progress_tag ();
        printf ("@@ XC                 = %16.9f Ha\n", ct.XC - xcstate);
        progress_tag ();
        printf ("@@ TOTAL ENERGY       = %16.9f Ha\n", ct.TOTAL);

    }


    /* Release our memory */
    my_free (exc);

    rmg_timings (GET_TE_TIME, (my_crtc () - time1), 0);

}                               /* end get_te */
