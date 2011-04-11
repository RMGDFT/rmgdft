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
 *   rho:  total charge density in spin-pairwised calculation; while in spin 
 *         polarized calculation, it's the processor's own spin charge density
 *   rho_oppo: the charge density for the opposite spin      
 *   rhocore: charge density of core electrons, only useful when we 
 *            include non-linear core correction for pseudopotential.
 *   rhoc:    compensating charge density
 *   vh:  Hartree potential
 *   vxc: exchange-correlation potential for the processor's own spin
 *   states: point to orbital structure
 * OUTPUT
 *   total energy is printed out
 * PARENTS
 *   cdfastrlx.c fastrlx.c moldyn.c quench.c
 * CHILDREN
 *   get_xc.c
 * SOURCE
 */



#include <math.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include "main.h"


void get_te (REAL * rho, REAL * rho_oppo, REAL * rhocore, REAL * rhoc, REAL * vh, REAL * vxc, STATE * states)
{
    int state, kpt, idx, i, j, three = 3, two = 2, one = 1, nspin = (pct.spin_flag + 1);
    REAL r, esum[3], t1, eigsum, xcstate, xtal_r[3], mag;
    REAL vel;
    REAL *exc, *nrho, *nrho_oppo;
    ION *iptr1, *iptr2;
    REAL time1, time2;

    time1 = my_crtc ();

    vel = ct.vel_f;

    /* Grab some memory */
    if (pct.spin_flag)
    {
    	my_malloc (exc, 3 * FP0_BASIS, REAL);
    	nrho_oppo = exc + 2 * FP0_BASIS;
    }
    else
    	my_malloc (exc, 2 * FP0_BASIS, REAL);
    
    nrho = exc + FP0_BASIS;


    /* Loop over states and get sum of the eigenvalues */
    eigsum = 0.0;

    for (idx = 0; idx < nspin; idx++)
    {
    	for (kpt = 0; kpt < ct.num_kpts; kpt++)
    	{
        	t1 = 0.0;
        	for (state = 0; state < ct.num_states; state++)
        	{

            		t1 += (states[state + kpt * ct.num_states].occupation[idx] *
                   		states[state + kpt * ct.num_states].eig[idx]);

        	}
        	eigsum += t1 * ct.kp[kpt].kweight;
    	}
    }


    /* Evaluate electrostatic energy correction terms */
    esum[0] = 0.0;
    if (pct.spin_flag)
    {
	/* Add the compensating charge to total charge to calculation electrostatic energy */    
    	for (idx = 0; idx < FP0_BASIS; idx++)
	    	esum[0] += (rho[idx] + rho_oppo[idx] + rhoc[idx]) * vh[idx];

    }
    else 
    {
    	for (idx = 0; idx < FP0_BASIS; idx++)
        	esum[0] += (rho[idx] + rhoc[idx]) * vh[idx];
    }


    time2 = my_crtc (); 

    /* Add the nonlinear core correction charge if there is any */
    if (pct.spin_flag)
    {
    	for (idx = 0; idx < FP0_BASIS; idx++)
    	{    
	    	nrho[idx] = rhocore[idx] * 0.5 + rho[idx];
	    	nrho_oppo[idx] = rhocore[idx] * 0.5 + rho_oppo[idx];
    	}
    }
    else
    {
    	for (idx = 0; idx < FP0_BASIS; idx++)
        	nrho[idx] = rhocore[idx] + rho[idx];
    }


    /* Evaluate XC energy and potential */
    get_xc(nrho, nrho_oppo, vxc, exc, ct.xctype);


    esum[1] = 0.0;
    esum[2] = 0.0;

    if (pct.spin_flag)
    {
	mag = 0.0;    
    	for (idx = 0; idx < FP0_BASIS; idx++)
	{
        	esum[1] += (rho[idx] + rho_oppo[idx] + rhocore[idx]) * (exc[idx]);
		mag += ( rho[idx] - rho_oppo[idx] );       /* calculation the magnetization */
        }
    }
    else
    {
    	for (idx = 0; idx < FP0_BASIS; idx++)
        	esum[1] += (rhocore[idx] + rho[idx]) * exc[idx];
    }


    for (idx = 0; idx < FP0_BASIS; idx++)
        esum[2] += rho[idx] * vxc[idx];

    rmg_timings (GET_TE_XC_TIME, (my_crtc () - time2), 0);


    /*Sum emergies over all processors */
    if (pct.spin_flag)
    {
    	global_sums (esum, &two, pct.grid_comm);
    	global_sums (&esum[2], &one, pct.img_comm);  
    	global_sums (&mag, &one, pct.grid_comm); 
    }
    else
    	global_sums (esum, &three, pct.grid_comm);


    /*Electrostatic E */
    ct.ES = 0.5 * vel * esum[0];

    /* XC E */
    ct.XC = vel * esum[1];


    /*XC potential energy */
    xcstate = vel * esum[2];
    mag *= vel;


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
   
    
    /* Print contributions to total energies into output file */
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
        
    if (pct.spin_flag)
    {
	/* Print the total magetization and absolute magnetization into output file */
	progress_tag ();
       	printf ("@@ TOTAL MAGNETIZATION    = %8.4f Bohr mag/cell\n", mag );
       	progress_tag ();
       	printf ("@@ ABSOLUTE MAGNETIZATION = %8.4f Bohr mag/cell\n", fabs(mag) );
    }

   

    /* Release our memory */
    my_free (exc);

    rmg_timings (GET_TE_TIME, (my_crtc () - time1), 0);

}                               /* end get_te */
