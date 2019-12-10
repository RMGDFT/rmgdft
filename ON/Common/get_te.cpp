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
 *   void get_te(double *rho, double *rhocore, double *rhoc, double *vh, double *vxc,
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

#include "common_prototypes.h"
#include "main.h"
#include "transition.h"
#include "prototypes_on.h"
#include "LdaU_on.h"
#include "Functional.h"


void get_te (double * rho, double * rho_oppo, double * rhocore, double * rhoc, double * vh, double * vxc, STATE * states, int ii_flag)
{
    int state, kpt, idx, three = 3, two = 2, one = 1, nspin = (ct.spin_flag + 1), FP0_BASIS;
    double esum[3], t1, eigsum, xcstate, mag, absmag;
    double vel;
    double *exc, *nrho, *nrho_oppo=NULL;


    FP0_BASIS = get_FP0_BASIS();

    vel = get_vel_f();

    /* Grab some memory */
    if (ct.spin_flag)
    {
    	my_malloc (exc, 3 * FP0_BASIS, double);
    	nrho_oppo = exc + 2 * FP0_BASIS;
    }
    else
    	my_malloc (exc, 2 * FP0_BASIS, double);
    
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
    if (ct.spin_flag)
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



    /* Add the nonlinear core correction charge if there is any */
    if (ct.spin_flag)
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
    //get_vxc_exc(nrho, nrho_oppo, vxc, exc, ct.xctype);


    esum[1] = 0.0;
    esum[2] = 0.0;

    if (ct.spin_flag)
    {
        mag = 0.0;    
        absmag = 0.0;    
        for (idx = 0; idx < FP0_BASIS; idx++)
        {
            esum[1] += (rho[idx] + rho_oppo[idx] + rhocore[idx]) * (exc[idx]);
            mag += ( rho[idx] - rho_oppo[idx] );       /* calculation the magnetization */
            absmag += fabs( rho[idx] - rho_oppo[idx] );       /* calculation the magnetization */
        }
    }
    else
    {
        for (idx = 0; idx < FP0_BASIS; idx++)
            esum[1] += (rhocore[idx] + rho[idx]) * exc[idx];
    }


    for (idx = 0; idx < FP0_BASIS; idx++)
        esum[2] += rho[idx] * vxc[idx];



    /*Sum emergies over all processors */
    if (ct.spin_flag)
    {
        global_sums (esum, &two, pct.grid_comm);
        global_sums (&esum[2], &one, pct.img_comm);  
        global_sums (&mag, &one, pct.grid_comm); 
        global_sums (&absmag, &one, pct.grid_comm); 
    }
    else
        global_sums (esum, &three, pct.grid_comm);


    /*Electrostatic E */
    ct.ES = 0.5 * vel * esum[0];

    /* XC E */
    //ct.XC = vel * esum[1];


    /*XC potential energy */
    xcstate = vel * esum[2];
    mag *= vel;
    absmag *= vel;

    if(ii_flag) {
        ct.II = IonIonEnergy_Ewald();
    }


    /* Sum them all up */
    ct.TOTAL = eigsum - ct.ES - xcstate + ct.XC + ct.II;

    if(ct.xc_is_hybrid && Functional::is_exx_active()) ct.TOTAL -= ct.FOCK;


    if((ct.ldaU_mode != LDA_PLUS_U_NONE) && (ct.num_ldaU_ions > 0))
    {
        ct.TOTAL += ldaU_on->Ecorrect;
    }

    /* Print contributions to total energies into output file */
    printf ("\n\n");
    printf ("@@ EIGENVALUE SUM     = %16.9f Ha\n", eigsum);
    printf ("@@ ION_ION            = %16.9f Ha\n", ct.II);
    printf ("@@ ELECTROSTATIC      = %16.9f Ha\n", -ct.ES);
    printf ("@@ VXC                 = %16.9f Ha\n",  xcstate);
    printf ("@@ EXC                 = %16.9f Ha\n", ct.XC );
    if(ct.xc_is_hybrid && Functional::is_exx_active())
        rmg_printf ("@@ FOCK               = %16.9f Ha\n", ct.FOCK);
    if((ct.ldaU_mode != LDA_PLUS_U_NONE) && (ct.num_ldaU_ions > 0))
    {
        printf ("@@ HUBBARD ENERGY     = %15.6f Ha\n", ldaU_on->Ehub);
    }

    printf ("@@ TOTAL ENERGY       = %16.9f Ha\n", ct.TOTAL);

    if (ct.spin_flag)
    {
        /* Print the total magetization and absolute magnetization into output file */
        progress_tag ();
        printf ("@@ TOTAL MAGNETIZATION    = %8.4f Bohr mag/cell\n", mag );
        progress_tag ();
        printf ("@@ ABSOLUTE MAGNETIZATION = %8.4f Bohr mag/cell\n", absmag );
    }



    /* Release our memory */
    my_free (exc);


}                               /* end get_te */
