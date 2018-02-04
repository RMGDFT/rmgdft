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



#include <complex>
#include "const.h"
#include "params.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "rmg_error.h"
#include "State.h"
#include "Kpoint.h"
#include "Functional.h"
#include "GlobalSums.h"
#include "transition.h"
#include "RmgSumAll.h"
#include "RmgParallelFft.h"

template void GetTe (double *, double *, double *, double *, double *, double *, Kpoint<double> **, int);
template void GetTe (double *, double *, double *, double *, double *, double *, Kpoint<std::complex<double> > **, int);

template <typename KpointType>
void GetTe (double * rho, double * rho_oppo, double * rhocore, double * rhoc, double * vh_in, double * vxc_in, Kpoint<KpointType> **Kptr, int ii_flag)
{
    int state, kpt, idx, nspin = (ct.spin_flag + 1), FP0_BASIS;
    double esum[3], t1, eigsum, xcstate, mag, absmag = 0.0;
    double vel;
    double *exc, *nrho, *nrho_oppo=NULL;
    Kpoint<KpointType> *kptr;

    FP0_BASIS = get_FP0_BASIS();

    double *vh = new double[FP0_BASIS];
    double *vxc = new double[FP0_BASIS];
    for(int i=0;i < FP0_BASIS;i++)vh[i] = vh_in[i];
    for(int i=0;i < FP0_BASIS;i++)vxc[i] = vxc_in[i];
    if(ct.filter_dpot) FftFilter(vh, *fine_pwaves, sqrt(ct.filter_factor) / (double)ct.FG_RATIO, LOW_PASS);
    if(ct.filter_dpot) FftFilter(vxc, *fine_pwaves, sqrt(ct.filter_factor) / (double)ct.FG_RATIO, LOW_PASS);


    vel = get_vel_f();

    /* Grab some memory */
    if (ct.spin_flag)
    {
        exc = new double[3 * FP0_BASIS];
    	nrho_oppo = exc + 2 * FP0_BASIS;
    }
    else
        exc = new double[2 * FP0_BASIS];
    
    nrho = exc + FP0_BASIS;


    /* Loop over states and get sum of the eigenvalues */
    eigsum = 0.0;

    for (idx = 0; idx < nspin; idx++)
    {
    	for (kpt = 0; kpt < ct.num_kpts_pe; kpt++)
    	{
            
            kptr = Kptr[kpt];
            t1 = 0.0;
            for (state = 0; state < ct.num_states; state++)
            {

                t1 += (kptr->Kstates[state].occupation[idx] *
                        kptr->Kstates[state].eig[idx]);

            }
            eigsum += t1 * kptr->kweight;
    	}
    }


    eigsum = RmgSumAll(eigsum, pct.kpsub_comm);
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
                absmag += fabs(rho[idx] - rho_oppo[idx]);
        }
    }


    for (idx = 0; idx < FP0_BASIS; idx++)
        esum[2] += rho[idx] * vxc[idx];



    /*Sum emergies over all processors */
    if (ct.spin_flag)
    {
    	GlobalSums (esum, 2, pct.grid_comm);
    	GlobalSums (&esum[2], 1, pct.img_comm);  
    	GlobalSums (&mag, 1, pct.grid_comm); 
    	GlobalSums (&absmag, 1, pct.grid_comm); 
    }
    else
    	GlobalSums (esum, 3, pct.grid_comm);

    /*Electrostatic E */
    ct.ES = 0.5 * vel * esum[0];

    /*XC potential energy */
    xcstate = vel * esum[2];
    mag *= vel;
    absmag *= vel;

    if(ii_flag) {

        /* Evaluate total ion-ion energy */
        ct.II = IonIonEnergy_Ewald();

    }


    /* Sum them all up */
    //ct.scf_correction = 0.0;
    //printf("TTT %12.8e  %12.8e  %12.8e  %12.8e  %12.8e  %12.8e\n",eigsum,ct.ES,xcstate,ct.XC,ct.II,ct.scf_correction); 
    ct.TOTAL = eigsum - ct.ES - xcstate + ct.XC + ct.II + ct.scf_correction;
    
    /* Print contributions to total energies into output file */
    double efactor = ct.energy_output_conversion[ct.energy_output_units];
    char *eunits = ct.energy_output_string[ct.energy_output_units];
    rmg_printf ("@@ EIGENVALUE SUM     = %15.6f %s\n", efactor*eigsum, eunits);
    rmg_printf ("@@ ION_ION            = %15.6f %s\n", efactor*ct.II, eunits);
    rmg_printf ("@@ ELECTROSTATIC      = %15.6f %s\n", -efactor*ct.ES, eunits);
    rmg_printf ("@@ VXC                = %15.6f %s\n",  efactor*xcstate, eunits);
    rmg_printf ("@@ EXC                = %15.6f %s\n", efactor*ct.XC, eunits);
    rmg_printf ("@@ TOTAL ENERGY       = %15.6f %s\n", efactor*ct.TOTAL, eunits);
    if(ct.scf_steps != 0) {
        rmg_printf ("@@ estimated error    =       %9.2e %s\n", efactor*ct.scf_accuracy, eunits);
    }
    else {
        rmg_printf ("@@ estimated error    =   ****************\n");
    }
        
    if (ct.spin_flag)
    {
	/* Print the total magetization and absolute magnetization into output file */
       	rmg_printf ("@@ TOTAL MAGNETIZATION    = %8.4f Bohr mag/cell\n", mag );
       	rmg_printf ("@@ ABSOLUTE MAGNETIZATION = %8.4f Bohr mag/cell\n", absmag );
    }

   

    /* Release our memory */
    delete [] exc;

    delete [] vxc;
    delete [] vh;

}                               /* end get_te */
