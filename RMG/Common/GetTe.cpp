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


double  vdw_d2_energy(Lattice &, std::vector<ION> &);
template void GetTe (double *, double *, double *, double *, double *, double *, Kpoint<double> **, int);
template void GetTe (double *, double *, double *, double *, double *, double *, Kpoint<std::complex<double> > **, int);

extern double *vnuc;

template <typename KpointType>
void GetTe (double * rho, double * rho_oppo, double * rhocore, double * rhoc, double * vh_in, double * vxc_in, Kpoint<KpointType> **Kptr, int ii_flag)
{
    int state, kpt, idx, FP0_BASIS, P0_BASIS;
    double t1;
    double vel, cvel, ES_pa = 0.0;
    bool potential_acceleration = (ct.potential_acceleration_constant_step > 0.0);
    if(Verify ("kohn_sham_solver","davidson", Kptr[0]->ControlMap)) potential_acceleration = false;
    Kpoint<KpointType> *kptr;

    FP0_BASIS = get_FP0_BASIS();
    P0_BASIS = get_P0_BASIS();

    double *vh = new double[FP0_BASIS];
    double *vxc = new double[ct.nspin * FP0_BASIS];
    double *vxc_up = vxc;
    double *vxc_down = vxc + FP0_BASIS;

    for(int i=0;i < FP0_BASIS;i++)vh[i] = vh_in[i];
    for(int i=0;i < ct.nspin*FP0_BASIS;i++)vxc[i] = vxc_in[i];

    vel = get_vel_f();
    cvel = get_vel();

    /* Loop over states and get sum of the eigenvalues and any LDA+U corrections */
    double eigsum = 0.0;
    double ldaU_E = 0.0;
    double ldaU_H = 0.0;

    int nspin = 1;
    if(ct.nspin == 2 && !ct.AFM) nspin = 2;
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
            eigsum += t1 * kptr->kp.kweight;
        }
    }

    eigsum = RmgSumAll(eigsum, pct.kpsub_comm);
    if(ct.AFM) eigsum *= 2.0;

    if((ct.ldaU_mode != LDA_PLUS_U_NONE) && (ct.num_ldaU_ions > 0))
    {
        ldaU_E = Kptr[0]->ldaU->Ecorrect;
        ldaU_H = Kptr[0]->ldaU->Ehub;
    }

    /* Evaluate electrostatic energy correction terms */
    ct.ES = 0.0;
    if (ct.nspin==2)
    {
        /* Add the compensating charge to total charge to calculation electrostatic energy */    
        for (idx = 0; idx < FP0_BASIS; idx++)
            ct.ES += (rho[idx] + rho_oppo[idx] + rhoc[idx]) * vh[idx];

    }
    else 
    {
        for (idx = 0; idx < FP0_BASIS; idx++)
            ct.ES += (rho[idx] + rhoc[idx]) * vh[idx];
    }
    ct.ES = 0.5 * vel * RmgSumAll(ct.ES, pct.grid_comm);

    if(potential_acceleration)
    {
        double *vf_tmp = new double[FP0_BASIS]();
        double *vc_tmp = new double[P0_BASIS]();

        for(int ix = 0;ix < FP0_BASIS;ix++) vf_tmp[ix] = vxc[ix] + vnuc[ix] + vh[ix];
        GetVtotPsi (vc_tmp, vf_tmp, Rmg_G->default_FG_RATIO);
        for (kpt = 0; kpt < ct.num_kpts_pe; kpt++)
        {
            kptr = Kptr[kpt];
            int my_pe_x, my_pe_y, my_pe_z;
            kptr->G->pe2xyz(pct.gridpe, &my_pe_x, &my_pe_y, &my_pe_z);
            int my_pe_offset = my_pe_x % pct.coalesce_factor;

            for (int is = 0; is < ct.num_states; is++)
            {
                int offset = (is / kptr->dvh_skip) * kptr->pbasis;
                double *veff = &kptr->dvh[offset*pct.coalesce_factor + my_pe_offset*kptr->pbasis];
                KpointType *psi = kptr->Kstates[is].psi;
                for(int ix = 0;ix < kptr->pbasis;ix++)
                {
                    ES_pa += std::real(psi[ix] * std::conj(psi[ix])) * kptr->kp.kweight *
                             kptr->Kstates[is].occupation[0] * (veff[ix] - vc_tmp[ix]);
                }
            }

        }
        delete [] vc_tmp;
        delete [] vf_tmp;
        ES_pa = 0.5 * cvel * RmgSumAll(ES_pa, pct.grid_comm);
        ES_pa = RmgSumAll(ES_pa, pct.spin_comm);
        ES_pa = RmgSumAll(ES_pa, pct.kpsub_comm);

    }

    ct.ES += ES_pa;

    double mag = 0.0;    
    double absmag = 0.0;
    double xcstate = 0.0;

    if (ct.nspin == 2)
    {
        for (idx = 0; idx < FP0_BASIS; idx++)
        {
            xcstate += rho[idx]*vxc_up[idx] + rho_oppo[idx]*vxc_down[idx];
            mag += ( rho[idx] - rho_oppo[idx] );       /* calculate the magnetization */
            absmag += fabs(rho[idx] - rho_oppo[idx]);
        }
        mag = vel * RmgSumAll(mag, pct.grid_comm);
        absmag = vel * RmgSumAll(absmag, pct.grid_comm);
    }
    else if(ct.nspin == 4)
    {
        for(int is = 0; is < 4; is++)
            for (idx = 0; idx < FP0_BASIS; idx++)
                xcstate += rho[idx + is* FP0_BASIS] * vxc_in[idx + is*FP0_BASIS];
    }
    else
    {
        //for (idx = 0; idx < FP0_BASIS; idx++) xcstate += rho[idx] * vxc_in[idx];
        for (idx = 0; idx < FP0_BASIS; idx++) xcstate += rho[idx] * vxc[idx];
        //for (idx = 0; idx < FP0_BASIS; idx++) xcstate += rhovxc[idx] * vxc[idx];
    }

    /*XC potential energy */
    xcstate = vel * RmgSumAll(xcstate, pct.grid_comm);


    if(ii_flag) {

        /* Evaluate total ion-ion energy */
        ct.II = IonIonEnergy_Ewald();
        ct.Evdw = 0.0;
        if(ct.force_vdw == NULL) ct.force_vdw = new double[3 * ct.num_ions];

        if(ct.vdw_corr == DFT_D2)
        {
            ct.Evdw = vdw_d2_energy(Rmg_L, Atoms);
        }

        if(ct.vdw_corr == DFT_D3)
        {
            RmgDftd3(&ct.Evdw, ct.force_vdw, ct.stress_vdw, ct.dftd3_version);
        }

    }


    /* Sum them all up */
    ct.ldaU_E = ldaU_E;
    ct.TOTAL = eigsum - ct.ES - xcstate + ct.XC + ct.II + ct.ldaU_E + ct.scf_correction + ct.Evdw;
    if(ct.xc_is_hybrid && Functional::is_exx_active()) ct.TOTAL -= ct.FOCK;
    // AFM case requires counting FOCK energy twice
    if(ct.xc_is_hybrid && Functional::is_exx_active() && ct.AFM) ct.TOTAL -= ct.FOCK;

    ct.xcstate = xcstate;
    /* Print contributions to total energies into output file */
    double efactor = ct.energy_output_conversion[ct.energy_output_units];
    const char *eunits = ct.energy_output_string[ct.energy_output_units].c_str();
    //rmg_printf ("@@ SCF CORRECTION     = %15.6f %s\n", efactor*ct.scf_correction, eunits);
    rmg_printf ("@@ EIGENVALUE SUM     = %15.6f %s\n", efactor*eigsum, eunits);
    rmg_printf ("@@ ION_ION            = %15.6f %s\n", efactor*ct.II, eunits);
    rmg_printf ("@@ ELECTROSTATIC      = %15.6f %s\n", -efactor*ct.ES, eunits);
    rmg_printf ("@@ VXC                = %15.6f %s\n",  efactor*xcstate, eunits);
    rmg_printf ("@@ EXC                = %15.6f %s\n", efactor*ct.XC, eunits);
    if(ct.vdw_corr)
        rmg_printf ("@@ vdw_corr           = %15.6f %s\n", efactor*ct.Evdw, eunits);
    if(ct.xc_is_hybrid && Functional::is_exx_active())
        rmg_printf ("@@ FOCK               = %15.6f %s\n", efactor*ct.FOCK, eunits);

    if((ct.ldaU_mode != LDA_PLUS_U_NONE) && (ct.num_ldaU_ions > 0))
        rmg_printf ("@@ HUBBARD ENERGY     = %15.6f %s\n", efactor*ldaU_H, eunits);

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
        rmg_printf ("@@ TOTAL MAGNETIZATION    = %12.8f Bohr mag/cell\n", mag );
        rmg_printf ("@@ ABSOLUTE MAGNETIZATION = %12.8f Bohr mag/cell\n", absmag );
    }

    //rmg_printf("CHECK  %12.6f  %12.6f\n", xcstate, ct.vtxc);


    /* Release our memory */
    delete [] vxc;
    delete [] vh;

}                               /* end get_te */
