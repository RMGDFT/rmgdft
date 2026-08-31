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


#include <math.h>
#include <float.h>
#include <complex>
#include <iostream>
#include <fstream>
#include <fcntl.h>
#include <sys/stat.h>

#include "const.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "RmgException.h"
#include "rmg_sum_all.h"
#include "Atomic.h"
#include "transition.h"
#include "prototypes_tddft.h"
#include "rmg_tddft.h"
#include "rmg_ortho.h"


/* Local function prototypes */
void center_of_mass_velocity(double &vx, double &vy, double &vz);

void rms_disp (double *, double *);


template void EhrenfestDynamics<std::complex<double> >(Kpoint<std::complex<double>> **, spinobj<double> &vxc,
             fgobj<double> &vh, fgobj<double> &vnuc, spinobj<double> &rho, fgobj<double> &rhoc, fgobj<double> &rhocore);

template void EhrenfestDynamics (Kpoint<double> **Kptr, spinobj<double> &vxc, fgobj<double> &vh,
                        fgobj<double> &vnuc, spinobj<double> &rho, fgobj<double> &rhoc, fgobj<double> &rhocore);


template <typename KpointType>
void EhrenfestDynamics (Kpoint<KpointType> **Kptr, spinobj<double> &vxc, fgobj<double> &vh, fgobj<double> &vnuc,
             spinobj<double> &rho, fgobj<double> &rhoc, fgobj<double> &rhocore)
{
    double rms[3], trms;
    double nosekin, nosepot;
    double iontemp;
    int N;
    std::vector<double> RMSdV;
    static double *rhodiff;

    // Save T=0 basis
    for(int kpt=0;kpt < ct.num_kpts_pe;kpt++) Kptr[kpt]->save_wavefunctions();

    // If tddft dynamics create tddft object
    rmg::tddft<KpointType, KpointType> *tddftobj = NULL;
    rmg::tddft<KpointType, std::complex<double>> *tddftobj_vp = NULL;
    if(ct.tddft_mode == VECTOR_POT)
    {
        tddftobj_vp = new rmg::tddft<KpointType, std::complex<double>>(vxc, vh, vnuc, rho, rhocore, rhoc, Kptr);
    }
    else
    {
        tddftobj = new rmg::tddft<KpointType, KpointType>(vxc, vh, vnuc, rho, rhocore, rhoc, Kptr);
    }


    if(ct.tddft_mode == VECTOR_POT)
    {
        tddftobj_vp->tddft_md();
    }
    else
    {
        tddftobj->tddft_md();
    }
    //        WriteRestart (ct.outfile, vh.data(), rho.up.data(), rho.dw.data(), vxc.data(), Kptr);


    /*Get some memory */




    rmg::printlog ("\n ==============================================");

    /* print out the title */
    rmg::printlog ("\n TDDFT Integration with Beeman");


    ct.ionke = 0.0;

    /* define Boltzmann param */
    double kB = 1.0 / (11605.0 * Ha_eV);

    /* count up moving atoms */
    N = 0;
    for (size_t ion = 0, i_end = Atoms.size(); ion < i_end; ++ion)
    {
        if (Atoms[ion].movable[0] || Atoms[ion].movable[1] || Atoms[ion].movable[2] ) N++;
    }
    ct.nose.N = N;

    /* check to see if random velocities are needed */
    if (ct.nose.randomvel && (ct.runflag != RESTART))
    {
        if (pct.gridpe == 0)
        {
            rmg::printlog ("\n\n Initializing temperature to %14.10f K\n", ct.nose.temp);
        }
        ranv ();
    }


    /* zero out the non-moving atoms */
    allatoms.zero_forces();
    allatoms.zero_velocities();

    rmg::printlog ("\n ==============================================\n\n");

    /*Reset timers, so that they do not count preceeding quench run if any */
    /*reset_timers (); */

    /*Also reset number of scf steps */
    ct.total_scf_steps = 0;

    int FP0_BASIS = Rmg_G->get_P0_BASIS(Rmg_G->default_FG_RATIO);
    fgobj<double> arho, trho;

    double_2d_array crds;
    crds.resize(boost::extents[Atoms.size()][3]);
    /* enforce periodic boundary conditions on the ions */
    allatoms.enforce_pbc();

    /* Save coordinates r(t) */
    for (size_t ion = 0, i_end = Atoms.size(); ion < i_end; ++ion)
    {
        ION &Atom = Atoms[ion];
        crds[ion][0] = Atom.crds[0];
        crds[ion][1] = Atom.crds[1];
        crds[ion][2] = Atom.crds[2];
    }

    int beeman_order = 4;
    double r_pred_coeff[4]{323.0/360.0, -264.0/360.0, 159.0/360.0, -38.0/360.0 };
    double r_corr_coeff[4]{38.0/360.0, 171.0/360.0, -36.0/360.0, 7.0/360.0 };
    double v_corr_coeff[4]{135.0/360.0, 285.0/360, -75.0/360.0, 15.0/360.0};
//    double r_pred_coeff[4]{342.0/480.0, -114.0/480.0, 12.0/480.0, 0.0 };
//    double r_corr_coeff[4]{60.0/360.0, 150.0/360.0, -30.0/360.0, 0.0 };
//    double v_corr_coeff[4]{9.0/24.0, 19.0/24.0, -5.0/24.0, 1.0/24.0};
//
//    3-rd order
//    double r_pred_coeff[4]{4.0/6.0, -1.0/6.0, 0.0, 0.0 };
//    double r_corr_coeff[4]{1.0/6.0, 2.0/6.0, 0.0,0.0};
//    double v_corr_coeff[4]{1.0/2.0, 1.0/2.0, 0.0, 0.0};
    /* begin the md loop */
    for (ct.md_steps = 1; ct.md_steps < ct.max_md_steps; ct.md_steps++)
    {

        
        //// --- STEP 1: UPGRADED PREDICTOR STAGE ---
        /* r_predic the positions a full timestep */
        // 6th-order Position Predictor expansion
        ////    r_pred[i] = p.r + p.v * dt + (dt * dt / 360.0) * (342.0 * p.a - 114.0 * p.a_prev1 + 12.0 * p.a_prev2);
        ////    Atoms.crds = r_pred 
        ////    p.r  = crds 
        for (size_t ion = 0, i_end = Atoms.size(); ion < i_end; ++ion)
        {
            ION &Atom = Atoms[ion];
            if (Atom.movable[0] || Atom.movable[1] || Atom.movable[2])
            {
                double mass = Species[Atom.species].atomic_mass * mu_me;
                //double t1 = ct.iondt *ct.iondt/360.0/ mass;
                double t1 = ct.iondt *ct.iondt/mass;
                for(int ic = 0; ic < 3; ic++)
                {
                    Atom.crds[ic] = crds[ion][ic] + ct.iondt * Atom.velocity[ic];
                    for(int iq = 0; iq < beeman_order; iq++)
                    {
                        Atom.crds[ic] += t1 * r_pred_coeff[iq] * Atom.force[iq][ic];
                    }
                }

                /* move atoms back to supercell */
                to_crystal (Atom.xtal, Atom.crds);
                to_cartesian (Atom.xtal, Atom.crds);

            }
        }

        // --- STEP 2: FORCE EVALUATION STAGE at position r_pred --
        /* Update items that change when the ionic coordinates change */
        RmgTimer *RT1 = new RmgTimer("1-TOTAL: run: ReinitIonicPotentials");
        ReinitIonicPotentials (Kptr, vnuc.data(), rhocore.data(), rhoc.data());
        delete RT1;

        //ct.mix = 0.0;
        spinobj<double> rho_save;
        rho_save = *Kptr[0]->rho;
        //ct.forceflag = NSCF;
        //ct.max_scf_steps = 10;
        Quench (Kptr, false);

        //ct.forceflag = TDDFT_CVE;
        *Kptr[0]->rho = rho_save;

        // Reset mixing
        MixRho(NULL, NULL, NULL, NULL, NULL, NULL, Kptr[0]->ControlMap, true);

        // tddft at position r_pred
        if(ct.tddft_mode == VECTOR_POT)
        {
            tddftobj_vp->tddft_md();
        }
        else
        {
            tddftobj->tddft_md();
        }
        write_force();

        /* zero out the non-moving atoms */
        allatoms.zero_forces();
        allatoms.zero_velocities();

        // --- STEP 3: UPGRADED CORRECTOR & HISTORY UPDATE STAGE ---

         // Multi-step implicit refinement utilizing the fresh a_pred
         //   p.r = p.r + p.v * dt + (dt * dt / 360.0) * (60.0 * a_pred[i] + 150.0 * p.a - 30.0 * p.a_prev1);
         //   p.v = p.v + (dt / 24.0) * (9.0 * a_pred[i] + 19.0 * p.a - 5.0 * p.a_prev1 + p.a_prev2);
        for (size_t ion = 0, i_end = Atoms.size(); ion < i_end; ++ion)
        {
            ION &Atom = Atoms[ion];
            if (Atom.movable[0] || Atom.movable[1] || Atom.movable[2])
            {
                double mass = Species[Atom.species].atomic_mass * mu_me;
                double t1 = ct.iondt *ct.iondt/ mass;
                double t2 = ct.iondt/mass;
                for(int ic = 0; ic < 3; ic++)
                {
                    crds[ion][ic] = crds[ion][ic] + ct.iondt * Atom.velocity[ic];
                    for(int iq = 0; iq < beeman_order; iq++)
                    {
                        crds[ion][ic] += t1 * r_corr_coeff[iq] * Atom.force[iq][ic];
                        Atom.velocity[ic] += t2 * v_corr_coeff[iq] * Atom.force[iq][ic];
                    }
                }
            }
        }


        /* get the total T of the system */
        ct.ionke = allatoms.kinetic_energy();
        iontemp = ct.ionke * 2.0 / (3.0 * (double) N * kB);

        // get center of mass velocities
        double vx, vy, vz;
        center_of_mass_velocity(vx, vy, vz);

        rmg::printlog ("\n @CVE %5d  %15.10f  %15.10f  %15.10f  %15.10f  %10.8e",
                ct.md_steps, ct.TOTAL, ct.ionke, ct.TOTAL + ct.ionke, iontemp, trms);
        rmg::printlog("\n Center of mass velocity (%15.10f, %15.10f, %15.10f)", vx, vy, vz);



        rmg::printlog ("\n Total number of SCF steps so far %d", ct.total_scf_steps);


        fflush (NULL);
    }                           /* end for ct.md_steps */

    WriteRestart (ct.outfile, vh.data(), rho.up.data(), rho.dw.data(), vxc.data(), Kptr);


    if (pct.gridpe == 0)
        rmg::printlog ("\n Total number of SCF steps %d", ct.total_scf_steps);


}

