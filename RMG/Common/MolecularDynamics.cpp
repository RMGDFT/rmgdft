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
#include "RmgSumAll.h"
#include "Atomic.h"
#include "transition.h"


/* Local function prototypes */
void init_nose (void);
void nose_velup1 (void);
void nose_posup (void);
void nose_velup2 (void);
void nose_energy (double *, double *);
void center_of_mass_velocity(double &vx, double &vy, double &vz);

void velup1 (void);
//void posup (void);
void velup2 (void);

void rms_disp (double *, double *);


template void MolecularDynamics<double>(Kpoint<double> **, double *, double *, double *, double *, double *, double *, double *);
template void MolecularDynamics<std::complex<double> >(Kpoint<std::complex<double>> **, double *, double *, double *, double *, double *, double *, double *);


template <typename KpointType>
void MolecularDynamics (Kpoint<KpointType> **Kptr, double * vxc, double * vh, double * vnuc,
             double * rho, double * rho_oppo, double * rhoc, double * rhocore)
{
    double rms[3], trms;
    double nosekin, nosepot;
    double iontemp;
    int N;
    std::vector<double> RMSdV;
    static double *rhodiff;

    // Prime the pump for extrapolations if this is an initial run.
    if(ct.runflag != RESTART)
    {
        Quench (vxc, vh, vnuc, rho, rho_oppo, rhocore, rhoc, Kptr, false);
        WriteRestart (ct.outfile, vh, rho, rho_oppo, vxc, Kptr);
    }


    /*Get some memory */
    double *crdsx = new double[Atoms.size()];
    double *crdsy = new double[Atoms.size()];
    double *crdsz = new double[Atoms.size()];



    if (pct.gridpe == 0)
    {
        rmg_printf ("\n ==============================================");

        /* print out the title */
        switch (ct.forceflag)
        {
            case MD_CVE:
                rmg_printf ("\n Constant Volume & Energy Molecular Dynamics");
                break;
            case MD_CVT:
                if (ct.tcontrol == T_NOSE_CHAIN)
                {
                    rmg_printf ("\n Finite Temperature MD with Nose-Hoover Chains");
                }
                else
                {
                    rmg_printf ("\n Finite Temperature MD with Anderson Rescaling");
                }                   /* end of if */
                break;
            case MD_CPT:
                rmg_printf ("\n Constant Pressure and Temperature Molecular Dynamics");
                break;
        }                       /* end of switch */

        switch (ct.mdorder)
        {
            case ORDER_2:
                rmg_printf ("\n Integration with Velocity Verlet");
                break;
            case ORDER_3:
                rmg_printf ("\n Integration with 3rd Order Beeman Velocity Verlet");
                break;
            case ORDER_5:
                rmg_printf ("\n Integration with 5th Order Beeman Velocity Verlet");
                break;
        }                       /* end of switch */
    }                           /* end of if pe */


    ct.ionke = 0.0;

    /* define Boltzmann param */
    double kB = 1.0 / (11605.0 * Ha_eV);

    /* count up moving atoms */
    N = 0;
    for (size_t ion = 0, i_end = Atoms.size(); ion < i_end; ++ion)
    {
        if (Atoms[ion].movable) N++;
    }
    ct.nose.N = N;

    /* check to see if random velocities are needed */
    if (ct.nose.randomvel && (ct.runflag != RESTART))
    {
        if (pct.gridpe == 0)
        {
            rmg_printf ("\n\n Initializing temperature to %14.10f K\n", ct.nose.temp);
        }
        ranv ();
    }


    /* init nose variables */
    if (ct.forceflag == MD_CVT && ct.tcontrol == T_NOSE_CHAIN)
        init_nose ();

    /* zero out the non-moving atoms */
    for (size_t ion = 0, i_end = Atoms.size(); ion < i_end; ++ion)
    {
        if (Atoms[ion].movable == 0)
        {
            ION &Atom = Atoms[ion];
            Atom.ZeroForces();
            Atom.ZeroVelocity();
        }
    }

    if (pct.gridpe == 0)
    {
        rmg_printf ("\n ==============================================\n\n");
    }

    /*Reset timers, so that they do not count preceeding quench run if any */
    /*reset_timers (); */

    /*Also reset number of scf steps */
    ct.total_scf_steps = 0;


    /* begin the md loop */
    for (ct.md_steps = 1; ct.md_steps < ct.max_md_steps; ct.md_steps++)
    {

        /* enforce periodic boundary conditions on the ions */
        for (size_t ion = 0, i_end = Atoms.size(); ion < i_end; ++ion)
        {
            ION &Atom = Atoms[ion];

            /* to_crystal enforces periodic boundary conditions */
            to_crystal (Atom.xtal, Atom.crds);
            to_cartesian (Atom.xtal, Atom.crds);
        }

        /* Save coordinates */
        for (size_t ion = 0, i_end = Atoms.size(); ion < i_end; ++ion)
        {
            ION &Atom = Atoms[ion];
            crdsx[ion] = Atom.crds[0];
            crdsy[ion] = Atom.crds[1];
            crdsz[ion] = Atom.crds[2];
        }

        /* Do a halfstep update of the velocities */
        velup1 ();

        // Get atomic rho for this ionic configuration and subtract from current rho
        int FP0_BASIS = Rmg_G->get_P0_BASIS(Rmg_G->default_FG_RATIO);
        double *arho = new double[FP0_BASIS];
        LcaoGetAtomicRho(arho);

        for(int idx = 0;idx < FP0_BASIS;idx++) rho[idx] -= arho[idx];

        if(rhodiff == NULL)
        {
            rhodiff = new double[FP0_BASIS];
            for(int idx = 0;idx < FP0_BASIS;idx++) rhodiff[idx] = rho[idx];
        }
        else
        {
            double *trho = new double[FP0_BASIS];
            for(int idx = 0;idx < FP0_BASIS;idx++) trho[idx] = rho[idx];
            for(int idx = 0;idx < FP0_BASIS;idx++) rho[idx] = 2.0*rho[idx] - rhodiff[idx];
            for(int idx = 0;idx < FP0_BASIS;idx++) rhodiff[idx] = trho[idx];
            delete [] trho;
        }


        /* Update the positions a full timestep */
        //posup ();
	move_ions (ct.iondt);

        // Get atomic rho for new configuration and add back to rho
        LcaoGetAtomicRho(arho);
        for(int idx = 0;idx < FP0_BASIS;idx++) rho[idx] += arho[idx];
        delete [] arho;

        /* Update items that change when the ionic coordinates change */
        RmgTimer *RT1 = new RmgTimer("1-TOTAL: run: ReinitIonicPotentials");
        ReinitIonicPotentials (Kptr, vnuc, rhocore, rhoc);
        delete RT1;

        // Reset mixing
        MixRho(NULL, NULL, NULL, NULL, NULL, NULL, Kptr[0]->ControlMap, true);

        /* converge to the ground state at the final positions */
        Quench (vxc, vh, vnuc, rho, rho_oppo, rhocore, rhoc, Kptr, true);


        /* zero out the non-moving atoms */
        for (size_t ion = 0, i_end = Atoms.size(); ion < i_end; ++ion)
        {
            if (Atoms[ion].movable == 0)
            {
                ION &Atom = Atoms[ion];
                Atom.ZeroForces();
                Atom.ZeroVelocity();
            }
        }

        /* Do another halfstep update of the velocities */
        /* to update them to the next time step */
        velup2 ();

        if (pct.gridpe == 0) rmg_printf ("\n Writing data to output file ...\n");
        WriteRestart (ct.outfile, vh, rho, rho_oppo, vxc, Kptr);

        // Extrapolate orbitals after first step
        ExtrapolateOrbitals(ct.outfile, Kptr);

        /* calculate the nose thermostat energies */
        if (ct.forceflag == MD_CVT && ct.tcontrol == T_NOSE_CHAIN)
        {
            nose_energy (&nosekin, &nosepot);
        }
        else
        {
            nosekin = 0.0;
            nosepot = 0.0;
        }

        /* calculate rms displacement */
        rms_disp (&rms[0], &trms);

        /* get the total T of the system */
        iontemp = ct.ionke * 2.0 / (3.0 * (double) N * kB);

        // get center of mass velocities
        double vx, vy, vz;
        center_of_mass_velocity(vx, vy, vz);


        if (pct.gridpe == 0)
        {
            switch (ct.forceflag)
            {

            case MD_CVE:
                rmg_printf ("\n @CVE %5d  %15.10f  %15.10f  %15.10f  %15.10f  %10.8e",
                        ct.md_steps, ct.TOTAL, ct.ionke, ct.TOTAL + ct.ionke, iontemp, trms);
                rmg_printf("\n Center of mass velocity (%15.10f, %15.10f, %15.10f)", vx, vy, vz);
                break;
            case MD_CVT:
                if (ct.tcontrol == T_NOSE_CHAIN)
                {
                    rmg_printf ("\n @CVT-NOSE %5d  %15.10f  %15.10f  %15.10f  %15.10f  %15.10f  %10.8e",
                            ct.md_steps, ct.TOTAL, ct.ionke, nosekin + nosepot,
                            ct.TOTAL + ct.ionke + nosekin + nosepot, iontemp, trms);
                }
                if (ct.tcontrol == T_AND_SCALE)
                {
                    rmg_printf ("\n @CVT-ANDERSON %5d  %15.10f  %15.10f  %15.10f  %15.10f  %10.8e",
                            ct.md_steps, ct.TOTAL, ct.ionke, ct.TOTAL + ct.ionke, iontemp, trms);
                }
                break;
            case MD_CPT:
                rmg_error_handler (__FILE__, __LINE__,"MD_CPT is not programmed yet.");
                break;

            default:
                throw RmgFatalException() << "Undefined output type" << " in " << __FILE__ << " at line " << __LINE__ << "\n";

            }                   /* end of switch */


            rmg_printf ("\n Total number of SCF steps so far %d", ct.total_scf_steps);

        }                       /* end of if */

        fflush (NULL);
    }                           /* end for ct.md_steps */


    // Final quench at these ionic positions without computing forces so the
    // final wavefunctions are converged. Write the restart file then compute
    Quench (vxc, vh, vnuc, rho, rho_oppo, rhocore, rhoc, Kptr, false);
    WriteRestart (ct.outfile, vh, rho, rho_oppo, vxc, Kptr);


    if (pct.gridpe == 0)
        rmg_printf ("\n Total number of SCF steps %d", ct.total_scf_steps);



    delete [] crdsx;
    delete [] crdsy;
    delete [] crdsz;

}


void init_nose ()
{
    int N, jc;
    double wNose, tau_nose, kB, step;
    double inittemp, nosesteps;

    step = ct.iondt;

    /* define Boltzmann param and number of moving atoms */
    kB = 1.0 / (11605.0 * Ha_eV);
    N = ct.nose.N;

    /* initialize nose parameters */
    wNose = ct.nose.fNose * 2.0 * PI * 2.418e-5;
    tau_nose = 1.0 / wNose;
    ct.nose.k0 = 1.5 * N * kB * ct.nose.temp;

    /* init thermostat masses */
    ct.nose.xq[0] = 2.0 * ct.nose.k0 * (tau_nose * tau_nose);
    for (jc = 1; jc < ct.nose.m; jc++)
    {

        ct.nose.xq[jc] = ct.nose.temp * kB * (tau_nose * tau_nose);

    }

    /* calculate ion KE */
    ct.ionke = 0.0;
    for (size_t ion = 0, i_end = Atoms.size(); ion < i_end; ++ion)
    {
        ct.ionke += Atoms[ion].GetKineticEnergy();
    }

    inittemp = ct.ionke * 2.0 / (3.0 * (double) N * kB);

    /* init thermostat forces */
    ct.nose.xf[ct.fpt[0]][0] = 2.0 * (ct.ionke - ct.nose.k0) / ct.nose.xq[0];
    for (jc = 1; jc < ct.nose.m; jc++)
    {

        ct.nose.xf[ct.fpt[0]][jc] = (ct.nose.xq[jc - 1] * ct.nose.xv[jc - 1] * ct.nose.xv[jc - 1]
                                     - ct.nose.temp * kB) / ct.nose.xq[jc];

    }

    nosesteps = 2.0 * PI / (wNose * step);

    if (pct.gridpe == 0)
    {
        rmg_printf ("\n Nose Frequency (THz) = %14.10f ", ct.nose.fNose);
        rmg_printf ("\n Steps/Oscillation    = %14.10f ", nosesteps);
        rmg_printf ("\n Target Temp          = %14.10f ", ct.nose.temp);
        rmg_printf ("\n Initial Temp         = %14.10f ", inittemp);
        rmg_printf ("\n i        xx        xv         xq     xf ");
        for (jc = 0; jc < ct.nose.m; jc++)
        {
            rmg_printf ("\n %d %14.10f %14.10f %16.10f %14.10f", jc,
                    ct.nose.xx[jc], ct.nose.xv[jc], ct.nose.xq[jc], ct.nose.xf[ct.fpt[0]][jc]);
        }
        rmg_printf ("\n ==============================================");
    }

}                               /* end of init_nose */


void velup1 ()
{
    int ic;
    double step, mass, kB;
    double t1, t2;
    double scale = 1.0;
    double temperature;

    step = ct.iondt;

    /* define Boltzmann param */
    kB = 1.0 / (11605.0 * Ha_eV);

    switch (ct.forceflag)
    {

    case MD_CVE:
        scale = 1.0;
        break;
    case MD_CVT:
        if (ct.tcontrol == T_NOSE_CHAIN)
        {
            nose_velup1 ();
            scale = exp (-0.5 * step * ct.nose.xv[0]);
        }
        if (ct.tcontrol == T_AND_SCALE)
        {

            ct.ionke = 0.0;
            /* Loop over ions */
            for (size_t ion = 0, i_end = Atoms.size(); ion < i_end; ++ion)
            {
                ct.ionke += Atoms[ion].GetKineticEnergy();
            }

            temperature = ct.ionke * 2.0 / (3.0 * ct.nose.N * kB);

            if (ct.ionke > 0.0)
            {
                scale = sqrt (ct.nose.temp / temperature);
                if (pct.gridpe == 0)
                    rmg_printf ("\ntscale=%f\n", scale);
            }
            else
            {
                scale = 1.0;
            }

        }
        break;
    }

    /* Loop over ions */
    ct.ionke = 0.0;
    for (size_t ion = 0, i_end = Atoms.size(); ion < i_end; ++ion)
    {

        ION &Atom = Atoms[ion];

        if (Atom.movable)
        {

            /* Get ionic mass */
            mass = Species[Atom.species].atomic_mass * mu_me;

            /* update velocity by one-half timestep */
            switch (ct.mdorder)
            {

                case ORDER_2:
                    t1 = step / mass;
                    t2 = t1 / 2.0;

                    /* step velocities forward */
                    for (ic = 0; ic < 3; ic++)
                    {

                        Atom.velocity[ic] *= scale;
                        Atom.velocity[ic] += t2 * Atom.force[0][ic];

                    }
                    break;

                case ORDER_3:
                    t1 = step / mass;
                    t2 = t1 / 6.0;

                    /* step velocities forward */
                    for (ic = 0; ic < 3; ic++)
                    {

                        Atom.velocity[ic] *= scale;
                        Atom.velocity[ic] += t2 * (4.0 * Atom.force[ct.fpt[0]][ic]
                                                    - 1.0 * Atom.force[ct.fpt[1]][ic]);

                    }
                    break;

                case ORDER_5:
                    t1 = step / mass;
                    t2 = t1 / 360.0;

                    /* step velocities forward */
                    for (ic = 0; ic < 3; ic++)
                    {

                        Atom.velocity[ic] *= scale;
                        Atom.velocity[ic] += t2 * (323.0 * Atom.force[ct.fpt[0]][ic]
                                                    - 264.0 * Atom.force[ct.fpt[1]][ic]
                                                    + 159.0 * Atom.force[ct.fpt[2]][ic]
                                                    - 38.0 * Atom.force[ct.fpt[3]][ic]);

                    }
                    break;

            }                   /* end of switch */

            ct.ionke += Atoms[ion].GetKineticEnergy();

        }                       /* if */

    }                           /* end for */

}                               /* end of velup1 */


#if 0
/* Function move_ions takes care of mving ions, it also keeps ionic history and other stuff*/
void posup ()
{
    int ion;
    ION *iptr;
    double step;

    step = ct.iondt;

    /* update nose thermostats */
    if (ct.forceflag == MD_CVT && ct.tcontrol == T_NOSE_CHAIN)
        nose_posup ();

    /* Loop over ions */
    for (size_t ion = 0, i_end = Atoms.size(); ion < i_end; ++ion)
    {

        /* Get ion pointer */
        iptr = &Atoms[ion];

        /* update positions by one full timestep */
        if (iptr->movable)
        {

            iptr->crds[0] += step * iptr->velocity[0];
            iptr->crds[1] += step * iptr->velocity[1];
            iptr->crds[2] += step * iptr->velocity[2];

            /* move atoms back to supercell */
            to_crystal (iptr->xtal, iptr->crds);
            to_cartesian (iptr->xtal, iptr->crds);

        }


    }

}                               /* end of posup */
#endif


void velup2 ()
{
    int ic;
    double step, mass;
    double t1, t2;
    double scale = 1.0;

    step = ct.iondt;

    switch (ct.forceflag)
    {

    case MD_CVE:
        scale = 1.0;
        break;
    case MD_CVT:
        if (ct.tcontrol == T_NOSE_CHAIN)
        {
            scale = exp (-0.5 * step * ct.nose.xv[0]);
        }
        if (ct.tcontrol == T_AND_SCALE)
        {
            scale = 1.0;
        }
        break;

    }                           /* end of switch */

    ct.ionke = 0.0;

    /* Loop over ions */
    for (size_t ion = 0, i_end = Atoms.size(); ion < i_end; ++ion)
    {

        ION &Atom = Atoms[ion];

        /* Get ionic mass */
        mass = Species[Atom.species].atomic_mass * mu_me;

        if (Atom.movable)
        {

            /* update velocity by one-half step */
            switch (ct.mdorder)
            {

            case ORDER_2:
                t1 = step / mass;
                t2 = t1 / 2.0;

                for (ic = 0; ic < 3; ic++)
                {
                    Atom.velocity[ic] += t2 * Atom.force[0][ic];
                    Atom.velocity[ic] *= scale;
                }
                break;

            case ORDER_3:
                t1 = step / mass;
                t2 = t1 / 6.0;

                for (ic = 0; ic < 3; ic++)
                {
                    Atom.velocity[ic] += t2 * (2.0 * Atom.force[ct.fpt[0]][ic]
                                                + 1.0 * Atom.force[ct.fpt[1]][ic]);
                    Atom.velocity[ic] *= scale;
                }
                break;

            case ORDER_5:
                t1 = step / mass;
                t2 = t1 / 360.0;

                for (ic = 0; ic < 3; ic++)
                {
                    Atom.velocity[ic] += t2 * (97.0 * Atom.force[ct.fpt[0]][ic]
                                                + 114.0 * Atom.force[ct.fpt[1]][ic]
                                                - 39.0 * Atom.force[ct.fpt[2]][ic]
                                                + 8.0 * Atom.force[ct.fpt[3]][ic]);
                    Atom.velocity[ic] *= scale;
                }
                break;

            }                   /* end of switch */

        }                       /* end of if */

        ct.ionke += Atom.GetKineticEnergy();

    }                           /* end for */

    /* perform second half update of nose coords */
    if (ct.forceflag == MD_CVT && ct.tcontrol == T_NOSE_CHAIN)
    {
        nose_velup2 ();
    }

}                               /* end of velup2 */





void rms_disp (double * rms, double * trms)
{

    double t1, t2, t3;

    *trms = 0.0;
    rms[0] = rms[1] = rms[2] = 0.0;
    for (size_t ion = 0, i_end = Atoms.size(); ion < i_end; ++ion)
    {

        ION &Atom = Atoms[ion];

        t1 = Atom.xtal[0] - Atom.ixtal[0];
        if (t1 < 0.0)
            t1 *= -1.0;
        if (t1 > 0.5)
            t1 = 1.0 - t1;

        t2 = Atom.xtal[1] - Atom.ixtal[1];
        if (t2 < 0.0)
            t2 *= -1.0;
        if (t2 > 0.5)
            t2 = 1.0 - t2;

        t3 = Atom.xtal[2] - Atom.ixtal[2];
        if (t3 < 0.0)
            t3 *= -1.0;
        if (t3 > 0.5)
            t3 = 1.0 - t3;

        rms[0] = t1;
        rms[1] = t2;
        rms[2] = t3;

        *trms += metric (rms) / ((double) ct.nose.N);

    }

}                               /* end of rms_disp */

void nose_velup1 ()
{
    int jc;
    double step;
    double scale;
    int mn;

    step = ct.iondt;

    /* do nose thermostat velocity update */
    mn = ct.nose.m - 1;
    switch (ct.mdorder)
    {

    case ORDER_2:
        ct.nose.xv[mn] += 0.5 * step * ct.nose.xf[0][mn];
        for (jc = mn - 1; jc >= 0; jc--)
        {
            scale = exp (-0.5 * step * ct.nose.xv[jc + 1]);
            ct.nose.xv[jc] = scale * ct.nose.xv[jc] + 0.5 * step * ct.nose.xf[0][jc];
        }
        break;

    case ORDER_3:
        ct.nose.xv[mn] +=
            step / 6.0 * (4.0 * ct.nose.xf[ct.fpt[0]][mn] - ct.nose.xf[ct.fpt[1]][mn]);
        for (jc = mn - 1; jc >= 0; jc--)
        {
            scale = exp (-0.5 * step * ct.nose.xv[jc + 1]);
            ct.nose.xv[jc] = scale * ct.nose.xv[jc] +
                step / 6.0 * (4.0 * ct.nose.xf[ct.fpt[0]][jc] - ct.nose.xf[ct.fpt[1]][jc]);
        }
        break;

    case ORDER_5:
        ct.nose.xv[mn] += step / 360.0 * (323.0 * ct.nose.xf[ct.fpt[0]][mn]
                                          - 264.0 * ct.nose.xf[ct.fpt[1]][mn]
                                          + 159.0 * ct.nose.xf[ct.fpt[2]][mn]
                                          - 38.0 * ct.nose.xf[ct.fpt[3]][mn]);
        for (jc = mn - 1; jc >= 0; jc--)
        {
            scale = exp (-0.5 * step * ct.nose.xv[jc + 1]);
            ct.nose.xv[jc] = scale * ct.nose.xv[jc] +
                step / 360.0 * (323.0 * ct.nose.xf[ct.fpt[0]][jc]
                                - 264.0 * ct.nose.xf[ct.fpt[1]][jc]
                                + 159.0 * ct.nose.xf[ct.fpt[2]][jc]
                                - 38.0 * ct.nose.xf[ct.fpt[3]][jc]);
        }
        break;
    }

}                               /* end of nose_velup1 */


void nose_posup ()
{

    int jc;
    double step;

    step = ct.iondt;

    /* loop over thermostat positions */
    for (jc = 0; jc < ct.nose.m; jc++)
    {

        ct.nose.xx[jc] += step * ct.nose.xv[jc];

    }

}                               /* end of nose_posup */

void nose_velup2 ()
{
    int jc;
    double kB, step;
    double scale;
    int mn;

    step = ct.iondt;

    /* define Boltzmann param */
    kB = 1.0 / (11605.0 * Ha_eV);

    /* calc new nose velocities and forces */
    if (ct.nose.m == 1)
    {
        ct.nose.xf[ct.fpt[0]][0] = 2.0 * (ct.ionke - ct.nose.k0) / ct.nose.xq[0];
        switch (ct.mdorder)
        {
        case ORDER_2:
            ct.nose.xv[0] += step / 2.0 * ct.nose.xf[0][0];
            break;
        case ORDER_3:
            ct.nose.xv[0] +=
                step / 6.0 * (2.0 * ct.nose.xf[ct.fpt[0]][0] + ct.nose.xf[ct.fpt[1]][0]);
            break;
        case ORDER_5:
            ct.nose.xv[0] += step / 360.0 * (97.0 * ct.nose.xf[ct.fpt[0]][0]
                                             + 114.0 * ct.nose.xf[ct.fpt[1]][0]
                                             - 39.0 * ct.nose.xf[ct.fpt[2]][0]
                                             + 8.0 * ct.nose.xf[ct.fpt[3]][0]);
            break;
        }                       /* end of switch */
    }
    else
    {
        scale = exp (-0.5 * step * ct.nose.xv[1]);
        ct.nose.xf[ct.fpt[0]][0] = 2.0 * (ct.ionke - ct.nose.k0) / ct.nose.xq[0];

        switch (ct.mdorder)
        {
        case ORDER_2:
            ct.nose.xv[0] = scale * (ct.nose.xv[0] + step / 2.0 * ct.nose.xf[0][0]);
            break;
        case ORDER_3:
            ct.nose.xv[0] =
                scale * (ct.nose.xv[0] +
                         step / 6.0 * (2.0 * ct.nose.xf[ct.fpt[0]][0] + ct.nose.xf[ct.fpt[1]][0]));
            break;
        case ORDER_5:
            ct.nose.xv[0] = scale * (ct.nose.xv[0] + step / 360.0 * (97.0 * ct.nose.xf[ct.fpt[0]][0]
                                                                     +
                                                                     114.0 *
                                                                     ct.nose.xf[ct.fpt[1]][0] -
                                                                     39.0 *
                                                                     ct.nose.xf[ct.fpt[2]][0] +
                                                                     8.0 *
                                                                     ct.nose.xf[ct.fpt[3]][0]));
            break;
        }                       /* end of switch */

        switch (ct.mdorder)
        {
        case ORDER_2:
            for (jc = 1; jc < ct.nose.m - 1; jc++)
            {
                scale = exp (-0.5 * step * ct.nose.xv[jc + 1]);
                ct.nose.xf[ct.fpt[0]][jc] = (ct.nose.xq[jc - 1] * ct.nose.xv[jc - 1] *
                                             ct.nose.xv[jc - 1] - ct.nose.temp * kB) /
                    ct.nose.xq[jc];
                ct.nose.xv[jc] = scale * (ct.nose.xv[jc] + step / 2.0 * ct.nose.xf[ct.fpt[0]][jc]);
            }
            mn = ct.nose.m - 1;
            ct.nose.xf[ct.fpt[0]][mn] =
                (ct.nose.xq[mn - 1] * ct.nose.xv[mn - 1] * ct.nose.xv[mn - 1] -
                 ct.nose.temp * kB) / ct.nose.xq[mn];
            ct.nose.xv[mn] = ct.nose.xv[mn] + step / 2.0 * ct.nose.xf[ct.fpt[0]][mn];
            break;
        case ORDER_3:
            for (jc = 1; jc < ct.nose.m - 1; jc++)
            {
                scale = exp (-0.5 * step * ct.nose.xv[jc + 1]);
                ct.nose.xf[ct.fpt[0]][jc] = (ct.nose.xq[jc - 1] * ct.nose.xv[jc - 1] *
                                             ct.nose.xv[jc - 1] - ct.nose.temp * kB) /
                    ct.nose.xq[jc];
                ct.nose.xv[jc] =
                    scale * (ct.nose.xv[jc] +
                             step / 6.0 * (2.0 * ct.nose.xf[ct.fpt[0]][jc] +
                                           ct.nose.xf[ct.fpt[1]][jc]));
            }
            mn = ct.nose.m - 1;
            ct.nose.xf[ct.fpt[0]][mn] =
                (ct.nose.xq[mn - 1] * ct.nose.xv[mn - 1] * ct.nose.xv[mn - 1] -
                 ct.nose.temp * kB) / ct.nose.xq[mn];
            ct.nose.xv[mn] =
                ct.nose.xv[mn] + step / 6.0 * (2.0 * ct.nose.xf[ct.fpt[0]][mn] +
                                               ct.nose.xf[ct.fpt[1]][mn]);
            break;
        case ORDER_5:
            for (jc = 1; jc < ct.nose.m - 1; jc++)
            {
                scale = exp (-0.5 * step * ct.nose.xv[jc + 1]);
                ct.nose.xf[ct.fpt[0]][jc] = (ct.nose.xq[jc - 1] * ct.nose.xv[jc - 1] *
                                             ct.nose.xv[jc - 1] - ct.nose.temp * kB) /
                    ct.nose.xq[jc];
                ct.nose.xv[jc] =
                    scale * (ct.nose.xv[jc] +
                             step / 360.0 * (97.0 * ct.nose.xf[ct.fpt[0]][jc] +
                                             114.0 * ct.nose.xf[ct.fpt[1]][jc] -
                                             39.0 * ct.nose.xf[ct.fpt[2]][jc] +
                                             8.0 * ct.nose.xf[ct.fpt[3]][jc]));
            }
            mn = ct.nose.m - 1;
            ct.nose.xf[ct.fpt[0]][mn] =
                (ct.nose.xq[mn - 1] * ct.nose.xv[mn - 1] * ct.nose.xv[mn - 1] -
                 ct.nose.temp * kB) / ct.nose.xq[mn];
            ct.nose.xv[mn] =
                ct.nose.xv[mn] + step / 360.0 * (97.0 * ct.nose.xf[ct.fpt[0]][mn] +
                                                 114.0 * ct.nose.xf[ct.fpt[1]][mn] -
                                                 39.0 * ct.nose.xf[ct.fpt[2]][mn] +
                                                 8.0 * ct.nose.xf[ct.fpt[3]][mn]);
            break;
        }                       /* end of switch */

    }                           /* end of if m */

}                               /* end of nose_velup2 */


void nose_energy (double * nosekin, double * nosepot)
{
    int jc;
    double kB;

    /* define Boltzmann param */
    kB = 1.0 / (11605.0 * Ha_eV);

    /* calculate total nose energy */
    *nosekin = 0.5 * ct.nose.xq[0] * ct.nose.xv[0] * ct.nose.xv[0];
    *nosepot = 2.0 * ct.nose.k0 * ct.nose.xx[0];
#if 0
    if (pct.gridpe == 0)
        rmg_printf ("\n @therm%d %14.10f %14.10f %14.10f %14.10f %14.10f %14.10f",
                0, ct.nose.xx[0], ct.nose.xv[0], ct.nose.xf[ct.fpt[0]][0],
                ct.nose.xf[ct.fpt[1]][0], ct.nose.xf[ct.fpt[2]][0], ct.nose.xf[ct.fpt[3]][0]);
#endif
    for (jc = 1; jc < ct.nose.m; jc++)
    {
        *nosekin += 0.5 * ct.nose.xq[jc] * ct.nose.xv[jc] * ct.nose.xv[jc];
        *nosepot += ct.nose.temp * kB * ct.nose.xx[jc];
#if 0
        if (pct.gridpe == 0)
            rmg_printf ("\n @therm%d %14.10f %14.10f %14.10f %14.10f %14.10f %14.10f",
                    jc, ct.nose.xx[jc], ct.nose.xv[jc], ct.nose.xf[ct.fpt[0]][jc],
                    ct.nose.xf[ct.fpt[1]][jc],
                    ct.nose.xf[ct.fpt[2]][jc], ct.nose.xf[ct.fpt[3]][jc]);
#endif

    }

#if 0
    if (pct.gridpe == 0)
        rmg_printf ("\n nose_energy: %20.10f %20.10f %20.10f", *nosekin, *nosepot, *nosekin + *nosepot);
#endif

}                               /* end of nose_energy */


void center_of_mass_velocity(double &vx, double &vy, double &vz)
{
    double px = 0.0, py = 0.0, pz = 0.0;
    double total_mass = 0.0;
    for (size_t ion = 0, i_end = Atoms.size(); ion < i_end; ++ion)
    {
        ION &Atom = Atoms[ion];
        double mass = Species[Atom.species].atomic_mass * mu_me;
        px += mass * Atom.velocity[0];
        py += mass * Atom.velocity[1];
        pz += mass * Atom.velocity[1];
        total_mass += mass;
    }
    
    vx = px / total_mass;
    vy = py / total_mass;
    vz = pz / total_mass;

}



void ranv (void)
{
    int i, ii, N;
    double c, p, tmass, vtot[3], ek, scale;
    double kB, mass;

    /* init the random number generator */
    srand (123);

    /* count up moving atoms */
    N = 0;
    for (size_t ion = 0, i_end = Atoms.size(); ion < i_end; ++ion)
    {

        if (Atoms[ion].movable)
            N++;

    }                           /* end for */

    /* define Boltzmann param */
    kB = 1.0 / (11605.0 * Ha_eV);

    if (ct.nose.temp > 0.0)
    {

        tmass = 0;
        vtot[0] = vtot[1] = vtot[2] = 0.0;

        /* randomize the velocity for each atom */
        for (size_t ion = 0, i_end = Atoms.size(); ion < i_end; ++ion)
        {

            ION &Atom = Atoms[ion];

            /* Get ionic mass */
            mass = Species[Atom.species].atomic_mass * mu_me;

            if (Atom.movable)
            {
                c = -kB * ct.nose.temp / mass;
                tmass += mass;
                for (i = 0; i < 3; i++)
                {

                    ii = rand ();
                    p = (double) ii / (double) RAND_MAX;

                    Atom.velocity[i] = sqrt (c * log (p));

                    if (p > 0.5) Atom.velocity[i] = -Atom.velocity[i];

                    vtot[i] += mass * Atom.velocity[i];

                }
            }                   /* end of if */

        }                       /* end for ion */

        /* scale set velocity by total mass */
        vtot[0] = vtot[0] / tmass;
        vtot[1] = vtot[1] / tmass;
        vtot[2] = vtot[2] / tmass;

        /* subtract out the center of mass velocity */
        for (size_t ion = 0, i_end = Atoms.size(); ion < i_end; ++ion)
        {

            ION &Atom = Atoms[ion];

            if (Atom.movable)
            {
                Atom.velocity[0] -= vtot[0];
                Atom.velocity[1] -= vtot[1];
                Atom.velocity[2] -= vtot[2];
            }
        }                       /* end of for ion */

        /* find kinetic energy */
        ek = 0.0;
        for (size_t ion = 0, i_end = Atoms.size(); ion < i_end; ++ion)
        {
            ek += Atoms[ion].GetKineticEnergy();
        }                       /* end of for ion */
        scale = sqrt (1.5 * N * kB * ct.nose.temp / ek);

        ek = 0.0;
        /* rescale to correct temperature */
        for (size_t ion = 0, i_end = Atoms.size(); ion < i_end; ++ion)
        {

            ION &Atom = Atoms[ion];

            /* Get ionic mass */
            mass = Species[Atom.species].atomic_mass * mu_me;

            if (Atom.movable)
            {
                Atom.velocity[0] *= scale;
                Atom.velocity[1] *= scale;
                Atom.velocity[2] *= scale;
                ek += Atom.velocity[0] * Atom.velocity[0] * mass;
                ek += Atom.velocity[1] * Atom.velocity[1] * mass;
                ek += Atom.velocity[2] * Atom.velocity[2] * mass;
            }

        }                       /* end of for ion */

        ct.ionke = 0.5 * ek;

    }                           /* end of if */

}                               /* end ranv */


