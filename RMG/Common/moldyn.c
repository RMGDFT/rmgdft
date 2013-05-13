/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

/****f* QMD-MGDFT/moldyn.c *****
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
 *   void moldyn(STATE *states, rmg_double_t *vxc, rmg_double_t *vh, rmg_double_t *vnuc,
 *               rmg_double_t *rho, rmg_double_t *rhoc, rmg_double_t *rhocore)
 *   Molecular Dynamics driver routine
 * INPUTS
 *   states: point to orbital structure (see main.h)
 *   vxc:    exchange-correlation potential
 *   vh:     Hartree potential
 *   vnuc:   pseudopotential
 *   rho:    total charge density
 *   rhoc:   compensating charge density
 *   rhocore: core charge density
 * OUTPUT
 *   coordinates are updated, and all above are also updated
 * PARENTS
 *   main.c
 * CHILDREN
 *   to_crystal.c to_cartesian.c get_nlop.c scf.c sortpsi.c subdiag.c get_te.c force.c
 * SOURCE
 */


#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "main.h"

/* Local function prototypes */
void init_nose (void);
void nose_velup1 (void);
void nose_posup (void);
void nose_velup2 (void);
void nose_energy (rmg_double_t *, rmg_double_t *);

void velup1 (void);
//void posup (void);
void velup2 (void);

void rms_disp (rmg_double_t *, rmg_double_t *);

int stepcount = 0;

void moldyn (STATE * states, rmg_double_t * vxc, rmg_double_t * vh, rmg_double_t * vnuc,
             rmg_double_t * rho, rmg_double_t * rho_oppo, rmg_double_t * rhoc, rmg_double_t * rhocore)
{

    rmg_double_t target;
    rmg_double_t *crdsx, *crdsy, *crdsz, rsteps;
    rmg_double_t kB, step;
    rmg_double_t tmix, tprjmix;
    rmg_double_t rms[3], trms;
    rmg_double_t nosekin, nosepot;
    rmg_double_t iontemp;
    int ion, it, steps1, isteps, nsteps = 1, nmsteps, N;
    int ic;
    ION *iptr;
    FILE *mfp = NULL, *dxmfp;
    FILE *xbsfp1 = NULL;
    char filename[60];
    char xbs_filename[60];

    bool CONVERGED = false;

    target = 0.0;


    /*Get some memory */
    my_malloc (crdsx, ct.num_ions, rmg_double_t);
    my_malloc (crdsy, ct.num_ions, rmg_double_t);
    my_malloc (crdsz, ct.num_ions, rmg_double_t);



    if (pct.gridpe == 0)
    {
        printf ("\n ==============================================");

        /* print out the title */
        switch (ct.forceflag)
        {
        case MD_CVE:
            printf ("\n Constant Volume & Energy Molecular Dynamics");
            break;
        case MD_CVT:
            if (ct.tcontrol == T_NOSE_CHAIN)
            {
                printf ("\n Finite Temperature MD with Nose-Hoover Chains");
            }
            else
            {
                printf ("\n Finite Temperature MD with Anderson Rescaling");
            }                   /* end of if */
            break;
        case MD_CPT:
            printf ("\n Constant Pressure and Temperature Molecular Dynamics");
            break;
        }                       /* end of switch */

        switch (ct.mdorder)
        {
        case ORDER_2:
            printf ("\n Integration with Velocity Verlet");
            break;
        case ORDER_3:
            printf ("\n Integration with 3rd Order Beeman Velocity Verlet");
            break;
        case ORDER_5:
            printf ("\n Integration with 5th Order Beeman Velocity Verlet");
            break;

        }                       /* end of switch */
    }                           /* end of if pe */




    /* number of substeps used to move charge density */
    isteps = 8;
    nmsteps = 4;

    ct.ionke = 0.0;
    step = ct.iondt;

    /* define Boltzmann param */
    kB = 1.0 / (11605.0 * Ha_eV);

    /* count up moving atoms */
    N = 0;
    for (ion = 0; ion < ct.num_ions; ion++)
    {

        if (ct.ions[ion].movable)
            N++;

    }
    ct.nose.N = N;

    /* check to see if random velocities are needed */
    if (ct.nose.randomvel)
    {
        if (pct.gridpe == 0)
        {
            printf ("\n\n Initializing temperature to %14.10f K\n", ct.nose.temp);
        }
        ranv ();
    }


    /* init nose variables */
    if (ct.forceflag == MD_CVT && ct.tcontrol == T_NOSE_CHAIN)
        init_nose ();

    /* zero out the non-moving atoms */
    for (ion = 0; ion < ct.num_ions; ion++)
    {

        if (ct.ions[ion].movable == 0)
        {

            iptr = &ct.ions[ion];

            for (ic = 0; ic < 4; ic++)
            {
                iptr->force[ic][0] = 0.0;
                iptr->force[ic][1] = 0.0;
                iptr->force[ic][2] = 0.0;
            }
            iptr->velocity[0] = 0.0;
            iptr->velocity[1] = 0.0;
            iptr->velocity[2] = 0.0;
        }

    }

    if (pct.gridpe == 0)
    {
        printf ("\n ==============================================");
    }

    /*Reset timers, so that they do not count preceeding quench run if any */
    /*reset_timers (); */

    /*Also reset number of scf steps */
    ct.total_scf_steps = 0;


    /* begin the md loop */
    for (ct.md_steps = 0; ct.md_steps < ct.max_md_steps; ct.md_steps++)
    {

        /* enforce periodic boundary conditions on the ions */
        for (it = 0; it < ct.num_ions; it++)
        {

            iptr = &ct.ions[it];

            /* to_crystal enforces periodic boundary conditions */
            to_crystal (iptr->xtal, iptr->crds);
            to_cartesian (iptr->xtal, iptr->crds);

        }

        /* Save coordinates */
        for (it = 0; it < ct.num_ions; it++)
        {

            iptr = &ct.ions[it];
            crdsx[it] = iptr->crds[0];
            crdsy[it] = iptr->crds[1];
            crdsz[it] = iptr->crds[2];

        }

        /* Do a halfstep update of the velocities */
        velup1 ();

        rsteps = 1.0 / (rmg_double_t) isteps;

        /* do nmsteps iterations to move the hamiltonian smoothly */
        /* to the next positions                            */
        for (steps1 = 0; steps1 < nmsteps; steps1++)
        {

            if (steps1 == 0)
                nsteps = 4;
            if (steps1 == 1)
                nsteps = 2;
            if (steps1 == 2)
                nsteps = 1;
            if (steps1 == 3)
                nsteps = 1;

            /* Step the ions forward */
            for (it = 0; it < ct.num_ions; it++)
            {

                iptr = &ct.ions[it];

                iptr->crds[0] += iptr->velocity[0] * step * rsteps * nsteps;
                iptr->crds[1] += iptr->velocity[1] * step * rsteps * nsteps;
                iptr->crds[2] += iptr->velocity[2] * step * rsteps * nsteps;

                to_crystal (iptr->xtal, iptr->crds);
                to_cartesian (iptr->xtal, iptr->crds);

            }

            /* Update items that change when the ionic coordinates change */
            reinit_ionic_pp (states, vnuc, rhocore, rhoc);

            /* Do an scf step */
            CONVERGED = scf (states, vxc, vh, vnuc, rho, rho_oppo, rhocore, rhoc);
            sortpsi (states);

        }                       /* end for */

        /* Reset coordinates to their original positions */
        for (it = 0; it < ct.num_ions; it++)
        {

            iptr = &ct.ions[it];

            iptr->crds[0] = crdsx[it];
            iptr->crds[1] = crdsy[it];
            iptr->crds[2] = crdsz[it];
            to_crystal (iptr->xtal, iptr->crds);
            to_cartesian (iptr->xtal, iptr->crds);

        }


        /* Update the positions a full timestep */
        //posup ();
	move_ions (ct.iondt);



        tmix = ct.mix;
        tprjmix = ct.prjmix;
        /* converge to the ground state at the final positions */



#if 0
        for (ct.scf_steps = 0, CONVERGED = false;
             ct.scf_steps < ct.max_scf_steps && !CONVERGED; ct.scf_steps++, ct.total_scf_steps++)
        {

            /* Perform a single self-consistent step */
            scf (states, vxc, vh, vnuc, rho, rhocore, rhoc);

            if (pct.gridpe == 0)
                printf ("\nThis is SCF step # %d", ct.scf_steps);


            /* Do diagonalizations if requested */
            if (((ct.scf_steps % ct.diag) == 0) && (ct.end_diag > ct.scf_steps))
                subdiag (states, vh, vnuc, vxc);

            /* Get the total energy */
            get_te (rho, rhocore, rhoc, vh, vxc, states);
        }


        if (pct.gridpe == 0)
            printf ("\n %d SCF steps were needed to reach the converegence criterion",
                    ct.scf_steps);
#endif

        quench (states, vxc, vh, vnuc, rho, rho_oppo, rhocore, rhoc);



#if 0
        if (ct.meanres > 1.0e-7 && (ct.md_steps % 10) == 0)
        {
            int steps;
            sortpsi (states);
#if GAMMA_PT
            subdiag_gamma (states, vh, vnuc, vxc);
#else
            subdiag_nongamma (states, vh, vnuc, vxc);
#endif
            for (steps = 0; steps < 7; steps++)
            {
                scf (states, vxc, vh, vnuc, rho, rhocore, rhoc);


                /* Get the total energy */
                get_te (rho, rhocore, rhoc, vh, vxc, states);
            }
        }
#endif

        ct.mix = tmix;
        ct.prjmix = tprjmix;



        /* zero out the non-moving atoms */
        for (ion = 0; ion < ct.num_ions; ion++)
        {

            if (ct.ions[ion].movable == 0)
            {

                iptr = &ct.ions[ion];

                for (ic = 0; ic < 4; ic++)
                {
                    iptr->force[ic][0] = 0.0;
                    iptr->force[ic][1] = 0.0;
                    iptr->force[ic][2] = 0.0;
                }
                iptr->velocity[0] = 0.0;
                iptr->velocity[1] = 0.0;
                iptr->velocity[2] = 0.0;
            }

        }


        /* Do another halfstep update of the velocities */
        /* to update them to the next time step */
        velup2 ();

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
        iontemp = ct.ionke * 2.0 / (3.0 * (rmg_double_t) N * kB);


        /*write data to output file */
        if (ct.checkpoint)
            if ( ct.md_steps % ct.checkpoint == 0 )
            {
                write_restart (ct.outfile, vh, rho, rho_oppo, vxc, states);
                if (pct.gridpe == 0)
                    printf ("\n Writing data to output file ...\n");
            }

        if (pct.gridpe == 0)
        {
            switch (ct.forceflag)
            {

            case MD_CVE:
                printf ("\n @CVE %5d  %15.10f  %15.10f  %15.10f  %15.10f  %10.8e",
                        ct.md_steps, ct.TOTAL, ct.ionke, ct.TOTAL + ct.ionke, iontemp, trms);
                break;
            case MD_CVT:
                if (ct.tcontrol == T_NOSE_CHAIN)
                {
                    printf ("\n @CVT-NOSE %5d  %15.10f  %15.10f  %15.10f  %15.10f  %15.10f  %10.8e",
                            ct.md_steps, ct.TOTAL, ct.ionke, nosekin + nosepot,
                            ct.TOTAL + ct.ionke + nosekin + nosepot, iontemp, trms);
                }
                if (ct.tcontrol == T_AND_SCALE)
                {
                    printf ("\n @CVT-ANDERSON %5d  %15.10f  %15.10f  %15.10f  %15.10f  %10.8e",
                            ct.md_steps, ct.TOTAL, ct.ionke, ct.TOTAL + ct.ionke, iontemp, trms);
                }
                break;
            case MD_CPT:
                printf ("\n not programed yet");
                break;

            default:
                error_handler ("Undefined output type");
            }                   /* end of switch */

            stepcount++;


            printf ("\n Total number of SCF steps so far %d", ct.total_scf_steps);

        }                       /* end of if */

        fflush (NULL);
    }                           /* end for ct.md_steps */


    if (pct.gridpe == 0)
        printf ("\n Total number of SCF steps %d", ct.total_scf_steps);



    /*Get some memory */
    my_free (crdsx);
    my_free (crdsy);
    my_free (crdsz);


}                               /* end moldyn */


void init_nose ()
{
    int ion, N, jc;
    ION *iptr;
    rmg_double_t wNose, tau_nose, kB, mass, step;
    rmg_double_t inittemp, nosesteps;
    rmg_double_t v1, v2, v3;

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
    for (ion = 0; ion < ct.num_ions; ion++)
    {

        /* Get ion pointer */
        iptr = &ct.ions[ion];

        /* Get ionic mass */
        mass = ct.sp[iptr->species].atomic_mass * mu_me;

        v1 = iptr->velocity[0];
        v2 = iptr->velocity[1];
        v3 = iptr->velocity[2];

        ct.ionke += 0.5 * mass * (v1 * v1 + v2 * v2 + v3 * v3);

    }

    inittemp = ct.ionke * 2.0 / (3.0 * (rmg_double_t) N * kB);

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
        printf ("\n Nose Frequency (THz) = %14.10f ", ct.nose.fNose);
        printf ("\n Steps/Oscillation    = %14.10f ", nosesteps);
        printf ("\n Target Temp          = %14.10f ", ct.nose.temp);
        printf ("\n Initial Temp         = %14.10f ", inittemp);
        printf ("\n i        xx        xv         xq     xf ");
        for (jc = 0; jc < ct.nose.m; jc++)
        {
            printf ("\n %d %14.10f %14.10f %16.10f %14.10f", jc,
                    ct.nose.xx[jc], ct.nose.xv[jc], ct.nose.xq[jc], ct.nose.xf[ct.fpt[0]][jc]);
        }
        printf ("\n ==============================================");
    }

}                               /* end of init_nose */


void velup1 ()
{
    int ion, ic;
    ION *iptr;
    rmg_double_t step, mass, kB;
    rmg_double_t t1, t2, v1, v2, v3;
    rmg_double_t scale = 1.0;
    rmg_double_t temperature;

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
            for (ion = 0; ion < ct.num_ions; ion++)
            {

                /* Get ion pointer */
                iptr = &ct.ions[ion];

                /* Get ionic mass */
                mass = ct.sp[iptr->species].atomic_mass * mu_me;

                /* calculate kinetic energy */
                v1 = iptr->velocity[0];
                v2 = iptr->velocity[1];
                v3 = iptr->velocity[2];

                ct.ionke += 0.5 * mass * (v1 * v1 + v2 * v2 + v3 * v3);

            }

            temperature = ct.ionke * 2.0 / (3.0 * ct.nose.N * kB);

            if (ct.ionke > 0.0)
            {
                scale = sqrt (ct.nose.temp / temperature);
                if (pct.gridpe == 0)
                    printf ("\ntscale=%f\n", scale);
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
    for (ion = 0; ion < ct.num_ions; ion++)
    {

        /* Get ion pointer */
        iptr = &ct.ions[ion];

        if (iptr->movable)
        {

            /* Get ionic mass */
            mass = ct.sp[iptr->species].atomic_mass * mu_me;

            /* update velocity by one-half timestep */
            switch (ct.mdorder)
            {

            case ORDER_2:
                t1 = step / mass;
                t2 = t1 / 2.0;

                /* step velocities forward */
                for (ic = 0; ic < 3; ic++)
                {

                    iptr->velocity[ic] *= scale;
                    iptr->velocity[ic] += t2 * iptr->force[0][ic];

                }
                break;

            case ORDER_3:
                t1 = step / mass;
                t2 = t1 / 6.0;

                /* step velocities forward */
                for (ic = 0; ic < 3; ic++)
                {

                    iptr->velocity[ic] *= scale;
                    iptr->velocity[ic] += t2 * (4.0 * iptr->force[ct.fpt[0]][ic]
                                                - 1.0 * iptr->force[ct.fpt[1]][ic]);

                }
                break;

            case ORDER_5:
                t1 = step / mass;
                t2 = t1 / 360.0;

                /* step velocities forward */
                for (ic = 0; ic < 3; ic++)
                {

                    iptr->velocity[ic] *= scale;
                    iptr->velocity[ic] += t2 * (323.0 * iptr->force[ct.fpt[0]][ic]
                                                - 264.0 * iptr->force[ct.fpt[1]][ic]
                                                + 159.0 * iptr->force[ct.fpt[2]][ic]
                                                - 38.0 * iptr->force[ct.fpt[3]][ic]);

                }
                break;

            }                   /* end of switch */

            v1 = iptr->velocity[0];
            v2 = iptr->velocity[1];
            v3 = iptr->velocity[2];

            ct.ionke += 0.5 * mass * (v1 * v1 + v2 * v2 + v3 * v3);

        }                       /* if */

    }                           /* end for */

}                               /* end of velup1 */


#if 0
/* Function move_ions takes care of mving ions, it also keeps ionic history and other stuff*/
void posup ()
{
    int ion;
    ION *iptr;
    rmg_double_t step;

    step = ct.iondt;

    /* update nose thermostats */
    if (ct.forceflag == MD_CVT && ct.tcontrol == T_NOSE_CHAIN)
        nose_posup ();

    /* Loop over ions */
    for (ion = 0; ion < ct.num_ions; ion++)
    {

        /* Get ion pointer */
        iptr = &ct.ions[ion];

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
    int ion, ic;
    ION *iptr;
    rmg_double_t step, mass;
    rmg_double_t t1, t2;
    rmg_double_t v1, v2, v3;
    rmg_double_t scale = 1.0, kB;

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
    for (ion = 0; ion < ct.num_ions; ion++)
    {

        /* Get ion pointer */
        iptr = &ct.ions[ion];

        /* Get ionic mass */
        mass = ct.sp[iptr->species].atomic_mass * mu_me;


        if (iptr->movable)
        {


            /* update velocity by one-half step */
            switch (ct.mdorder)
            {

            case ORDER_2:
                t1 = step / mass;
                t2 = t1 / 2.0;

                for (ic = 0; ic < 3; ic++)
                {
                    iptr->velocity[ic] += t2 * iptr->force[0][ic];
                    iptr->velocity[ic] *= scale;
                }
                break;

            case ORDER_3:
                t1 = step / mass;
                t2 = t1 / 6.0;

                for (ic = 0; ic < 3; ic++)
                {
                    iptr->velocity[ic] += t2 * (2.0 * iptr->force[ct.fpt[0]][ic]
                                                + 1.0 * iptr->force[ct.fpt[1]][ic]);
                    iptr->velocity[ic] *= scale;
                }
                break;

            case ORDER_5:
                t1 = step / mass;
                t2 = t1 / 360.0;

                for (ic = 0; ic < 3; ic++)
                {
                    iptr->velocity[ic] += t2 * (97.0 * iptr->force[ct.fpt[0]][ic]
                                                + 114.0 * iptr->force[ct.fpt[1]][ic]
                                                - 39.0 * iptr->force[ct.fpt[2]][ic]
                                                + 8.0 * iptr->force[ct.fpt[3]][ic]);
                    iptr->velocity[ic] *= scale;
                }
                break;

            }                   /* end of switch */

        }                       /* end of if */

        v1 = iptr->velocity[0];
        v2 = iptr->velocity[1];
        v3 = iptr->velocity[2];

        ct.ionke += 0.5 * mass * (v1 * v1 + v2 * v2 + v3 * v3);

    }                           /* end for */

    /* perform second half update of nose coords */
    if (ct.forceflag == MD_CVT && ct.tcontrol == T_NOSE_CHAIN)
    {
        nose_velup2 ();
    }

}                               /* end of velup2 */





void rms_disp (rmg_double_t * rms, rmg_double_t * trms)
{

    int it;
    ION *iptr;
    rmg_double_t t1, t2, t3;

    *trms = 0.0;
    rms[0] = rms[1] = rms[2] = 0.0;
    for (it = 0; it < ct.num_ions; it++)
    {

        iptr = &ct.ions[it];

        t1 = iptr->xtal[0] - iptr->ixtal[0];
        if (t1 < 0.0)
            t1 *= -1.0;
        if (t1 > 0.5)
            t1 = 1.0 - t1;

        t2 = iptr->xtal[1] - iptr->ixtal[1];
        if (t2 < 0.0)
            t2 *= -1.0;
        if (t2 > 0.5)
            t2 = 1.0 - t2;

        t3 = iptr->xtal[2] - iptr->ixtal[2];
        if (t3 < 0.0)
            t3 *= -1.0;
        if (t3 > 0.5)
            t3 = 1.0 - t3;

        rms[0] = t1;
        rms[1] = t2;
        rms[2] = t3;

        *trms += metric (rms) / ((rmg_double_t) ct.nose.N);

    }

}                               /* end of rms_disp */

void nose_velup1 ()
{
    int jc;
    rmg_double_t step;
    rmg_double_t scale;
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
    rmg_double_t step;

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
    rmg_double_t kB, step;
    rmg_double_t scale;
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


void nose_energy (rmg_double_t * nosekin, rmg_double_t * nosepot)
{
    int jc;
    rmg_double_t kB;

    /* define Boltzmann param */
    kB = 1.0 / (11605.0 * Ha_eV);

    /* calculate total nose energy */
    *nosekin = 0.5 * ct.nose.xq[0] * ct.nose.xv[0] * ct.nose.xv[0];
    *nosepot = 2.0 * ct.nose.k0 * ct.nose.xx[0];
#if 0
    if (pct.gridpe == 0)
        printf ("\n @therm%d %14.10f %14.10f %14.10f %14.10f %14.10f %14.10f",
                0, ct.nose.xx[0], ct.nose.xv[0], ct.nose.xf[ct.fpt[0]][0],
                ct.nose.xf[ct.fpt[1]][0], ct.nose.xf[ct.fpt[2]][0], ct.nose.xf[ct.fpt[3]][0]);
#endif
    for (jc = 1; jc < ct.nose.m; jc++)
    {
        *nosekin += 0.5 * ct.nose.xq[jc] * ct.nose.xv[jc] * ct.nose.xv[jc];
        *nosepot += ct.nose.temp * kB * ct.nose.xx[jc];
#if 0
        if (pct.gridpe == 0)
            printf ("\n @therm%d %14.10f %14.10f %14.10f %14.10f %14.10f %14.10f",
                    jc, ct.nose.xx[jc], ct.nose.xv[jc], ct.nose.xf[ct.fpt[0]][jc],
                    ct.nose.xf[ct.fpt[1]][jc],
                    ct.nose.xf[ct.fpt[2]][jc], ct.nose.xf[ct.fpt[3]][jc]);
#endif

    }

#if 0
    if (pct.gridpe == 0)
        printf ("\n nose_energy: %20.10f %20.10f %20.10f", *nosekin, *nosepot, *nosekin + *nosepot);
#endif

}                               /* end of nose_energy */

/******/
