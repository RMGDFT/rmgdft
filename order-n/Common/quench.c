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
 *   states: point to orbital structure (see md.h)
 *   vxc:    exchange correlation potential
 *   vh:     Hartree potential
 *   vnuc:   Pseudopotential 
 *   rho:    total valence charge density
 *   rhocore: core chare density only for non-linear core correction
 *   rhoc:   Gaussian compensating charge density
 * OUTPUT
 *   states, vxc, vh, rho are updated
 * PARENTS
 *   cdfastrlx.c fastrlx.c md.c
 * CHILDREN
 *   scf.c force.c get_te.c subdiag.c get_ke.c
 * SOURCE
 */

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "md.h"

void quench(STATE * states, STATE * states1, REAL * vxc, REAL * vh,
            REAL * vnuc, REAL * vh_old, REAL * vxc_old, REAL * rho, REAL * rhoc, REAL * rhocore)
{
    int outcount = 0;
    static int CONVERGENCE = FALSE;
    double time1, time2;

    time1 = my_crtc();


    for (ct.scf_steps = 0; ct.scf_steps < ct.max_scf_steps; ct.scf_steps++)
    {
        if (pct.gridpe == 0)
            printf("\n\n\n ITERATION     %d\n", ct.scf_steps);

        /* Perform a single self-consistent step */
        if (!CONVERGENCE)
            scf(states, states1, vxc, vh, vnuc, rho, rhoc, rhocore, vxc_old, vh_old, &CONVERGENCE);

        if (pct.gridpe == 0)
            write_eigs(states);

        if (CONVERGENCE)
        {
            if (pct.gridpe == 0)
                printf ("\n\n Convergence has been achieved. stopping ...\n");
            break;
        }


        /* Check if we need to output intermediate results */
        if (outcount >= ct.outcount)
        {
            if (pct.gridpe == 0)
                printf("\n TOTAL ENERGY = %14.7f\n", ct.TOTAL);
            outcount = 0;
        }
        outcount++;
    }

    get_mat_Omega(states, mat_Omega);
    /* Calculate the force */
    force(rho, rhoc, vh, vxc, vnuc, states); 
    /* write out the force */
   if (pct.gridpe == 0) write_force();


    time2 = my_crtc();
    rmg_timings(QUENCH_TIME, time2 - time1);


}
