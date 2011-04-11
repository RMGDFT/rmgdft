/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

/****f* QMD-MGDFT/init_wf.c *****
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
 *   void init_wf(STATE *states)
 *   Generates an initial set of wavefunctions randomly. 
 *   This function only generates one set of initial wavefunctions
 *   so if multiple k-points are being using in the calculations
 *   then it must be called explicitly for each k-point. 
 *   If ct.occflag =/ 0, start with each state equally occupied  
 * INPUTS
 *   nothing
 * OUTPUT
 *   states: random initial wave function
 * PARENTS
 *   init.c
 * CHILDREN
 *   scatter_psi.c norm_psi.c pe2xyz.c pack_ptos.c app_cir.c
 * SOURCE
 */



#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "main.h"


void init_wf (STATE * states)
{

    int idx, state, ix, iy, iz, pbasis, nspin=(pct.spin_flag+1);
    int xoff, yoff, zoff;
    REAL *tmp_psiR, *tmp_psiI;
    STATE *sp;
    REAL xrand[2 * NX_GRID], yrand[2 * NY_GRID], zrand[2 * NZ_GRID];
    long idum;

    pbasis = PX0_GRID * PY0_GRID * PZ0_GRID;

    my_malloc (tmp_psiR, pbasis, REAL);
#if !GAMMA_PT
    my_malloc (tmp_psiI, pbasis, REAL);
#endif

    /* Set state 0 to a constant */
    sp = states;
    for (idx = 0; idx < P0_BASIS; idx++)
        tmp_psiR[idx] = ONE;
    scatter_psi (tmp_psiR, NULL, sp, 0);

    /* If not a gamma point calculation do the imaginary part */
#if !GAMMA_PT
    scatter_psi (NULL, tmp_psiR, sp, 0);
#endif

    
    /* If random start and Fermi occupation, start with
       each state equally occupied  */

    if (ct.occ_flag && (ct.runflag == 0))
    {  
	/* Set occupation for the first state */    
	for (idx = 0; idx < nspin; idx++)
	    sp->occupation[idx] = ct.nel / (nspin * ct.num_states);	
    }


    pe2xyz (pct.gridpe, &ix, &iy, &iz);
    xoff = ix * PX0_GRID;
    yoff = iy * PY0_GRID;
    zoff = iz * PZ0_GRID;


    /* Initialize the random number generator */
    idum = 3356;
    rand0 (&idum);
    for (state = 1; state < ct.num_states; state++)
    {


        /* Generate x, y, z random number sequences */
        for (idx = 0; idx < ct.psi_nxgrid; idx++)
            xrand[idx] = rand0 (&idum) - 0.5;
        for (idx = 0; idx < ct.psi_nygrid; idx++)
            yrand[idx] = rand0 (&idum) - 0.5;
        for (idx = 0; idx < ct.psi_nzgrid; idx++)
            zrand[idx] = rand0 (&idum) - 0.5;

#if !GAMMA_PT
        for (idx = ct.psi_nxgrid; idx < 2 * ct.psi_nxgrid; idx++)
            xrand[idx] = rand0 (&idum) - 0.5;
        for (idx = ct.psi_nygrid; idx < 2 * ct.psi_nygrid; idx++)
            yrand[idx] = rand0 (&idum) - 0.5;
        for (idx = ct.psi_nzgrid; idx < 2 * ct.psi_nzgrid; idx++)
            zrand[idx] = rand0 (&idum) - 0.5;
#endif

        sp = &states[state];
        
	/* If random start and Fermi occupation, start with
           each state equally occupied  */
    	
	if (ct.occ_flag && (ct.runflag == 0))
    	{  
		for (idx = 0; idx < nspin; idx++)
	    		sp->occupation[idx] = ct.nel / ( nspin * ct.num_states );	
    	}



        idx = 0;
        for (ix = 0; ix < PX0_GRID; ix++)
        {

            for (iy = 0; iy < PY0_GRID; iy++)
            {

                for (iz = 0; iz < PZ0_GRID; iz++)
                {

                    tmp_psiR[idx] = xrand[xoff + ix] * yrand[yoff + iy] * zrand[zoff + iz];
                    tmp_psiR[idx] = tmp_psiR[idx] * tmp_psiR[idx];
#if !GAMMA_PT
                    tmp_psiI[idx] =
                        xrand[ct.psi_nxgrid + xoff + ix] * yrand[ct.psi_nygrid + yoff +
                                                                 iy] * zrand[ct.psi_nzgrid + zoff +
                                                                             iz];
                    tmp_psiI[idx] = tmp_psiI[idx] * tmp_psiI[idx];

#endif
                    idx++;

                }               /* end for */
            }                   /* end for */
        }                       /* end for */

        scatter_psi (tmp_psiR, NULL, sp, 0);

#if !GAMMA_PT
        scatter_psi (NULL, tmp_psiI, sp, 0);
#endif

    }                           /* end for */

    my_free (tmp_psiR);
#if !GAMMA_PT
    my_free (tmp_psiI);
#endif

}                               /* end init_wf */

/******/
