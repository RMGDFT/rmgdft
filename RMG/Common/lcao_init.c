/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

/****f* QMD-MGDFT/init_psp.c *****
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
 *   void init_psp(void)
 *   Initializes radial Kleinman-Bylander projectors and
 *   calculates the normalization coefficient.
 *
 *   Also subtracts off the potential due to the compensating charges from
 *   the local potential.
 * INPUTS
 *   nothing
 * OUTPUT
 *   projectors are stored in structure SPECIES (see main.h)
 * PARENTS
 *   init.c
 * CHILDREN
 *   rft1.c radiff.c 
 * SOURCE
 */


#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "grid.h"
#include "common_prototypes.h"
#include "main.h"

#define SMALL 1.e-35

void lcao_init (void)
{

    int isp, idx, ip, it1, write_flag, ibrav;
    SPECIES *sp;
    rmg_double_t *work, rfil;
    rmg_double_t t1, t2, rcut, exp_fac;
    char newname[MAX_PATH];
    FILE *psp = NULL;

    ibrav = get_ibrav_type();

    write_flag = 0;
    if (verify ("write_pseudopotential_plots", &SET))
        write_flag = 1;



    /* Loop over species */
    for (isp = 0; isp < ct.num_species; isp++)
    {
        if (pct.gridpe == 0 && write_flag)
        {
            snprintf (newname, MAX_PATH, "%s.awf%d.xmgr", ct.basename, isp);
            my_fopen (psp, newname, "w+");
        }
        sp = &ct.sp[isp];

        work = sp->atomic_rho;

	/*Find radii of atomic wavefunctions and charge in terms of number of grid points*/
	sp->adim_rho  = radius2grid (sp->aradius, ct.hmingrid/ (rmg_double_t) get_FG_NX());
        if ((sp->adim_rho >= ct.psi_fnxgrid) || (sp->adim_rho >= ct.psi_fnygrid)
            || (sp->adim_rho >= ct.psi_fnzgrid))
            error_handler ("LCAO charge radius exceeds global grid size");
	
	sp->adim_wave = radius2grid (sp->aradius, ct.hmingrid);
	if ((sp->adim_wave >= ct.psi_nxgrid) || (sp->adim_wave >= ct.psi_nygrid)
            || (sp->adim_wave >= ct.psi_nzgrid))
            error_handler ("LCAO wavefunctions radius exceeds global grid size");
	

	sp->drlig_awave = sqrt (3.0) * (sp->adim_wave + 1.0) * ct.hmaxgrid / 2.0;
	if (ibrav == HEXAGONAL)
	     sp->drlig_awave *= 2.0;
	t1 = (rmg_double_t) MAX_LOCAL_LIG;
	sp->drlig_awave /= t1;

	sp->drlig_arho = sqrt (3.0) * (sp->adim_rho/get_FG_NX() + 1.0) * ct.hmaxgrid / 2.0;
	if (ibrav == HEXAGONAL)
	     sp->drlig_arho *= 2.0;
	t1 = (rmg_double_t) MAX_LOCAL_LIG;
	sp->drlig_arho /= t1;

        if (pct.gridpe == 0 && write_flag)
        {
            for (idx = 0; idx < sp->rg_points; idx++)
            {
                    fprintf (psp, "%15.8f  %15.8f\n", sp->r[idx], work[idx]);
            }
            fprintf (psp, "\n&&\n");
        }

	filter_potential(work, &sp->r[0], sp->rg_points, sp->aradius, 0.25, ct.cparm, &sp->arho_lig[0],
		&sp->rab[0], 0, sp->drlig_arho, sp->agwidth, MAX_LOCAL_LIG, sp->acut, sp->arwidth, NULL);

#if 0
        /* Transform to g-space and filter it */
        /*May need ct.qcparm here */
        rft1 (ct.cparm, work, &sp->r[0], sp->arho_lig, &sp->rab[0], sp->rg_points, 0, sp->drlig_arho,
              sp->agwidth, MAX_LOCAL_LIG);



        /* Fix up the first point */
        sp->arho_lig[0] = 2.0 * sp->arho_lig[1] - sp->arho_lig[2];

        /* Now cut it off smoothly in real space */
        rcut = sp->acut;

        rfil = 0.0;
        for (idx = 0; idx < MAX_LOCAL_LIG; idx++)
        {

            if (rfil > rcut)
            {

                t1 = (rfil - rcut) / rcut;
		exp_fac = exp (-sp->arwidth * t1 * t1);
                
		sp->arho_lig[idx] *= exp_fac;

                if (fabs (sp->arho_lig[idx]) < SMALL)
                    sp->arho_lig[idx] = 0.; 

            }                   /* end if */

            /* output local projector */
            if (pct.gridpe == 0 && write_flag)
            {
                    fprintf (psp, "%15.8f  %15.8f\n", rfil, sp->arho_lig[idx]);

            }                   /* endif */

            rfil += sp->drlig_arho;
        }                       /* end for idx */
#endif
            
	if (pct.gridpe == 0 && write_flag)
	{
	    rfil = 0.0;
	    for (idx = 0; idx < MAX_LOCAL_LIG; idx++)
	    {
		fprintf (psp, "%15.8f  %15.8f\n", rfil, sp->arho_lig[idx]);
		rfil += sp->drlig_arho;
	    }
	}


        /* output xmgr data separator */
        if (pct.gridpe == 0 && write_flag)
            fprintf (psp, "\n&&\n");




	/*Now get wavefunctions on linear grid*/
        for (ip = 0; ip < sp->num_atomic_waves; ip++)
        {

            if (pct.gridpe == 0 && write_flag)
            {
                for (idx = 0; idx < sp->rg_points; idx++)
                    fprintf (psp, "%15.8f  %15.8f\n", sp->r[idx], sp->atomic_wave[ip][idx]);
                fprintf (psp, "\n&&\n");
            }
	    
	    filter_potential(&sp->atomic_wave[ip][0], &sp->r[0], sp->rg_points, sp->aradius, 0.25, ct.betacparm, &sp->awave_lig[ip][0],
		    &sp->rab[0], sp->atomic_wave_l[ip], sp->drlig_awave, sp->agwidth, MAX_LOCAL_LIG, sp->acut, sp->arwidth, NULL);

#if 0
            /* Transform to g-space and filter it */
            rft1 (ct.betacparm, &sp->atomic_wave[ip][0], &sp->r[0], &sp->awave_lig[ip][0],
                  &sp->rab[0], sp->rg_points, sp->atomic_wave_l[ip], sp->drlig_awave, sp->agwidth,
                  MAX_LOCAL_LIG);


            /* Fix up the first point */
            sp->awave_lig[ip][0] = 2.0 * sp->awave_lig[ip][1] - sp->awave_lig[ip][2];


            /* Now cut it off smoothly in real space */
            rcut = sp->acut;

            rfil = 0.0;
            for (idx = 0; idx < MAX_LOCAL_LIG; idx++)
            {

                if (rfil > rcut)
                {

                    t1 = (rfil - rcut) / rcut;
		    exp_fac = exp (-sp->arwidth * t1 * t1);
		    
                    sp->awave_lig[ip][idx] *= exp_fac;

                    if (fabs (sp->awave_lig[ip][idx]) < SMALL)
                        sp->awave_lig[ip][idx] = 0.;

                }               /* end if */

                /* output non-local projector */
                if (pct.gridpe == 0 && write_flag)
                    fprintf (psp, "%15.8f  %15.8f\n", rfil, sp->awave_lig[ip][idx]);

                rfil += sp->drlig_awave;

            }                   /* end for */

#endif
	if (pct.gridpe == 0 && write_flag)
	{
	    rfil = 0.0;
	    for (idx = 0; idx < MAX_LOCAL_LIG; idx++)
	    {
		fprintf (psp, "%15.8f  %15.8f\n", rfil, sp->awave_lig[ip][idx]);
		rfil += sp->drlig_arho;
	    }
	}

            /* output xmgr data separator */
            if (pct.gridpe == 0 && write_flag)
                fprintf (psp, "\n&&\n");

        }                       /* end for ip */

        if (pct.gridpe == 0 && write_flag)
            fclose (psp);

    }                           /* end for */


    //my_free (work);

}                               /* end init_kbr.c */

/******/
