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
#include "main.h"

#define SMALL 1.e-35

/*For a quantity localized around ionic positions, this function finds radius in number of grid points
 * given a radius in a.u.*/
int radius2grid (rmg_double_t radius, rmg_double_t mingrid_spacing)
{

    rmg_double_t scale, t1, t2;
    int it1, dim;
    
    /* Set the scaling factor for determining the radius of the local grids */
    scale = 1.0;
    if (ct.ibrav == CUBIC_BC)
	scale = 1.1;
    if (ct.ibrav == CUBIC_FC)
	scale = 1.3;
        
	t1 = 2.0 * scale * radius / mingrid_spacing;
        t1 = modf (t1, &t2);
        it1 = (int) t2;
        if (t1 > 0.5)
            it1++;
        if (!(it1 % 2))
            it1++;
	dim = it1;

	return dim;
}


void init_psp (void)
{

    int isp, idx, ip, write_flag;
    SPECIES *sp;
    rmg_double_t *work, *workr, Zv, rc, rfil, t1;
    rmg_double_t rcut, exp_fac;
    char newname[MAX_PATH];
    FILE *psp = NULL;
    FILE *psp2 = NULL;


    my_malloc (work, MAX_RGRID + MAX_LOCAL_LIG, rmg_double_t);
    workr = work + MAX_RGRID;

    write_flag = 0;
    if (verify ("write_pseudopotential_plots", &SET))
        write_flag = 1;

    /*Initialize max_nlpoints and max_nlfpoints */
    ct.max_nlpoints = 0;
    ct.max_nlfpoints = 0;


    /* Loop over species */
    for (isp = 0; isp < ct.num_species; isp++)
    {
        if (pct.gridpe == 0 && write_flag)
        {
            snprintf (newname, MAX_PATH, "%s.local%d.xmgr", ct.basename, isp);
            my_fopen (psp, newname, "w+");
        }
        sp = &ct.sp[isp];

        /*Get ldim */
	sp->ldim = radius2grid (sp->lradius, ct.hmingrid/ (rmg_double_t) FG_NX);
        if ((sp->ldim >= ct.psi_fnxgrid) || (sp->ldim >= ct.psi_fnygrid)
            || (sp->ldim >= ct.psi_fnzgrid))
            error_handler ("local potential radius exceeds global grid size");


        /*Get drnlig */
        /*sp->drlig = sqrt(3.0) * (sp->ldim + 1.0) * ct.hmaxgrid / 2.0 /(rmg_double_t)FG_NX; */
        t1 = sp->ldim / FG_NX + 1;
        sp->drlig = sqrt (3.0) * (t1 + 1.0) * ct.hmaxgrid / 2.0;
        if (ct.ibrav == HEXAGONAL)
            sp->drlig *= 2.0;
        t1 = (rmg_double_t) MAX_LOCAL_LIG;
        sp->drlig /= t1;


        /*Get nldim */
	sp->nldim = radius2grid (sp->nlradius, ct.hmingrid);
        sp->nlfdim = ct.nxfgrid * sp->nldim;
        
	if ((sp->nldim >= ct.psi_nxgrid) || (sp->nldim >= ct.psi_nygrid)
            || (sp->nldim >= ct.psi_nzgrid))
            error_handler ("Non-local potential radius exceeds global grid size");

        /*Get drnlig */
        sp->drnlig = sqrt (3.0) * (sp->nldim + 1.0) * ct.hmaxgrid / 2.0;
        if (ct.ibrav == HEXAGONAL)
            sp->drnlig *= 2.0;
        t1 = (rmg_double_t) MAX_LOCAL_LIG;
        sp->drnlig /= t1;



        /*ct.max_nlpoints is max of nldim*nldim*nldim for all species */
        if (ct.max_nlpoints < (sp->nldim * sp->nldim * sp->nldim))
            ct.max_nlpoints = sp->nldim * sp->nldim * sp->nldim;

        if (ct.max_nlfpoints < (sp->nlfdim * sp->nlfdim * sp->nlfdim))
            ct.max_nlfpoints = sp->nlfdim * sp->nlfdim * sp->nlfdim;


	/*Filter and interpolate local potential into fine linear grid*/
        Zv = sp->zvalence;
        rc = sp->rc;

        /* Generate the difference potential */
        for (idx = 0; idx < sp->rg_points; idx++)
            work[idx] = sp->vloc0[idx] + Zv * erf (sp->r[idx] / rc) / sp->r[idx];


        if (pct.gridpe == 0 && write_flag)
        {
            for (idx = 0; idx < sp->rg_points; idx++)
                    fprintf (psp, "%e  %e\n", sp->r[idx], work[idx]);
            
	    fprintf (psp, "\n&&\n");
        }



        /* Transform to g-space and filter it */
        /*rft1 (ct.cparm, work, &sp->r[0], sp->localig, &sp->rab[0], sp->rg_points, 0, sp->drlig,
              sp->gwidth, MAX_LOCAL_LIG);*/

	filter_potential(work, &sp->r[0], sp->rg_points, sp->lradius, 0.25, ct.cparm, 
		sp->localig, &sp->rab[0], 0, sp->drlig, sp->gwidth, MAX_LOCAL_LIG, sp->lrcut, sp->rwidth, sp->drlocalig);


	/*Write local projector into a file if requested*/
	if ((pct.gridpe == 0) && write_flag)
	{
	    rfil = 0.0;
	    for (idx = 0; idx < MAX_LOCAL_LIG; idx++)
	    {
                
		fprintf (psp, "%e  %e \n", rfil, sp->localig[idx]);
		rfil += sp->drlig;
	    }
            
	    /* output xmgr data separator */
	    fprintf (psp, "\n&&\n");
	    
	    rfil = 0.0;
	    for (idx = 0; idx < MAX_LOCAL_LIG; idx++)
	    {
                
		fprintf (psp, "%e  %e \n", rfil, sp->drlocalig[idx]);
		rfil += sp->drlig;
	    }
            
	    fclose (psp);
	}


	/*Open file for writing beta function*/
        if (pct.gridpe == 0 && write_flag)
	{
            snprintf (newname, MAX_PATH, "%s.beta%d.xmgr", ct.basename, isp);
            my_fopen (psp, newname, "w+");
            
	    snprintf (newname, MAX_PATH, "%s.drbeta%d.xmgr", ct.basename, isp);
            my_fopen (psp2, newname, "w+");
	}
	
	/* Write raw beta function into file if requested*/
        for (ip = 0; ip < sp->nbeta; ip++)
        {

            if (pct.gridpe == 0 && write_flag)
            {
                for (idx = 0; idx < sp->kkbeta; idx++)
                    fprintf (psp, "%e  %e\n", sp->r[idx], sp->beta[ip][idx]);
                fprintf (psp, "\n&&\n");
            }

	    filter_potential(&sp->beta[ip][0], &sp->r[0], sp->rg_points, sp->nlradius, 0, ct.betacparm, &sp->betalig[ip][0], 
		    &sp->rab[0], sp->llbeta[ip], sp->drnlig, sp->gwidth, MAX_LOCAL_LIG, sp->nlrcut[sp->llbeta[ip]], sp->rwidth, &sp->drbetalig[ip][0]);


	    /* Is this necessary ??? */
            if (sp->llbeta[ip])
                sp->betalig[ip][0] = 0.0;


	    /* output filtered non-local projector to a file  if requested */
	    if (pct.gridpe == 0 && write_flag)
	    {
		rfil = 0.0;
		for (idx = 0; idx < MAX_LOCAL_LIG; idx++)
		{
		    {
			fprintf (psp, "%e  %e\n", rfil, sp->betalig[ip][idx]);
			fprintf (psp2, "%e  %e\n", rfil, sp->drbetalig[ip][idx]);
		    }

		    rfil += sp->drnlig;
		}                   /* end for */
	    }

            /* output xmgr data separator */
            if (pct.gridpe == 0 && write_flag)
	    {
                fprintf (psp, "\n&&\n");
                fprintf (psp2, "\n&&\n");
	    }

        }                       /* end for ip */

        /* Now take care of the core charge if nlcc is being used */
        if (sp->nlccflag)
        {

            for (idx = 0; idx < sp->rg_points; idx++)
                work[idx] = sp->rspsco[idx] / (4.0 * PI);


	    /*Write raw (pre-filtered) data to a file if requested */
            if (pct.gridpe == 0 && write_flag)
            {
		for (idx = 0; idx < sp->rg_points; idx++)
		    fprintf (psp, "%e  %e\n", sp->r[idx], work[idx]);
		fprintf (psp, "\n&&\n");
            }

	
	    filter_potential(work, &sp->r[0], sp->rg_points, sp->lradius, 0.25, ct.cparm, &sp->rhocorelig[0], 
		    &sp->rab[0], 0, sp->drlig, sp->gwidth, MAX_LOCAL_LIG, sp->lrcut, sp->rwidth, NULL);

	    /*Oscilations at the tail end of filtered function may cause rhocore to be negative
	     * but I am not sure if this is the right solution, it may be better to fix charge density
	     * this rhocore less smooth*/
	    for (idx = 0; idx < MAX_LOCAL_LIG; idx++)
		if (sp->rhocorelig[idx] < 0.0)
		    sp->rhocorelig[idx] = 0.0;



	    /*Write filtered data to a file if requested */
	    if (pct.gridpe == 0 && write_flag)
	    {
		rfil = 0.0;
		
		for (idx = 0; idx < MAX_LOCAL_LIG; idx++)
		{
		    fprintf (psp, "%e  %e\n", rfil, sp->rhocorelig[idx]);
		    rfil += sp->drlig;
		}

		fprintf (psp, "\n&&\n");
	    }

        }                       /* end if */


        if (pct.gridpe == 0 && write_flag)
	{
            fclose (psp);
            fclose (psp2);
	}

    }                           /* end for */


    my_free (work);

}                               /* end init_kbr.c */

/******/
