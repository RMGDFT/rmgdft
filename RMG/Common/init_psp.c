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
int radius2grid (REAL radius, REAL mingrid_spacing)
{

    REAL scale, t1, t2;
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
    REAL *work, *workr, Zv, rc, rfil, t1;
    REAL rcut, exp_fac;
    char newname[MAX_PATH];
    FILE *psp = NULL;


    my_malloc (work, MAX_RGRID + MAX_LOCAL_LIG, REAL);
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
            snprintf (newname, MAX_PATH, "%s.beta%d.xmgr", ct.basename, isp);
            my_fopen (psp, newname, "w+");
        }
        sp = &ct.sp[isp];

        /*Get ldim */
	sp->ldim = radius2grid (sp->lradius, ct.hmingrid/ (REAL) FG_NX);
        if ((sp->ldim >= ct.psi_fnxgrid) || (sp->ldim >= ct.psi_fnygrid)
            || (sp->ldim >= ct.psi_fnzgrid))
            error_handler ("local potential radius exceeds global grid size");


        /*Get drnlig */
        /*sp->drlig = sqrt(3.0) * (sp->ldim + 1.0) * ct.hmaxgrid / 2.0 /(REAL)FG_NX; */
        t1 = sp->ldim / FG_NX + 1;
        sp->drlig = sqrt (3.0) * (t1 + 1.0) * ct.hmaxgrid / 2.0;
        if (ct.ibrav == HEXAGONAL)
            sp->drlig *= 2.0;
        t1 = (REAL) MAX_LOCAL_LIG;
        sp->drlig /= t1;


        /*Get nldim */
	sp->nldim = radius2grid (sp->nlradius, ct.hmingrid/ (REAL) FG_NX);
        sp->nlfdim = ct.nxfgrid * sp->nldim;
        
	if ((sp->nldim >= ct.psi_nxgrid) || (sp->nldim >= ct.psi_nygrid)
            || (sp->nldim >= ct.psi_nzgrid))
            error_handler ("Non-local potential radius exceeds global grid size");

        /*Get drnlig */
        sp->drnlig = sqrt (3.0) * (sp->nldim + 1.0) * ct.hmaxgrid / 2.0;
        if (ct.ibrav == HEXAGONAL)
            sp->drnlig *= 2.0;
        t1 = (REAL) MAX_LOCAL_LIG;
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
            {
                if (sp->r[idx] < sp->lradius)
                    fprintf (psp, "%15.8f  %15.8f\n", sp->r[idx], work[idx]);
            }
            fprintf (psp, "\n&&\n");
        }



        /* Transform to g-space and filter it */
        rft1 (ct.cparm, work, &sp->r[0], sp->localig, &sp->rab[0], sp->rg_points, 0, sp->drlig,
              sp->gwidth, MAX_LOCAL_LIG);


        /* Evaluate it's radial derivative */
        for (idx = 0; idx < MAX_LOCAL_LIG; idx++)
            workr[idx] = sp->drlig * ((REAL) idx);

        radiff (sp->localig, sp->drlocalig, workr, MAX_LOCAL_LIG, 0.0);


        /* Fix up the first point */
        sp->localig[0] = 2.0 * sp->localig[1] - sp->localig[2];
        sp->drlocalig[1] = 2.0 * sp->drlocalig[2] - sp->drlocalig[3];
        sp->drlocalig[0] = 2.0 * sp->drlocalig[1] - sp->drlocalig[2];


        /* Now cut it off smoothly in real space */
        rcut = sp->lrcut;

        rfil = 0.0;
        for (idx = 0; idx < MAX_LOCAL_LIG; idx++)
        {

            if (rfil > rcut)
            {

                t1 = (rfil - rcut) / rcut;
		exp_fac = exp (-sp->rwidth * t1 * t1);

                sp->localig[idx] *= exp_fac;
                sp->drlocalig[idx] *= exp_fac;

                if (fabs (sp->localig[idx]) < SMALL)
                    sp->localig[idx] = 0.;
                if (fabs (sp->drlocalig[idx]) < SMALL)
                    sp->drlocalig[idx] = 0.;

            }                   /* end if */

            /* output local projector */
            if (pct.gridpe == 0 && write_flag)
            {
                if (rfil < sp->lradius)
                    fprintf (psp, "%1.8f  %15.8f  %15.8f\n", rfil, sp->localig[idx],
                             sp->drlocalig[idx]);
            }                   /* endif */

            rfil += sp->drlig;
        }                       /* end for idx */


        /* output xmgr data separator */
        if (pct.gridpe == 0 && write_flag)
            fprintf (psp, "\n&&\n");



	/*Filter and interpolate beta functions into fine linear grid*/
        for (ip = 0; ip < sp->nbeta; ip++)
        {

            if (pct.gridpe == 0 && write_flag)
            {
                for (idx = 0; idx < sp->kkbeta; idx++)
                    fprintf (psp, "%15.8f  %15.8f\n", sp->r[idx], sp->beta[ip][idx]);
                fprintf (psp, "\n&&\n");
            }

            /* Transform to g-space and filter it */
            rft1 (ct.betacparm, &sp->beta[ip][0], &sp->r[0], &sp->betalig[ip][0],
                  &sp->rab[0], sp->rg_points, sp->llbeta[ip], sp->drnlig, sp->gwidth,
                  MAX_LOCAL_LIG);


            /* Evaluate it's radial derivative */
            for (idx = 0; idx < MAX_LOCAL_LIG; idx++)
                workr[idx] = sp->drnlig * ((REAL) idx);

            radiff (&sp->betalig[ip][0], &sp->drbetalig[ip][0], workr, MAX_LOCAL_LIG, 0.0);


            /* Fix up the first point */
            sp->betalig[ip][0] = 2.0 * sp->betalig[ip][1] - sp->betalig[ip][2];
            if (sp->llbeta[ip])
                sp->betalig[ip][0] = 0.0;

            sp->drbetalig[ip][1] = 2.0 * sp->drbetalig[ip][2] - sp->drbetalig[ip][3];
            sp->drbetalig[ip][0] = 2.0 * sp->drbetalig[ip][1] - sp->drbetalig[ip][2];


            /* Now cut it off smoothly in real space */
            rcut = sp->nlrcut[sp->llbeta[ip]];

            rfil = 0.0;
            for (idx = 0; idx < MAX_LOCAL_LIG; idx++)
            {

                if (rfil > rcut)
                {

                    t1 = (rfil - rcut) / rcut;
		    exp_fac = exp (-sp->rwidth * t1 * t1);
		    
                    sp->betalig[ip][idx] *= exp_fac;
                    sp->drbetalig[ip][idx] *= exp_fac;

                    if (fabs (sp->betalig[ip][idx]) < SMALL)
                        sp->betalig[ip][idx] = 0.;
                    if (fabs (sp->drbetalig[ip][idx]) < SMALL)
                        sp->drbetalig[ip][idx] = 0.;


                }               /* end if */

                /* output non-local projector */
                if (pct.gridpe == 0 && write_flag)
                    fprintf (psp, "%15.8f  %15.8f %15.8f\n", rfil, sp->betalig[ip][idx],
                             sp->drbetalig[ip][idx]);

                rfil += sp->drnlig;

            }                   /* end for */

            /* output xmgr data separator */
            if (pct.gridpe == 0 && write_flag)
                fprintf (psp, "\n&&\n");

        }                       /* end for ip */

        /* Now take care of the core charge if nlcc is being used */
        if (sp->nlccflag)
        {

            for (idx = 0; idx < sp->rg_points; idx++)
                work[idx] = sp->rspsco[idx] / (4.0 * PI);


            if (pct.gridpe == 0 && write_flag)
            {
                for (idx = 0; idx < sp->rg_points; idx++)
                {
                    if (sp->r[idx] < sp->lrcut)
                        fprintf (psp, "%15.8f  %15.8f\n", sp->r[idx], work[idx]);
                }
                fprintf (psp, "\n&&\n");
            }


            /* Transform to g-space and filter it */
            rft1 (ct.cparm, work, &sp->r[0], &sp->rhocorelig[0],
                  &sp->rab[0], sp->rg_points, 0, sp->drlig, sp->gwidth, MAX_LOCAL_LIG);


            /* Fix up the first point */
            sp->rhocorelig[0] = 2.0 * sp->rhocorelig[1] - sp->rhocorelig[2];


            /* Now cut it off smoothly in real space */
            rcut = sp->lrcut;

            rfil = 0.0;
            for (idx = 0; idx < MAX_LOCAL_LIG; idx++)
            {

                if (rfil > rcut)
                {

                    t1 = (rfil - rcut) / rcut;
		    exp_fac = exp (-sp->rwidth * t1 * t1);

                    sp->rhocorelig[idx] *= exp_fac;

                    if (fabs (sp->rhocorelig[idx]) < SMALL)
                        sp->rhocorelig[idx] = 0.;

                }               /* end if */

		/*Can this happen ???*/
                if (sp->rhocorelig[idx] < 0.0)
                    sp->rhocorelig[idx] = 0.0;

                if (pct.gridpe == 0 && write_flag)
                    if (rfil < rcut)
                        fprintf (psp, "%15.8f  %15.8f\n", rfil, sp->rhocorelig[idx]);

                rfil += sp->drlig;

            }                   /* end for */

            
	    if (pct.gridpe == 0 && write_flag)
                fprintf (psp, "\n&&\n");

        }                       /* end if */


        if (pct.gridpe == 0 && write_flag)
            fclose (psp);

    }                           /* end for */


    my_free (work);

}                               /* end init_kbr.c */

/******/
