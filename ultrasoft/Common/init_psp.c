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


void init_psp (void)
{

    int isp, idx, ip, it1;
    SPECIES *sp;
    REAL *work, *workr, Zv, rc, rfil;
    REAL t1, t2, rcut, scale;
    char name[] = "projectors";
    char newname[20];
    FILE *psp = NULL;

    if (pct.gridpe == 0 && verify ("do_write_pseudopotential_plots", &SET))
    {

        printf ("QMD status message: Opening projectors.xmgr\n");
        fflush (NULL);

    }                           /* end if */


    my_malloc (work, MAX_RGRID + MAX_LOCAL_LIG, REAL);
    workr = work + MAX_RGRID;

    /*Initialize max_nlpoints and max_nlfpoints */
    ct.max_nlpoints = 0;
    ct.max_nlfpoints = 0;

    /* Set the scaling factor for determining the radius of the local grids */
    scale = 1.0;
    if (ct.ibrav == CUBIC_BC)
        scale = 1.1;
    if (ct.ibrav == CUBIC_FC)
        scale = 1.3;


    /* Loop over species */
    for (isp = 0; isp < ct.num_species; isp++)
    {
        sprintf (newname, "%s%d.xmgr", name, isp);
        if (pct.gridpe == 0 && verify ("do_write_pseudopotential_plots", &SET))
            my_fopen (psp, newname, "w+");
        sp = &ct.sp[isp];

        t1 = 2.0 * scale * (REAL) FG_NX *sp->lradius / ct.hmingrid;
        t1 = modf (t1, &t2);
        it1 = (int) t2;
        if (t1 > 0.5)
            it1++;
        if (!(it1 % 2))
            it1++;
        sp->ldim = it1;

        if ((sp->ldim >= ct.psi_fnxgrid) || (sp->ldim >= ct.psi_fnygrid)
            || (sp->ldim >= ct.psi_fnzgrid))
            error_handler ("local potential radius exceeds global grid size");

        t1 = 2.0 * scale * sp->nlradius / ct.hmingrid;
        t1 = modf (t1, &t2);
        it1 = (int) t2;
        if (t1 > 0.5)
            it1++;
        if (!(it1 % 2))
            it1++;
        sp->nldim = it1;
        sp->nlfdim = ct.nxfgrid * it1;

        if ((sp->nldim >= ct.psi_nxgrid) || (sp->nldim >= ct.psi_nygrid)
            || (sp->nldim >= ct.psi_nzgrid))
            error_handler ("local potential radius exceeds global grid size");

        /*ct.max_nlpoints is max of nldim*nldim*nldim for all species */
        if (ct.max_nlpoints < (sp->nldim * sp->nldim * sp->nldim))
            ct.max_nlpoints = sp->nldim * sp->nldim * sp->nldim;

        if (ct.max_nlfpoints < (sp->nlfdim * sp->nlfdim * sp->nlfdim))
            ct.max_nlfpoints = sp->nlfdim * sp->nlfdim * sp->nlfdim;

        /* Check value to make sure a local potential is defined */
        sp->localidx = -1;

        Zv = sp->zvalence;
        rc = sp->rc;
        sp->localidx = sp->local;       /*thinking */


        /* Loop over radial grid points */
        for (idx = 0; idx < sp->rg_points; idx++)
        {

            /* Generate the difference potential */
            work[idx] = sp->vloc0[idx] + Zv * erf (sp->r[idx] / rc) / sp->r[idx];

        }                       /* end for */

        if (pct.gridpe == 0 && verify ("do_write_pseudopotential_plots", &SET))
        {
            for (idx = 0; idx < sp->rg_points; idx++)
            {
                if (sp->r[idx] < sp->lradius)
                    fprintf (psp, "%15.8f  %15.8f\n", sp->r[idx], work[idx]);
            }
            fprintf (psp, "\n&&\n");
        }



        /*sp->drlig = sqrt(3.0) * (sp->ldim + 1.0) * ct.hmaxgrid / 2.0 /(REAL)FG_NX; */
        t1 = sp->ldim / FG_NX + 1;
        sp->drlig = sqrt (3.0) * (t1 + 1.0) * ct.hmaxgrid / 2.0;
        if (ct.ibrav == HEXAGONAL)
            sp->drlig *= 2.0;
        t1 = (REAL) MAX_LOCAL_LIG;
        sp->drlig = sp->drlig / t1;


        /* Transform to g-space and filter it */
        rft1 (ct.cparm, work, &sp->r[0], sp->localig, &sp->rab[0], sp->rg_points, 0, sp->drlig,
              sp->gwidth, MAX_LOCAL_LIG);


        /* Evaluate it's radial derivative */
        for (idx = 0; idx < MAX_LOCAL_LIG; idx++)
        {

            workr[idx] = sp->drlig * ((REAL) idx);

        }                       /* end for */
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
                sp->localig[idx] = sp->localig[idx] * exp (-sp->rwidth * t1 * t1);
                sp->drlocalig[idx] = sp->drlocalig[idx] * exp (-sp->rwidth * t1 * t1);

                if (fabs (sp->localig[idx]) < 1.e-35)
                    sp->localig[idx] = 0.;
                if (fabs (sp->drlocalig[idx]) < 1.e-35)
                    sp->drlocalig[idx] = 0.;

            }                   /* end if */

            /* output local projector */
            if (pct.gridpe == 0 && verify ("do_write_pseudopotential_plots", &SET))
            {

		if (rfil < sp->lradius)
                    fprintf (psp, "%15.8f  %15.8f  %15.8f\n", rfil, sp->localig[idx], sp->drlocalig[idx]);

            }                   /* endif */

            rfil += sp->drlig;

        }                       /* end for idx */

        /* output xmgr data separator */
        if (pct.gridpe == 0 && verify ("do_write_pseudopotential_plots", &SET))
        {

            fprintf (psp, "\n&&\n");

        }                       /* endif */

        sp->drnlig = sqrt (3.0) * (sp->nldim + 1.0) * ct.hmaxgrid / 2.0;
        if (ct.ibrav == HEXAGONAL)
            sp->drnlig *= 2.0;
        t1 = (REAL) MAX_LOCAL_LIG;
        sp->drnlig = sp->drnlig / t1;

        for (ip = 0; ip < sp->nbeta; ip++)
        {

            if (pct.gridpe == 0 && verify ("do_write_pseudopotential_plots", &SET))
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
            {

                workr[idx] = sp->drnlig * ((REAL) idx);

            }                   /* end for */
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
                    sp->betalig[ip][idx] = sp->betalig[ip][idx] * exp (-sp->rwidth * t1 * t1);
                    sp->drbetalig[ip][idx] = sp->drbetalig[ip][idx] * exp (-sp->rwidth * t1 * t1);

                    if (fabs (sp->betalig[ip][idx]) < 1.e-35)
                        sp->betalig[ip][idx] = 0.;
                    if (fabs (sp->drbetalig[ip][idx]) < 1.e-35)
                        sp->drbetalig[ip][idx] = 0.;


                }               /* end if */

                /* output non-local projector */
                if (pct.gridpe == 0 && verify ("do_write_pseudopotential_plots", &SET))
                {

/*                        if(!(idx % 32))*/
                    fprintf (psp, "%15.8f  %15.8f %15.8f\n", rfil, sp->betalig[ip][idx], sp->drbetalig[ip][idx]);

                }               /* endif */

                rfil += sp->drnlig;

            }                   /* end for */

            /* output xmgr data separator */
            if (pct.gridpe == 0 && verify ("do_write_pseudopotential_plots", &SET))
            {

                fprintf (psp, "\n&&\n");

            }                   /* endif */

        }                       /* end for ip */

        /* Now take care of the core charge if nlcc is being used */
        if (sp->nlccflag)
        {

            for (idx = 0; idx < sp->rg_points; idx++)
            {

                work[idx] = sp->rspsco[idx] / (4.0 * PI);

            }                   /* end for */

            if (pct.gridpe == 0 && verify ("do_write_pseudopotential_plots", &SET))
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
                    sp->rhocorelig[idx] = sp->rhocorelig[idx] * exp (-sp->rwidth * t1 * t1);

                    if (fabs (sp->rhocorelig[idx]) < 1.e-35)
                        sp->rhocorelig[idx] = 0.;

                }               /* end if */

                if (sp->rhocorelig[idx] < 0.0)
                    sp->rhocorelig[idx] = 0.0;

                if (pct.gridpe == 0 && verify ("do_write_pseudopotential_plots", &SET))
                {

                    if (rfil < rcut)
                        fprintf (psp, "%15.8f  %15.8f\n", rfil, sp->rhocorelig[idx]);

                }               /* endif */

                rfil += sp->drlig;

            }                   /* end for */

            if (pct.gridpe == 0 && verify ("do_write_pseudopotential_plots", &SET))
                fprintf (psp, "\n&&\n");


        }                       /* end if */



        /* Make sure that a local potential was specified */
        if (sp->localidx < 0)
            error_handler ("No local potential defined");

        if (pct.gridpe == 0 && verify ("do_write_pseudopotential_plots", &SET))
            fclose (psp);

    }                           /* end for */


    my_free (work);

    if (pct.gridpe == 0 && verify ("do_write_pseudopotential_plots", &SET))
    {

        printf ("QMD status message: Closing projectors.xmgr\n");
        fflush (NULL);

    }                           /* end if */


}                               /* end init_kbr.c */

/******/
