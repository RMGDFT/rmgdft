/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

#include "portability.h"
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "grid.h"
#include "main.h"
#include "common_prototypes.h"


void init_qfunct (void)
{
    int isp, idx, i, j, k, num, il, jl, ll, it1, ibrav;
    double rfil, t1, t2, scale;
    double work[MAX_RGRID];
    double *qnmlig_tpr, *drqnmlig_tpr, *qnm_tpr, *workr;
    SPECIES *sp;
    char newname1[MAX_PATH], newname2[MAX_PATH];
    FILE *fqq = NULL;
    FILE *fdq = NULL;
    const bool SET = true;

    if(ct.norm_conserving_pp) return;
    ibrav = get_ibrav_type();

    my_malloc (workr, MAX_QLIG, double);

    scale = 1.0;
    if (ibrav == CUBIC_BC)
        scale = 1.1;
    if (ibrav == CUBIC_FC)
        scale = 1.3;
    for (isp = 0; isp < ct.num_species; isp++)
    {

        sp = &ct.sp[isp];
        if (verify_boolean ("write_pseudopotential_plots", &SET))
        {
            snprintf (newname1, MAX_PATH, "q_%s.xmgr", sp->atomic_symbol);
            snprintf (newname2, MAX_PATH, "drq_%s.xmgr", sp->atomic_symbol);
            if (pct.gridpe == 0)
            {
                my_fopen (fqq, newname1, "w+");
                my_fopen (fdq, newname2, "w+");
            }
        }


        t1 = 2.0 * scale * (double) get_FG_RATIO() *sp->qradius / ct.hmingrid;
        /*		t1=2.0 *  scale * sp->qradius / ct.hmingrid;*/
        t1 = modf (t1, &t2);
        it1 = (int) t2;
        if (t1 > 0.5)
            it1++;
        if (!(it1 % 2))
            it1++;
        sp->qdim = it1;
        /*		sp->qdim = 2 * FG_NX * (it1 / 2) + 1;*/

        //        if ((sp->qdim >= get_FNX_GRID()) || (sp->qdim >= get_FNY_GRID())
        //            || (sp->qdim >= get_FNZ_GRID()))
        //            error_handler ("nlocal potential radius exceeds global grid size");
        if (sp->qdim >= get_FNX_GRID()) sp->qdim = get_FNX_GRID();
        if (sp->qdim >= get_FNY_GRID()) sp->qdim = get_FNY_GRID();
        if (sp->qdim >= get_FNZ_GRID()) sp->qdim = get_FNZ_GRID();

        if (ct.max_Qpoints < (sp->qdim * sp->qdim * sp->qdim))
            ct.max_Qpoints = sp->qdim * sp->qdim * sp->qdim;

        num = (sp->nbeta * (sp->nbeta + 1)) * sp->nlc / 2;
        my_malloc (sp->qnmlig, num * MAX_QLIG, double);
        my_malloc (sp->drqnmlig, num * MAX_QLIG, double);
        idx = 0;

        t1 = sp->qdim / get_FG_RATIO() + 1;
        sp->drqlig = 2.0 * sqrt (THREE) * (t1 + 1.0) * ct.hmaxgrid / TWO;
        if (ibrav == HEXAGONAL)
            sp->drqlig *= 2.0;
        t1 = (double) MAX_QLIG;
        sp->drqlig = sp->drqlig / t1;

        for (i = 0; i < sp->nbeta; i++)
        {
            il = sp->llbeta[i];
            for (j = i; j < sp->nbeta; j++)
            {
                jl = sp->llbeta[j];
                idx = j * (j + 1) / 2 + i;
                qnm_tpr = sp->qnm + idx * MAX_RGRID;
                for (ll = fabs (il - jl); ll <= il + jl; ll = ll + 2)
                {
                    for (k = 0; k < MAX_RGRID; k++)
                        work[k] = 0.0;

                    for (k = 0; k < sp->rg_points; k++)
                    {
                        if (sp->r[k] >= sp->rinner[ll])
                            work[k] = qnm_tpr[k];
                        else
                            work[k] = get_QnmL (idx, ll, sp->r[k], sp);
                    }
                    qnmlig_tpr = sp->qnmlig + (idx * sp->nlc + ll) * MAX_QLIG;
                    drqnmlig_tpr = sp->drqnmlig + (idx * sp->nlc + ll) * MAX_QLIG;

                    if (pct.gridpe == 0 && verify_boolean ("write_pseudopotential_plots", &SET))
                    {
                        for (k = 0; k < sp->kkbeta; k++)
                            fprintf (fqq, "%e  %e\n", sp->r[k], work[k]);
                        fprintf (fqq, "&\n");

                    }

                    filter_potential(work, &sp->r[0], sp->rg_points, sp->nlradius, 0, ct.qcparm, qnmlig_tpr, 
                            &sp->rab[0], ll, sp->drqlig, sp->gwidth, MAX_QLIG, sp->nlrcut[sp->llbeta[i]], sp->rwidth, drqnmlig_tpr);

                    /*Is this necessary ???*/
                    if (ll)
                        qnmlig_tpr[0] = 0.0;


                    /*Write final filtered Q function if requested*/
                    if (pct.gridpe == 0 && verify_boolean ("write_pseudopotential_plots", &SET))
                    {

                        rfil = ZERO;
                        for (k = 0; k < MAX_QLIG; k++)
                        {
                            fprintf (fqq, "%e  %e\n", rfil, qnmlig_tpr[k]);
                            fprintf (fdq, "%e  %e\n", rfil, drqnmlig_tpr[k]);

                            rfil += sp->drqlig;
                        }

                        fprintf (fqq, "&\n");
                        fprintf (fdq, "&\n");
                    }
                }               /*end for ll */
            }                   /*end for j */
        }                       /*end for i */

        if (pct.gridpe == 0 && verify_boolean ("write_pseudopotential_plots", &SET))
        {
            fclose (fqq);
            fclose (fdq);
        }

    }                           /*end for isp */

    my_free (workr);

}
