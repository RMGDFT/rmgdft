/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "md.h"


void init_qfunct(void)
{
    int isp, idx, i, j, k, num, il, jl, ll, it1;
    REAL rcut, rfil, t1, t2, scale;
    REAL work[MAX_RGRID];
    REAL *qnmlig_tpr, *drqnmlig_tpr, *qnm_tpr, *workr;
    SPECIES *sp;
    char name1[] = "qfunction", name2[] = "drqfunction";
    char newname1[20], newname2[20];
    FILE *fqq = NULL;
    FILE *fdq = NULL;

    my_malloc_init( workr, MAX_QLIG, REAL );

    scale = 1.0;
    if (ct.ibrav == CUBIC_BC)
        scale = 1.1;
    if (ct.ibrav == CUBIC_FC)
        scale = 1.3;
    for (isp = 0; isp < ct.num_species; isp++)
    {

        sprintf(newname1, "%s%d.xmgr", name1, isp);
        sprintf(newname2, "%s%d.xmgr", name2, isp);
        if (pct.gridpe == 0)
        {
            fqq = fopen(newname1, "w+");
            fdq = fopen(newname2, "w+");
        }

        sp = &ct.sp[isp];

        t1 = 2.0 * scale * (REAL) RHO_NX *sp->qradius / ct.hmingrid;
/*		t1=2.0 *  scale * sp->qradius / ct.hmingrid;*/
        t1 = modf(t1, &t2);
        it1 = (int) t2;
        if (t1 > 0.5)
            it1++;
        if (!(it1 % 2))
            it1++;
        sp->qdim = it1;
/*		sp->qdim = 2 * RHO_NX * (it1 / 2) + 1;*/

        if ((sp->qdim >= ct.vh_nxgrid) || (sp->qdim >= ct.vh_nygrid) || (sp->qdim >= ct.vh_nzgrid))
            error_handler("nlocal potential radius exceeds global grid size");
        if (ct.max_Qpoints < (sp->qdim * sp->qdim * sp->qdim))
            ct.max_Qpoints = sp->qdim * sp->qdim * sp->qdim;

        num = (sp->nbeta * (sp->nbeta + 1)) * sp->nlc / 2;
        my_malloc_init( sp->qnmlig, num * MAX_QLIG, REAL );
        my_malloc_init( sp->drqnmlig, num * MAX_QLIG, REAL );
        idx = 0;

        t1 = sp->qdim / RHO_NX + 1;
        sp->drqlig = sqrt(THREE) * (t1 + 1.0) * ct.hmaxgrid / TWO;
        if (ct.ibrav == HEXAGONAL)
            sp->drqlig *= 2.0;
        t1 = (REAL) MAX_QLIG;
        sp->drqlig = sp->drqlig / t1;

        for (i = 0; i < sp->nbeta; i++)
        {
            il = sp->llbeta[i];
            for (j = i; j < sp->nbeta; j++)
            {
                jl = sp->llbeta[j];
                idx = j * (j + 1) / 2 + i;
                qnm_tpr = sp->qnm + idx * MAX_RGRID;
                for (ll = fabs(il - jl); ll <= il + jl; ll = ll + 2)
                {
                    for (k = 0; k < MAX_RGRID; k++)
                        work[k] = 0.0;

                    for (k = 0; k < MAX_RGRID; k++)
                    {
                        if (sp->r[k] >= sp->rinner[ll])
                            work[k] = qnm_tpr[k];
                        else
                            work[k] = get_QnmL(idx, ll, sp->r[k], sp);
                    }
                    qnmlig_tpr = sp->qnmlig + (idx * sp->nlc + ll) * MAX_QLIG;
                    drqnmlig_tpr = sp->drqnmlig + (idx * sp->nlc + ll) * MAX_QLIG;

                    if (pct.gridpe == 0)
                    {
                        for (k = 0; k < sp->kkbeta; k++)
                            fprintf(fqq, "%e  %e\n", sp->r[k], work[k]);
                        fprintf(fqq, "&&\n");

                    }

                    rft1(ct.qcparm, &work[0], &sp->r[0], qnmlig_tpr,
                         &sp->rab[0], sp->rg_points, ll, sp->drqlig, sp->gwidth, MAX_QLIG);


                    for (k = 0; k < MAX_QLIG; k++)
                    {
                        workr[k] = sp->drqlig * ((REAL) k);
                    }

                    radiff(qnmlig_tpr, drqnmlig_tpr, workr, MAX_QLIG, 0.0);

                    qnmlig_tpr[0] = TWO * qnmlig_tpr[1] - qnmlig_tpr[2];
                    if (ll)
                        qnmlig_tpr[0] = 0.0;
                    drqnmlig_tpr[1] = TWO * drqnmlig_tpr[2] - drqnmlig_tpr[3];
                    drqnmlig_tpr[0] = TWO * drqnmlig_tpr[1] - drqnmlig_tpr[2];

                    rcut = sp->nlrcut[sp->llbeta[i]];
                    rfil = ZERO;
                    for (k = 0; k < MAX_QLIG; k++)
                    {
                        if (rfil > rcut)
                        {
                            t1 = (rfil - rcut) / rcut;
                            qnmlig_tpr[k] = qnmlig_tpr[k] * exp(-sp->rwidth * t1 * t1);
                            drqnmlig_tpr[k] = drqnmlig_tpr[k] * exp(-sp->rwidth * t1 * t1);
                            if (fabs(qnmlig_tpr[k]) < 1.e-35)
                                qnmlig_tpr[k] = 0.0;
                            if (fabs(drqnmlig_tpr[k]) < 1.e-35)
                                drqnmlig_tpr[k] = 0.0;
                        }       /*end for if */

                        if (pct.gridpe == 0)
                        {
                            fprintf(fqq, "%e  %e\n", rfil, qnmlig_tpr[k]);
                            fprintf(fdq, "%e  %e\n", rfil, drqnmlig_tpr[k]);
                        }

                        rfil += sp->drqlig;
                    }           /*end for k */

                    if (pct.gridpe == 0)
                    {
                        fprintf (fqq, "&&\n");
                        fprintf (fdq, "&&\n");
                    }
                }               /*end for ll */
            }                   /*end for j */
        }                       /*end for i */

        if (pct.gridpe == 0)
        {
            fclose(fqq);
            fclose(fdq);
        }

    }                           /*end for isp */

    my_free(workr);

}
