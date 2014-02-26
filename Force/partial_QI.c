/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "main.h"


void partial_QI (int ion, rmg_double_t * QI_R, ION * iptr)
{
    int idx, idx1, i, j, k, l, num;
    int ix, iy, iz, ih, size, *dvec;
    int lpx[9][9], lpl[9][9][9];
    int  nh, icount;
    rmg_double_t x[3], cx[3], r, invdr, ap[25][9][9];
    rmg_double_t ylm[25], ylm_x[25], ylm_y[25], ylm_z[25];
    rmg_double_t xc, yc, zc;
    rmg_double_t *QI_x, *QI_y, *QI_z;
    rmg_double_t *qnmlig_tpr, *drqnmlig_tpr;
    rmg_double_t hxxgrid, hyygrid, hzzgrid;
    SPECIES *sp;

    hxxgrid = get_hxxgrid();
    hyygrid = get_hyygrid();
    hzzgrid = get_hzzgrid();

    aainit (ct.max_l + 1, 2 * ct.max_l + 1, 2 * ct.max_l + 1, 4 * ct.max_l + 1, (ct.max_l + 1) * (ct.max_l + 1), ap, lpx,
            lpl);
    
    sp = &ct.sp[iptr->species];
    nh = sp->nh;

    size = nh * (nh + 1) / 2;
    QI_x = QI_R;
    QI_y = QI_x + size * pct.Qidxptrlen[ion];
    QI_z = QI_y + size * pct.Qidxptrlen[ion];

    idx = 0;
    icount = 0;
    dvec = pct.Qdvec[ion];

    invdr = ONE / sp->drqlig;
    qnmlig_tpr = sp->qnmlig;
    drqnmlig_tpr = sp->drqnmlig;

    /* Generate index arrays */
    xc = iptr->Qxcstart;
    for (ix = 0; ix < sp->qdim; ix++)
    {

        yc = iptr->Qycstart;
        for (iy = 0; iy < sp->qdim; iy++)
        {

            zc = iptr->Qzcstart;
            for (iz = 0; iz < sp->qdim; iz++)
            {

                if (dvec[idx])
                {
                    x[0] = xc - iptr->xtal[0];
                    x[1] = yc - iptr->xtal[1];
                    x[2] = zc - iptr->xtal[2];

                    r = metric (x);
                    to_cartesian (x, cx);
                    ylmr2 (cx, ylm);
                    ylmr2_x (cx, ylm_x);
                    ylmr2_y (cx, ylm_y);
                    ylmr2_z (cx, ylm_z);

                    num = 0;
                    for (i = 0; i < nh; i++)
                    {

                        for (j = i; j < nh; j++)
                        {

                            idx1 = num * pct.Qidxptrlen[ion] + icount;
                            qval_R (i, j, r, cx, qnmlig_tpr, drqnmlig_tpr, invdr, sp->nhtol,
                                    sp->nhtom, sp->indv, ylm, ylm_x, ylm_y, ylm_z, ap, lpx, lpl,
                                    &QI_x[idx1], &QI_y[idx1], &QI_z[idx1], sp);
                            ++num;

                        }       /*end for j */

                    }           /* end for i */

                    icount++;
                }               /* end if */

                ++idx;
                zc += hzzgrid;

            }                   /* end for iz */
            yc += hyygrid;

        }                       /* end for iy */
        xc += hxxgrid;

    }                           /* end for ix */

}
