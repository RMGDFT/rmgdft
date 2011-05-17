/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "main.h"


void partial_QI (int ion, REAL * QI_R, ION * iptr)
{
    int idx, idx1, i, j, k, l, num;
    int ix, iy, iz, ih, lmax, size, *dvec;
    int lpx[9][9], lpl[9][9][9], nhtol[18];
    int nhtom[18], indv[18], nh, icount;
    REAL x[3], cx[3], r, invdr, ap[25][9][9];
    REAL ylm[25], ylm_x[25], ylm_y[25], ylm_z[25];
    REAL xc, yc, zc;
    REAL *QI_x, *QI_y, *QI_z;
    REAL *qnmlig_tpr, *drqnmlig_tpr;
    SPECIES *sp;


    for (k = 0; k < 9; k++)
    {
        for (j = 0; j < 9; j++)
        {
            for (i = 0; i < 25; i++)
            {
                ap[i][j][k] = 0.0;
                if (i < 9)
                    lpl[i][j][k] = 0;
            }
            lpx[k][j] = 0;
        }
    }

    for (j = 0; j < 18; j++)
    {
        nhtol[j] = 0;
        nhtom[j] = 0;
        indv[j] = 0;
    }
    nh = 0;

    sp = &ct.sp[iptr->species];
    lmax = sp->llbeta[sp->nbeta - 1];

    ih = 0;
    for (i = 0; i < sp->nbeta; i++)
    {
        l = sp->llbeta[i];
        for (j = 0; j < 2 * l + 1; j++)
        {
            nhtol[ih] = l;
            nhtom[ih] = j;
            indv[ih] = i;
            ih = ih + 1;
        }
    }
    nh = ih;

    aainit (lmax + 1, 2 * lmax + 1, 2 * lmax + 1, 4 * lmax + 1, (lmax + 1) * (lmax + 1), ap, lpx,
            lpl);

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
                            qval_R (i, j, r, cx, qnmlig_tpr, drqnmlig_tpr, invdr, nhtol,
                                    nhtom, indv, ylm, ylm_x, ylm_y, ylm_z, ap, lpx, lpl,
                                    &QI_x[idx1], &QI_y[idx1], &QI_z[idx1], sp);
                            ++num;

                        }       /*end for j */

                    }           /* end for i */

                    icount++;
                }               /* end if */

                ++idx;
                zc += ct.hzzgrid;

            }                   /* end for iz */
            yc += ct.hyygrid;

        }                       /* end for iy */
        xc += ct.hxxgrid;

    }                           /* end for ix */

}
