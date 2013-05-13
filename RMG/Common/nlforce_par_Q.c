/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "main.h"

void nlforce_par_Q (rmg_double_t * veff, rmg_double_t * gamma, int ion, ION * iptr, int nh, rmg_double_t * forces)
{
    int idx1, idx2, n, m, count, icount, size;
    int *pidx;
    rmg_double_t tmp[3];
    rmg_double_t *QnmI_R, *QnmI_x, *QnmI_y, *QnmI_z;

    /*Forces array is assumed to be already initialized */
    /*for(idx1=0;idx1<3;idx1++)  forces[idx1]=0.0; */

    count = pct.Qidxptrlen[ion];
    pidx = pct.Qindex[ion];

    if (count)
    {
        size = (nh * (nh + 1) / 2) * count;
        my_malloc (QnmI_R, 3 * size, rmg_double_t);
        QnmI_x = QnmI_R;
        QnmI_y = QnmI_x + size;
        QnmI_z = QnmI_y + size;

        for (idx1 = 0; idx1 < 3 * size; idx1++)
            QnmI_R[idx1] = 0.0;

        partial_QI (ion, QnmI_R, iptr);

        for (icount = 0; icount < count; icount++)
        {
            tmp[0] = 0.0;
            tmp[1] = 0.0;
            tmp[2] = 0.0;

            idx2 = 0;
            for (n = 0; n < nh; n++)
            {
                for (m = n; m < nh; m++)
                {
                    idx1 = idx2 * count + icount;
                    if (m == n)
                    {
                        tmp[0] += QnmI_x[idx1] * gamma[idx2];
                        tmp[1] += QnmI_y[idx1] * gamma[idx2];
                        tmp[2] += QnmI_z[idx1] * gamma[idx2];
                    }
                    else
                    {
                        tmp[0] += 2.0 * QnmI_x[idx1] * gamma[idx2];
                        tmp[1] += 2.0 * QnmI_y[idx1] * gamma[idx2];
                        tmp[2] += 2.0 * QnmI_z[idx1] * gamma[idx2];
                    }

                    ++idx2;
                }
            }
            forces[0] += veff[pidx[icount]] * tmp[0];
            forces[1] += veff[pidx[icount]] * tmp[1];
            forces[2] += veff[pidx[icount]] * tmp[2];
        }

        my_free (QnmI_R);
    }



}
