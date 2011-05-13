/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
#include <float.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "md.h"

void get_partial_ddd (REAL * QnmI_R, REAL * veff, int ion, int nh)
{
    int idx, i, j;
    int ncount, icount, size;
    REAL *dnmI_x, *dnmI_y, *dnmI_z;
    REAL *QnmI_x, *QnmI_y, *QnmI_z;
    int *ivec;

    ncount = pct.Qidxptrlen[ion];
    ivec = pct.Qindex[ion];

    size = (nh * (nh + 1) / 2) * ncount;

    QnmI_x = QnmI_R;
    QnmI_y = QnmI_x + size;
    QnmI_z = QnmI_y + size;

    if (pct.dnmI_x[ion] == NULL)
        my_calloc( pct.dnmI_x[ion], nh * nh, REAL );
    dnmI_x = pct.dnmI_x[ion];

    if (pct.dnmI_y[ion] == NULL)
        my_calloc( pct.dnmI_y[ion], nh * nh, REAL );
    dnmI_y = pct.dnmI_y[ion];

    if (pct.dnmI_z[ion] == NULL)
        my_calloc( pct.dnmI_z[ion], nh * nh, REAL );
    dnmI_z = pct.dnmI_z[ion];

    for (idx = 0; idx < nh * nh; idx++) 
    {
        dnmI_x[idx] = 0.0;
        dnmI_y[idx] = 0.0;
        dnmI_z[idx] = 0.0;
    }

    for (i = 0; i < nh; i++)
    {
        for (j = i; j < nh; j++)
        {
            if (ncount)
            {
                for (icount = 0; icount < ncount; icount++)
                {
                    dnmI_x[i * nh + j] += QnmI_x[icount] * veff[ivec[icount]];
                    dnmI_y[i * nh + j] += QnmI_y[icount] * veff[ivec[icount]];
                    dnmI_z[i * nh + j] += QnmI_z[icount] * veff[ivec[icount]];
                }

                if (i != j)
                {
                    dnmI_x[j * nh + i] = dnmI_x[i * nh + j];
                    dnmI_y[j * nh + i] = dnmI_y[i * nh + j];
                    dnmI_z[j * nh + i] = dnmI_z[i * nh + j];
                }

                QnmI_x += ncount;
                QnmI_y += ncount;
                QnmI_z += ncount;
            }

        }
    }

    size = nh * nh;
    global_sums (dnmI_x, &size, pct.grid_comm);
    global_sums (dnmI_y, &size, pct.grid_comm);
    global_sums (dnmI_z, &size, pct.grid_comm);

    for (idx = 0; idx < size; idx++) 
    {
        dnmI_x[idx] *= ct.vel_f;
        dnmI_y[idx] *= ct.vel_f;
        dnmI_z[idx] *= ct.vel_f;
    }
}
