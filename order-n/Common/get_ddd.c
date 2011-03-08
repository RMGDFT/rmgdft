/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
#include <float.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "md.h"

void get_ddd(REAL * veff)
{
    int idx, i, j, ion;
    int nh, ncount, icount;
    REAL *qnmI, *dnmI, sum;
    int *ivec;
    ION *iptr;
    SPECIES *sp;
    int size;

    for (ion = 0; ion < ct.num_ions; ion++)
    {
        iptr = &ct.ions[ion];
        sp = &ct.sp[iptr->species];

        ivec = pct.Qindex[ion];
        nh = sp->num_projectors;
        ncount = pct.Qidxptrlen[ion];

        if (pct.dnmI[ion] == NULL)
            my_calloc( pct.dnmI[ion], nh * nh, REAL );
        dnmI = pct.dnmI[ion];

        idx = 0;
        sum = 0.0;
        for (i = 0; i < nh; i++)
        {
            for (j = i; j < nh; j++)
            {
                sum = 0.0;
                if (ncount)
                {
                    qnmI = pct.augfunc[ion] + idx * ncount;
                    for (icount = 0; icount < ncount; icount++)
                    {
                        sum += qnmI[icount] * veff[ivec[icount]];
                    }
                }
                sum = real_sum_all(sum);
                sum = sum * ct.vel_f;

                if (fabs(sum) < 1.0e-10)
                    sum = 0.0;
                dnmI[i * nh + j] = sp->ddd0[i][j] + sum;
                if (i != j)
                    dnmI[j * nh + i] = dnmI[i * nh + j];
                idx++;
            }                   /*end for j */
        }                       /*end for i */

    }                           /*end for ion */
}
