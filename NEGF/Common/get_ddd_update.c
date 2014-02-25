/************************** SVN Revision Information **************************
 **    $Id$    **
 ******************************************************************************/

#include <float.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "main.h"
#include "init_var.h"
#include "LCR.h"

void get_ddd_update (rmg_double_t * veff)
{
    int idx, i, j, ion;
    int nh, ncount, icount;
    rmg_double_t *qnmI, *dnmI, *sum;
    int *ivec, sum_dim, sum_idx;
    ION *iptr;
    SPECIES *sp;


    /*Count the number of elements in sum array */
    sum_dim = 0;
    for (ion = 0; ion < ct.num_ions; ion++)
    {
        iptr = &ct.ions[ion];
        sp = &ct.sp[iptr->species];
        
        nh = sp->nh;

        /*Number of elements is sum of 1+2+3+...+nh */
        sum_dim += nh * (nh + 1) / 2;
    }

    my_calloc (sum, sum_dim, rmg_double_t);


    sum_idx = 0;

    for (ion = 0; ion < ct.num_ions; ion++)
    {
        iptr = &ct.ions[ion];
        sp = &ct.sp[iptr->species];

        ivec = pct.Qindex[ion];
        nh = sp->nh;
        ncount = pct.Qidxptrlen[ion];

        if (pct.dnmI[ion] == NULL)
            my_malloc (pct.dnmI[ion], nh * nh, rmg_double_t);

        idx = 0;
        for (i = 0; i < nh; i++)
        {
            for (j = i; j < nh; j++)
            {
                if (ncount)
                {
                    qnmI = pct.augfunc[ion] + idx * ncount;
                    for (icount = 0; icount < ncount; icount++)
                    {
                        sum[sum_idx] += qnmI[icount] * veff[ivec[icount]];
                    }
                }               /*end if (ncount) */

                sum[sum_idx] *= ct.vel_f;

                idx++;
                sum_idx++;
            }                   /*end for (j = i; j < nh; j++) */
        }                       /*end for (i = 0; i < nh; i++) */

    }                           /*end for (ion = 0; ion < ct.num_ions; ion++) */



    if (sum_idx != sum_dim)
        error_handler ("Problem with sum index");

    global_sums (sum, &sum_dim, pct.grid_comm);

    sum_idx = 0;

    for (ion = 0; ion < ct.num_ions; ion++)
    {
        iptr = &ct.ions[ion];
        sp = &ct.sp[iptr->species];

        nh = sp->nh;

        dnmI = pct.dnmI[ion];

        for (i = 0; i < nh; i++)
        {
            for (j = i; j < nh; j++)
            {

                /*if (fabs (sum[sum_idx]) < 1.0e-10)
                   sum[sum_idx] = 0.0; */

                dnmI[i * nh + j] = sum[sum_idx];
                if (i != j)
                    dnmI[j * nh + i] = dnmI[i * nh + j];

                sum_idx++;

            }                   /*end for j */
        }                       /*end for i */
    }                           /*end for ion */

    my_free (sum);
}
