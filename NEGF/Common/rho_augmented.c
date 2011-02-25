/************************** SVN Revision Information **************************
 **    $Id: rho_augmented.c 1242 2011-02-02 18:55:23Z luw $    **
******************************************************************************/
 
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "md.h"

/* begin shuchun wang */
void rho_augmented(REAL * rho, REAL * global_mat_X)
{

    int istate, incx, idx;
    int *ivec, size, idx1, idx2;
    int nh, icount, ncount, i, j, ion;
    REAL *qnmI, *qtpr;
    REAL *product, *ptr_product;
    REAL time1, time2;
    ION *iptr;
    SPECIES *sp;
    double alfa;

    alfa = 2.0;


    size = ct.num_ions * ct.max_nl * ct.max_nl;
    my_malloc_init( product, size, REAL );
    for (idx = 0; idx < size; idx++)
        product[idx] = 0.0;

    rho_Qnm_mat(product, global_mat_X);
    for (idx = 0; idx < size; idx++)
        product[idx] *= alfa;
    global_sums(product, &size);

    for (ion = 0; ion < ct.num_ions; ion++)
    {
        iptr = &ct.ions[ion];
        sp = &ct.sp[iptr->species];
        nh = sp->num_projectors;
        ptr_product = product + ion * ct.max_nl * ct.max_nl;
        ivec = pct.Qindex[ion];
        ncount = pct.Qidxptrlen[ion];
        qnmI = pct.augfunc[ion];

        if (pct.Qidxptrlen[ion])
        {

            idx = 0;
            for (i = 0; i < nh; i++)
            {
                for (j = i; j < nh; j++)
                {
                    idx1 = i * ct.max_nl + j;
                    idx2 = j * ct.max_nl + i;
                    qtpr = qnmI + idx * ncount;
                    for (icount = 0; icount < ncount; icount++)
                    {
                        if (i != j)
                            rho[ivec[icount]] +=
                                qtpr[icount] * (ptr_product[idx1] + ptr_product[idx2]);
                        else
                            rho[ivec[icount]] += qtpr[icount] * ptr_product[idx1];
                    }           /*end for icount */
                    idx++;
                }               /*end for j */
            }                   /*end for i */

        }                       /*end if */

    }                           /*end for ion */

    /* release our memory */
    my_free(product);

}

/* end shuchun wang */
