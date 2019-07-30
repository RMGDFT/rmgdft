/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "main.h"
#include "prototypes_on.h"

/* begin shuchun wang */
void rho_augmented(double * rho, double * global_mat_X,
int *state_begin, int *state_end, int *num_nonlocal_ion, double *kbpsi,
int max_ion_nonlocal, double *kbpsi_comm, int *ionidx_allproc)

{

    int idx;
    int *ivec, size, idx1, idx2;
    int nh, icount, ncount, i, j, ion;
    double *qnmI, *qtpr;
    double *product, *ptr_product;
    ION *iptr;
    SPECIES *sp;



    size = ct.num_ions * ct.max_nl * ct.max_nl;
    my_malloc_init( product, size, double );
    for (idx = 0; idx < size; idx++)
        product[idx] = 0.0;

    

    RmgTimer *RT5 = new RmgTimer("3-get_new_rho: augmented_qnm");

    rho_Qnm_mat(product, global_mat_X, state_begin, state_end, num_nonlocal_ion, 
            kbpsi, max_ion_nonlocal, kbpsi_comm, ionidx_allproc);
    delete [] RT5;

    global_sums(product, &size, pct.grid_comm);

    RmgTimer *RT6 = new RmgTimer("3-get_new_rho: augmented_Q(r)");
    for (ion = 0; ion < ct.num_ions; ion++)
    {
        iptr = &Atoms[ion];
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

    delete [] RT6;
    /* release our memory */
    my_free(product);

}

/* end shuchun wang */
