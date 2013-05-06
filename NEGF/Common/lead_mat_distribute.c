/************************** SVN Revision Information **************************
 **    $Id: 
******************************************************************************/
 


#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <assert.h>
#include "main.h"
#include "pmo.h"





/*
*
*      distributes matrice a
*
*/
void lead_mat_distribute (double *a_local, int *desca, double *a_global,
int iprobe)
{
    int i, j, ii, jj, iii, jjj, li, lj, maxli;
    int iistart, jjstart, limb, ljnb;
    int mycol, myrow, nprow, npcol;
    int ictxt = desca[1], mb = desca[4], nb = desca[5], mxllda = desca[8];
    int idx_local, idx_global;



    Cblacs_gridinfo (ictxt, &nprow, &npcol, &myrow, &mycol);

    /* li, lj: the index in distributed matrix
     * i, j: related index in the global matrix
     */


    for (lj = 0; lj < pmo.mxlocc_lead[iprobe-1]; lj++)
    {
        for (li = 0; li < pmo.mxllda_lead[iprobe-1]; li++)
        {
            i = li/mb * nprow *mb + myrow * mb + li - li/mb * mb;
            j = lj/nb * npcol *nb + mycol * nb + lj - lj/nb * nb;

            idx_local = lj * pmo.mxllda_lead[iprobe-1] + li;
            idx_global = j * lcr[iprobe].num_states + i;
            a_local[idx_local] = a_global [idx_global];

        }

    }
}
