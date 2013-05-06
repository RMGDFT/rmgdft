/************************** SVN Revision Information **************************
 **    $Id: 
******************************************************************************/
 


#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <assert.h>
#include "main.h"
#include "pmo.h"





void H01_to_HCL (double *H01_global, double *HCL_local, int iprobe)
{
    int i, j, ii, jj, iii, jjj, li, lj, maxli;
    int iistart, jjstart, limb, ljnb;
    int mycol, myrow, nprow, npcol;
    int ictxt, mb, nb, *desca, n0;
    int idx_local, idx_global, item;

    n0 = cei.probe_in_block[iprobe - 1];  /* n0 => block index */
    desca = &pmo.desc_cond_lead[ (n0 + (iprobe - 1) * ct.num_blocks) * DLEN ]; /* (C,L) */

    ictxt = desca[1], mb = desca[4], nb = desca[5];

    Cblacs_gridinfo (ictxt, &nprow, &npcol, &myrow, &mycol);

    /* li, lj: the index in distributed matrix
     * i, j: related index in the global matrix
     */


    for (lj = 0; lj < pmo.mxlocc_lead[iprobe-1]; lj++)
    {
        for (li = 0; li < pmo.mxllda_cond[n0]; li++)
        {
            i = li/mb * nprow *mb + myrow * mb + li - li/mb * mb;
            j = lj/nb * npcol *nb + mycol * nb + lj - lj/nb * nb;


            idx_local = lj * pmo.mxllda_cond[n0] + li;
            HCL_local[idx_local] = 0.0;


            /* for probe 1,3,5,... copy H01 to the begining 
             */

            if( (iprobe != iprobe/2 *2) && i < lcr[iprobe].num_states)
            {
                idx_global = j * lcr[iprobe].num_states + i;
                HCL_local[idx_local] = H01_global [idx_global];
            }

            /* for probe 2,4,6,... copy H01 to the end 
             */

            item = ct.block_dim[n0]- lcr[iprobe].num_states;
            if( (iprobe == iprobe/2 *2) && i >= item)
            {
                idx_global = j * lcr[iprobe].num_states + i - item;
                HCL_local[idx_local] = H01_global [idx_global];
            }

        }

    }
}


