/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "main.h"
#include "init_var.h"
#include "LCR.h"
#include "pmo.h"


void row_to_tri_p (rmg_double_t * A_tri, rmg_double_t * Aii, int N, int *ni)
{
    /* Semi_tridiagonal matrix  
     *
     *    A00   A01   0    0   ...		0
     *    A10   A11   A12  0  ...		0	
     *    0     A21   A22  A23 ....		0
     *    .     .     .    .		.	
     *    .     .     .    .		.
     *    .     .     .    .     An-1,n-1	An-1,n
     *    0     0     0          An,n-1     Ann
     *
     *   A_tri: output store the input matrix in the order of A00, A01, A11, A12, .... Ann
     *   each block is distributed in a pmo.nrow * pmo.ncol processor grid.
     *
     *   Ai,i+1  = transpose (Ai+1, i)
     *
     *   Aii: store the whole matrix in row-distribution
     *   N:   number of blocks
     *   ni:  dimension of each block
     *   for example Aii is a ni[i] x ni[i] matrix
     *               Ai,i+1 is a ni[i] x ni[i+1] matrix 
     *
     */

    int i, j,  k;
    int ndim;
    int ictxt, mb, nprow, npcol, myrow, mycol;
    int istart, jj, kk, jjj, kkk;


    ictxt = pmo.ictxt[pmo.myblacs];
    mb = pmo.mblock;

    double sym_fold = 0.5;
    /* dimension of matrix Aii */
    ndim = 0;
    for (i = 0; i < N; i++)
    {
        ndim += ni[i];
    }

    Cblacs_gridinfo(ictxt, &nprow, &npcol, &myrow, &mycol);

    for(i = 0; i < pmo.ntot; i++) A_tri[i] = 0.0;

    /* for diagonal blocks */

    istart = 0;
    for(i = 0; i < N; i++)
    {

        for(j =0; j < pmo.mxllda_cond[i]; j++)
        {
            for(k=0; k < pmo.mxlocc_cond[i]; k++)
            {

                /* for each block, distributed local index (j,k) is (jj, kk) in block matrix
                 * and (jjj,kkk) in whole matrix
                 */

                jj = (j/mb ) * nprow * mb + myrow * mb + j - j/mb * mb; 
                kk = (k/mb ) * npcol * mb + mycol * mb + k - k/mb * mb; 

                jjj = jj + istart;
                kkk = kk + istart;

                if(kkk >= ct.state_begin && kkk < ct.state_end)
                {

                    A_tri[ pmo.diag_begin[i] + k * pmo.mxllda_cond[i] + j] 
                        += Aii[(kkk-ct.state_begin) * ndim + jjj]*sym_fold;
                }
                if(jjj >= ct.state_begin && jjj < ct.state_end)
                {
                    A_tri[ pmo.diag_begin[i] + k * pmo.mxllda_cond[i] + j] 
                        += Aii[(jjj-ct.state_begin) * ndim + kkk]*sym_fold;
                }

            }
        }

        istart += ct.block_dim[i];
    }



    /* for off-diagonal blocks */

    istart = 0;
    for(i = 1; i < N; i++)
    {

        for(j =0; j < pmo.mxllda_cond[i-1]; j++)
        {
            for(k=0; k < pmo.mxlocc_cond[i]; k++)
            {

                /* for each block, distributed local index (j,k) is (jj, kk) in block matrix
                 * and (jjj,kkk) in whole matrix
                 */

                jj = (j/mb ) * nprow * mb + myrow * mb + j - j/mb * mb; 
                kk = (k/mb ) * npcol * mb + mycol * mb + k - k/mb * mb; 

                jjj = jj + istart;
                kkk = kk + istart + ct.block_dim[i-1];

                if(kkk >= ct.state_begin && kkk < ct.state_end)
                {
                    A_tri[ pmo.offdiag_begin[i-1] + k * pmo.mxllda_cond[i-1] + j] 
                        += Aii[(kkk-ct.state_begin) * ndim + jjj] * sym_fold;
                }
                if(jjj >= ct.state_begin && jjj < ct.state_end)
                {
                   A_tri[ pmo.offdiag_begin[i-1] + k * pmo.mxllda_cond[i-1] + j] 
                      += Aii[(jjj-ct.state_begin) * ndim + kkk]*sym_fold;
                }

            }
        }

        istart += ct.block_dim[i-1];
    }


    global_sums(A_tri, &pmo.ntot, pct.grid_comm);

}

