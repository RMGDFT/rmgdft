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
    int n2, ndim;
    int ictxt, mb, nprow, npcol, myrow, mycol;
    int istart, jj, kk;


    ictxt = pmo.ictxt[pmo.myblacs];
    mb = pmo.mblock;

    /* dimension of matrix Aii */
    ndim = 0;
    for (i = 0; i < N; i++)
    {
        ndim += ni[i];
    }

    Cblacs_gridinfo(ictxt, &nprow, &npcol, &myrow, &mycol);


    /* for diagonal blocks */

    istart = 0;
    for(i = 0; i < N; i++)
    {

        n2 = ct.block_dim[i] * ct.block_dim[i];

        for(j = 0; j < n2; j++ ) work_matrix[j] = 0.0;

        for(j = 0; j < ct.block_dim[i]; j++)
            for(k = 0; k < ct.block_dim[i]; k++)
            {
                jj = istart + j;
                kk = istart + k;

                if(kk >= ct.state_begin && kk < ct.state_end)
                {
                    work_matrix[k * ct.block_dim[i] + j] += 0.5 *Aii[ (kk-ct.state_begin) * ndim + jj];
                    work_matrix[j * ct.block_dim[i] + k] += 0.5 *Aii[ (kk-ct.state_begin) * ndim + jj];
                }
            }


        global_sums(work_matrix, &n2, pct.grid_comm);


        for(j =0; j < pmo.mxllda_cond[i]; j++)
        {
            for(k=0; k < pmo.mxlocc_cond[i]; k++)
            {

                /* for each block, distributed local index (j,k) is (jj, kk) in block matrix
                 * and (jjj,kkk) in whole matrix
                 */

                jj = (j/mb ) * nprow * mb + myrow * mb + j - j/mb * mb; 
                kk = (k/mb ) * npcol * mb + mycol * mb + k - k/mb * mb; 


                A_tri[ pmo.diag_begin[i] + k * pmo.mxllda_cond[i] + j] 
                    += work_matrix[k * ct.block_dim[i] + j];
            }

        }

        istart += ct.block_dim[i];
    }



    /* for off-diagonal blocks */

    istart = 0;
    for(i = 1; i < N; i++)
    {
        n2 = ct.block_dim[i-1] * ct.block_dim[i];

        for(j = 0; j < n2; j++ ) work_matrix[j] = 0.0;

        for(j = 0; j < ct.block_dim[i-1]; j++)
            for(k = 0; k < ct.block_dim[i]; k++)
            {
                jj = istart + j;
                kk = istart + k + ct.block_dim[i-1];

                if(kk >= ct.state_begin && kk < ct.state_end)
                {
                    work_matrix[k * ct.block_dim[i-1] + j] += 0.5 *Aii[ (kk-ct.state_begin) * ndim + jj];
                }
                if(jj >= ct.state_begin && jj < ct.state_end)
                {
                    work_matrix[k * ct.block_dim[i-1] + j] += 0.5 *Aii[ (jj-ct.state_begin) * ndim + kk];
                }
            }

        global_sums(work_matrix, &n2, pct.grid_comm);

        for(j =0; j < pmo.mxllda_cond[i-1]; j++)
        {
            for(k=0; k < pmo.mxlocc_cond[i]; k++)
            {

                /* for each block, distributed local index (j,k) is (jj, kk) in block matrix
                 * and (jjj,kkk) in whole matrix
                 */

                jj = (j/mb ) * nprow * mb + myrow * mb + j - j/mb * mb; 
                kk = (k/mb ) * npcol * mb + mycol * mb + k - k/mb * mb; 

                A_tri[ pmo.offdiag_begin[i-1] + k * pmo.mxllda_cond[i-1] + j] 
                    += work_matrix[kk * ct.block_dim[i-1]  + jj] ;

            }
        }

        istart += ct.block_dim[i-1];
    }


}

