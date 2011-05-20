/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "md.h"
#include "pmo.h"


void tri_to_row (REAL * A_tri, REAL * Aii, int N, int *ni)
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
     *   Aii: store the whole matrix in input
     *   N:   number of blocks
     *   ni:  dimension of each block
     *   for example Aii is a ni[i] x ni[i] matrix
     *               Ai,i+1 is a ni[i] x ni[i+1] matrix 
     *
     *  output: Aii [ct.states_end - ct.state_begins, ct.num_states]
     */

    int i, j,  k;
    int ndim;
    int ictxt, mb, nprow, npcol, myrow, mycol;
    int istart, jj, kk, jjj, kkk, jjjj, kkkk;

    double spin_degenerate;
    
    // the density matrix multiply by 2.0 to count for spin degeneracy
    spin_degenerate = 2.0;

    ictxt = pmo.ictxt[pmo.myblacs];
    mb = pmo.mblock;

    /* dimension of matrix Aii */
    ndim = 0;
    for (i = 0; i < N; i++)
    {
        ndim += ni[i];
    }

    for(i=0; i< (ct.state_end - ct.state_begin ) * ndim; i++)
    {
        Aii[i] = 0.0;
    }

    Cblacs_gridinfo(ictxt, &nprow, &npcol, &myrow, &mycol);

    /* for diagonal blocks */

    istart = 0;
    for(i = 0; i < N; i++)
    {

        for(k=0; k < pmo.mxlocc_cond[i]; k++)
        {

            kk = (k/mb ) * npcol * mb + mycol * mb + k - k/mb * mb; 
            kkk = kk + istart;
            if(kkk >= ct.state_begin && kkk < ct.state_end)
            {

                kkkk = kkk - ct.state_begin;
                for(j =0; j < pmo.mxllda_cond[i]; j++)
                {

                    /* for each block, distributed local index (j,k) is (jj, kk) in block matrix
                     * and (jjj,kkk) in whole matrix
                     */

                    jj = (j/mb ) * nprow * mb + myrow * mb + j - j/mb * mb; 

                    jjj = jj + istart;
                    Aii[kkkk * ndim + jjj] = spin_degenerate *
                        A_tri[ pmo.diag_begin[i] + k * pmo.mxllda_cond[i] + j];

                }
            }
        }

        istart += ct.block_dim[i];
    }


    /* for upper off-diagonal blocks */

    istart = 0;
    for(i = 1; i < N; i++)
    {

            for(k=0; k < pmo.mxlocc_cond[i]; k++)
            {
                kk = (k/mb ) * npcol * mb + mycol * mb + k - k/mb * mb; 
                kkk = kk + istart + ct.block_dim[i-1];
                if(kkk >= ct.state_begin && kkk < ct.state_end)
                {

                    kkkk = kkk - ct.state_begin;
                    for(j =0; j < pmo.mxllda_cond[i-1]; j++)
                    {

                        /* for each block, distributed local index (j,k) is (jj, kk) in block matrix
                         * and (jjj,kkk) in whole matrix
                         */

                        jj = (j/mb ) * nprow * mb + myrow * mb + j - j/mb * mb; 

                        jjj = jj + istart;
                        Aii[kkkk * ndim + jjj] = spin_degenerate *
                            A_tri[ pmo.offdiag_begin[i-1] + k * pmo.mxllda_cond[i-1] + j];


                    }
                }
            }

            istart += ct.block_dim[i-1];
    }



    /* for lower off-diagonal blocks */

    istart = 0;
    for(i = 1; i < N; i++)
    {

        for(j =0; j < pmo.mxllda_cond[i-1]; j++)
        {

            jj = (j/mb ) * nprow * mb + myrow * mb + j - j/mb * mb; 
            jjj = jj + istart;
            if(jjj >= ct.state_begin && jjj < ct.state_end)
            {

                jjjj = jjj - ct.state_begin;
                for(k=0; k < pmo.mxlocc_cond[i]; k++)
                {
                    kk = (k/mb ) * npcol * mb + mycol * mb + k - k/mb * mb; 
                    kkk = kk + istart + ct.block_dim[i-1];

                    /* for each block, distributed local index (j,k) is (jj, kk) in block matrix
                     * and (jjj,kkk) in whole matrix
                     */
                    Aii[jjjj * ndim + kkk] = spin_degenerate *
                        A_tri[ pmo.offdiag_begin[i-1] + k * pmo.mxllda_cond[i-1] + j];


                }
            }
        }

        istart += ct.block_dim[i-1];
    }

    i = (ct.state_end - ct.state_begin) * ndim;
    comm_sums(Aii, &i, COMM_EN2);

}

