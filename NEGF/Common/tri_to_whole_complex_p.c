/************************** SVN Revision Information **************************
 **    $Id: tri_to_whole_complex_p.c 834 2007-08-07 17:32:44Z luw $    **
******************************************************************************/
 
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "md.h"
#include "pmo.h"


void tri_to_whole_complex_p (doublecomplex * A_tri, doublecomplex * Aii, int N, int *ni)
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
     */

    int i, j,  k;
    int ndim;
    int ictxt, mb, nprow, npcol, myrow, mycol;
    int istart, jj, kk, jjj, kkk;


    ictxt = pmo.ictxt[pmo.myblacs];
    mb = pmo.mblock;

    /* dimension of matrix Aii */
    ndim = 0;
    for (i = 0; i < N; i++)
    {
        ndim += ni[i];
    }

    for(i=0; i<ndim * ndim; i++)
    {
        Aii[i].r = 0.0;
        Aii[i].i = 0.0;
    }

    Cblacs_gridinfo(ictxt, &nprow, &npcol, &myrow, &mycol);

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
                Aii[kkk * ndim + jjj].r =
                A_tri[ pmo.diag_begin[i] + k * pmo.mxllda_cond[i] + j].r;
                Aii[kkk * ndim + jjj].i =
                A_tri[ pmo.diag_begin[i] + k * pmo.mxllda_cond[i] + j].i;

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
                Aii[kkk * ndim + jjj].r =
                    A_tri[ pmo.offdiag_begin[i-1] + k * pmo.mxllda_cond[i-1] + j].r;
                Aii[kkk * ndim + jjj].i =
                    A_tri[ pmo.offdiag_begin[i-1] + k * pmo.mxllda_cond[i-1] + j].i;
                
                Aii[jjj * ndim + kkk].r = Aii[kkk * ndim + jjj].r ;
                Aii[jjj * ndim + kkk].i = Aii[kkk * ndim + jjj].i ;

            }
        }

        istart += ct.block_dim[i-1];
    }

    i = 2 * ndim * ndim;
    comm_sums((REAL *)Aii, &i, COMM_EN2);

}

