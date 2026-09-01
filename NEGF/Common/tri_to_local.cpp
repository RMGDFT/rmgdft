/************************** SVN Revision Information **************************
 **    $Id: tri_to_local.c 3140 2015-08-06 15:48:24Z luw $    **
 ******************************************************************************/

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
    
#include "blacs.h"
#include "blas.h"
#include "main.h"
#include "init_var.h"
#include "LCR.h"
#include "pmo.h"
#include "LocalObject.h"
#include "GlobalSums.h"



void tri_to_local (double * A_tri, double * Aii_local, LocalObject<double> &Phi)
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
     *  output: Aii_local, depending on how Phi localized, it only  has the necessary matrix elements.
     *                     dimentsion Phi.num_thispe * Phi.num_thispe
     */

    int size;

    double spin_degenerate;

    int ione = 1;

    int i, j,  k;
    int ictxt, mb, nprow, npcol, myrow, mycol;
    int jj, kk;

    // the density matrix multiply by 2.0 to count for spin degeneracy
    spin_degenerate = 1.0;
    if(ct.nspin == 1) spin_degenerate = 2.0;

    ictxt = pmo.ictxt[pmo.myblacs];
    mb = pmo.mblock;

    Cblacs_gridinfo(ictxt, &nprow, &npcol, &myrow, &mycol);

    // in communictor COMM_EN2, there are nprow * npcol process, the
    // matrix are distributed over these process. depending their rank
    // in pct.grid_comm, each process have different orbitals, and

    //determing the starting and ending ranks in pct.grid_comm for these
    //nprow * npcol PE



    // change the tri-diagonal distributed matrix to global matrix, need
    // to be optimized later for memory issues with matrix_tem.
    //

    /* for diagonal blocks */
    for(int idx = 0; idx < Phi.num_thispe * Phi.num_thispe; idx++) Aii_local[idx] = 0.0;

    int max_block_size = *std::max_element(ct.block_dim, ct.block_dim+ct.num_blocks);
    double *matrix_tem = new double[max_block_size * max_block_size];

    for(i = 0; i < ct.num_blocks; i++)
    {
        for(j=0; j<ct.block_dim[i] * ct.block_dim[i]; j++) 
        {
            matrix_tem[j] = 0.0;
        }

        for(j =0; j < pmo.mxllda_cond[i]; j++)
        {
            for(k=0; k < pmo.mxlocc_cond[i]; k++)
            {

                /* for each block, distributed local index (j,k) is (jj, kk) in block matrix
                 * and (jjj,kkk) in whole matrix
                 */

                jj = (j/mb ) * nprow * mb + myrow * mb + j - j/mb * mb; 
                kk = (k/mb ) * npcol * mb + mycol * mb + k - k/mb * mb; 

                matrix_tem[kk * ct.block_dim[i] + jj] =
                    A_tri[ pmo.diag_begin[i] + k * pmo.mxllda_cond[i] + j];

            }
        }

        size = ct.block_dim[i] * ct.block_dim[i];
        GlobalSums(matrix_tem, size, COMM_EN2);

        for(j =0; j < ct.block_dim[i]; j++)
        {
            for(k=0; k < ct.block_dim[i]; k++)
            {
                jj = Phi.index_global_to_proj[pmo.orb_index[i] + j];
                kk = Phi.index_global_to_proj[pmo.orb_index[i] + k];

                if(jj >= 0 && kk >=0)
                    Aii_local[kk * Phi.num_thispe + jj] = matrix_tem[k * ct.block_dim[i] + j] ;


            }
        }


    }


    /* for off-diagonal blocks */

    for(i = 1; i < ct.num_blocks; i++)
    {
        for(j=0; j<ct.block_dim[i-1] * ct.block_dim[i]; j++) 
        {
            matrix_tem[j] = 0.0;
        }

        for(j =0; j < pmo.mxllda_cond[i-1]; j++)
        {
            for(k=0; k < pmo.mxlocc_cond[i]; k++)
            {

                /* for each block, distributed local index (j,k) is (jj, kk) in block matrix
                 * and (jjj,kkk) in whole matrix
                 */

                jj = (j/mb ) * nprow * mb + myrow * mb + j - j/mb * mb; 
                kk = (k/mb ) * npcol * mb + mycol * mb + k - k/mb * mb; 

                matrix_tem[kk * ct.block_dim[i-1] + jj] =
                    A_tri[ pmo.offdiag_begin[i-1] + k * pmo.mxllda_cond[i-1] + j];


            }
        }

        size = ct.block_dim[i-1] * ct.block_dim[i];
        GlobalSums(matrix_tem, size, COMM_EN2);


        for(j =0; j < ct.block_dim[i-1]; j++)
        {
            for(k=0; k < ct.block_dim[i]; k++)
            {
                jj = Phi.index_global_to_proj[pmo.orb_index[i-1] + j];
                kk = Phi.index_global_to_proj[pmo.orb_index[i] + k];

                if(jj >= 0 && kk >=0)
                {
                    Aii_local[kk * Phi.num_thispe + jj] = matrix_tem[k * ct.block_dim[i-1] + j] ;
                    Aii_local[jj * Phi.num_thispe + kk] = matrix_tem[k * ct.block_dim[i-1] + j] ;
                }


            }
        }

    }



    size = Phi.num_thispe * Phi.num_thispe;
    if(size > 0)
    {
        dscal(&size, &spin_degenerate, Aii_local, &ione);
    }

    delete [] matrix_tem;
}

