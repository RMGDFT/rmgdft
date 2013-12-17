/************************** SVN Revision Information **************************
 **    $Id$    **
 ******************************************************************************/

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "main.h"
#include "pmo.h"


void tri_to_local (STATE *states_distribute, REAL * A_tri, REAL * Aii_local)
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

    int size;
    int idx;

    double spin_degenerate;

    int st1, st2, st11, st22;
    int ione = 1;

    int i, j,  k;
    int ictxt, mb, nprow, npcol, myrow, mycol;
    int istart, jj, kk, jjj, kkk;

    // the density matrix multiply by 2.0 to count for spin degeneracy
    spin_degenerate = 2.0;

    ictxt = pmo.ictxt[pmo.myblacs];
    mb = pmo.mblock;

    Cblacs_gridinfo(ictxt, &nprow, &npcol, &myrow, &mycol);

    // in communictor COMM_EN2, there are nprow * npcol process, the
    // matrix are distributed over these process. depending their rank
    // in pct.grid_comm, each process have different orbitals, and
    // returned Aii_row only storeis (ct.state_end-ct.state_begin) * ndim, 
    // not the whole matrix ndim * ndim

    //determing the starting and ending ranks in pct.grid_comm for these
    //nprow * npcol PE



    // change the tri-diagonal distributed matrix to global matrix, need
    // to be optimized later for memory issues with work_matrix.
    //
    for(i=0; i<ct.num_states * ct.num_states; i++)
    {
        work_matrix[i] = 0.0;
    }



    /* for diagonal blocks */

    istart = 0;
    for(i = 0; i < ct.num_blocks; i++)
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
                work_matrix[kkk * ct.num_states + jjj] =
                    A_tri[ pmo.diag_begin[i] + k * pmo.mxllda_cond[i] + j];

            }
        }

        istart += ct.block_dim[i];
    }



    /* for off-diagonal blocks */

    istart = 0;
    for(i = 1; i < ct.num_blocks; i++)
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
                work_matrix[kkk * ct.num_states + jjj] =
                    A_tri[ pmo.offdiag_begin[i-1] + k * pmo.mxllda_cond[i-1] + j];

                work_matrix[jjj * ct.num_states + kkk] = work_matrix[kkk * ct.num_states + jjj] ;

            }
        }

        istart += ct.block_dim[i-1];
    }

    size = ct.num_states * ct.num_states;
    comm_sums((REAL *)work_matrix, &size, COMM_EN2);


    // global matrix map to local matrix(pct.num_local * pct.num_local);
    //
    for(st1 = 0; st1 < pct.num_local_orbit * pct.num_local_orbit; st1++) Aii_local[st1] = 0.0;


    for (st1 = 0; st1 < pct.num_local_orbit; st1++)
    {
        st11 = states_distribute[st1].istate;
        for (st2 = 0; st2 < pct.num_local_orbit; st2++)
        {
            st22 = states_distribute[st2].istate;

            idx = st11 * ct.num_states + st22;
            Aii_local[st1 * pct.num_local_orbit + st2] = work_matrix[idx];

        }
    }                           /* end for st1 = .. */




    size = pct.num_local_orbit * pct.num_local_orbit;
    if(size > 0)
    {
        dscal(&size, &spin_degenerate, Aii_local, &ione);
    }

}

