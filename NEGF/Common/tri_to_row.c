/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "main.h"
#include "pmo.h"


void tri_to_row (REAL * A_tri, REAL * Aii_row, int N, int *ni)
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
    
    int size, npe_en2, pe_start, pe_end, num_orbital_thisgroup;
    int idx, idx1;

    double spin_degenerate, *Aii;
    
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
    
    npe_en2 = nprow * npcol;
    pe_start = pct.gridpe / npe_en2 * npe_en2;
    pe_end = pe_start + npe_en2 -1 ;   


    // total number of orbitals in these PE group

    num_orbital_thisgroup = state_end[pe_end] - state_begin[pe_start];


    /* dimension of matrix Aii */
    ndim = 0;
    for (i = 0; i < N; i++)
    {
        ndim += ni[i];
    }

    size = num_orbital_thisgroup * ndim;
    my_malloc_init( Aii, size, double );


    /* for diagonal blocks */

    istart = 0;
    for(i = 0; i < N; i++)
    {

        for(k=0; k < pmo.mxlocc_cond[i]; k++)
        {

            kk = (k/mb ) * npcol * mb + mycol * mb + k - k/mb * mb; 
            kkk = kk + istart;
            if(kkk >= state_begin[pe_start] && kkk < state_end[pe_end])
            {

                kkkk = kkk - state_begin[pe_start];
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
                if(kkk >= state_begin[pe_start] && kkk < state_end[pe_end])
                {

                    kkkk = kkk - state_begin[pe_start];
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
            if(jjj >= state_begin[pe_start] && jjj < state_end[pe_end])
            {

                jjjj = jjj - state_begin[pe_start];
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

    comm_sums(Aii, &size, COMM_EN2);

    //assign Aii to Aii_row, the (ct.state_end - ct.state_begin) * ndim
    //

    for( i = ct.state_begin; i < ct.state_end; i++)
        for(j = 0; j < ndim; j++)
        {

            idx = (i - ct.state_begin ) * ndim +j;
            idx1 = (i - state_begin[pe_start]) * ndim +j;
            Aii_row [idx] = Aii [idx1];
        }

    my_free(Aii);

}

