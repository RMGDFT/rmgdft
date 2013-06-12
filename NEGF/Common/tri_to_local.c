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

    int ictxt, mb, nprow, npcol, myrow, mycol;
    
    int size;
    int idx2, idx1;

    double spin_degenerate;

    int st1, st2, block1, block2, st1_local, st2_local, col_proc, row_proc;
    int st1_index, st2_index;
    int ione = 1;
    
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
    

    for(st1 = 0; st1 < pct.num_local_orbit * pct.num_local_orbit; st1++) Aii_local[st1] = 0.0;

    for(st1 = 0; st1 < pct.num_local_orbit; st1++)
    {
        block1 = states_distribute[st1].whichblock;
        st1_index = states_distribute[st1].istate_in_block;

        col_proc = (st1_index/mb)%npcol;

        st1_local = (st1_index/mb)/npcol *mb + st1_index%mb;

        if(col_proc == mycol)
        {

            for(st2 = st1; st2 < pct.num_local_orbit; st2++)
            {
                block2 = states_distribute[st2].whichblock;
                st2_index = states_distribute[st2].istate_in_block;

                row_proc = (st2_index/mb)%nprow;

                st2_local = (st2_index/mb)/nprow *mb + st2_index%mb;

                if(row_proc == myrow)
                {

                    if(block1 == block2) 
                    {
                        idx1 = st1 * pct.num_local_orbit + st2;
                        idx2 = pmo.diag_begin[block1] + st1_local * pmo.mxllda_cond[block1] + st2_local;

                        Aii_local[idx1] = A_tri[idx2];
                        Aii_local[st2 * pct.num_local_orbit + st1 ] = Aii_local[idx1];

                    }
                    else if(block1 +1 == block2)
                    {
                        idx1 = st1 * pct.num_local_orbit + st2;
                        idx2 = pmo.offdiag_begin[block1] + st2_local * pmo.mxllda_cond[block1] + st1_local;

                        Aii_local[idx1] = A_tri[idx2];
                        Aii_local[st2 * pct.num_local_orbit + st1 ] = Aii_local[idx1];

                    }
                    else
                    {

                        Aii_local[st1 * pct.num_local_orbit + st2 ] = 0.0;
                        Aii_local[st2 * pct.num_local_orbit + st1 ] = 0.0;
                    }
                         
                }
            }
        }
    }




    size = pct.num_local_orbit * pct.num_local_orbit;
    if(size > 0)
    {
        comm_sums(Aii_local, &size, COMM_EN2);

        dscal(&size, &spin_degenerate, Aii_local, &ione);
    }

}

