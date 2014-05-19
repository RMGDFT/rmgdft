/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <assert.h>
#include "blacs.h"
#include "main.h"
#include "prototypes_on.h"

/*
 *
 *     set a matrix in distributed processor for
 *     parallel matrix operations.
 *
 */

void diag_eig_matrix(double *eig_matrix, double *eig, int *desca)
{
    int i, j, ii, jj, iii, jjj, li, lj;
    int iistart, jjstart, limb, ljnb;
    int mycol, myrow, nprow, npcol;
    int ictxt = desca[1], mb = desca[4], nb = desca[5], mxllda = desca[8];
    int m = desca[2], n = desca[3];
    int max_block_row, max_block_col;

    /* a_global: m x n matrix
     * mb, nb: block size in row and column
     */

    Cblacs_gridinfo (ictxt, &nprow, &npcol, &myrow, &mycol);

    max_block_row = m / (mb * nprow) +1;
    max_block_col = n / (nb * npcol) +1;


    for (li = 0; li < max_block_row; li++)
    {

        iistart = (li * nprow + myrow) * mb;
        limb = li * mb;

        for (lj = 0; lj < max_block_col; lj++)
        {

            jjstart = (lj * npcol + mycol) * nb;
            ljnb = lj * nb;

            for (i = 0; i < mb; i++)
            {

                ii = iistart + i;
                iii = i + limb;

                if (ii < m)
                {

                    for (j = 0; j < nb; j++)
                    {

                        jj = jjstart + j;
                        jjj = j + ljnb;

                        if ( jj < n )
                        {
                            eig_matrix[iii + jjj * mxllda]= 0.0;
                            if(ii == jj) eig_matrix[iii + jjj * mxllda] = eig[ii];
                        }
                    }
                }

            }
        }
    }

}
