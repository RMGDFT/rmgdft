/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "main.h"


double trace_AB_tri(REAL * H_tri, REAL * par_D_tri, int N, int *ni)
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
 *   H_tri and par_D_tri store the matrix in the order of A00, A01, A11, A12, .... Ann
 *   Ai,i+1  = transpose (Ai+1, i)
 *
 *   H_tri: tridiagonal matrix storing Halmitonian
 *   par_D_tri: tridiagonal matrix storing partial charge density matrix
 *   N:   number of blocks
 *   ni:  dimension of each block
 *   for example Aii is a ni[i] x ni[i] matrix
 *               Ai,i+1 is a ni[i] x ni[i+1] matrix 
 */

    int i, j, idx;
    int *n_begin;
    double product;


    my_malloc_init( n_begin, N, int );

/*  n_begin: starting address of each diagonal block in A_tri and G_tri
 *  the Hi,i+1 block will start at n_begin[i] + ni[i] * ni[i]
 */

    n_begin[0] = 0;
    for (i = 1; i < N; i++)
    {
        n_begin[i] = n_begin[i - 1] + ni[i - 1] * ni[i - 1] + ni[i - 1] * ni[i];
    }


    product = 0.0;
    for (i = 1; i < N; i++)
    {
        idx = n_begin[i-1];
        for (j = 0; j < ni[i-1] * ni[i-1]; j++)
           product += H_tri[idx + j] * par_D_tri[idx + j];
    
        idx += ni[i-1] * ni[i-1];
        for (j = 0; j < ni[i-1] * ni[i]; j++)
            product += 2.0 * H_tri[idx + j] * par_D_tri[idx + j];
    }

    idx = n_begin[N-1];
    for(j=0; j < ni[N-1] * ni[N-1]; j++)
        product += H_tri[idx + j] * par_D_tri[idx + j];

    return (product);

    my_free( n_begin );
}
