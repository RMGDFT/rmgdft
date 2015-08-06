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


void whole_to_tri_real (double * A_tri, double * Aii, int N, int *ni)
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
 *   Ai,i+1  = transpose (Ai+1, i)
 *
 *   Aii: store the whole matrix in input
 *   N:   number of blocks
 *   ni:  dimension of each block
 *   for example Aii is a ni[i] x ni[i] matrix
 *               Ai,i+1 is a ni[i] x ni[i+1] matrix 
 */

    int i, j, ntem, k;
    int *n_begin;
    int ione = 1;
    int ntot, ndim;


/*  find the maximum dimension of the blocks  */


    ntot = 0;
    for (i = 0; i < N; i++)
    {
        ntot += ni[i];
    }

    ndim = 0;
    for (i = 0; i < N - 1; i++)
    {
        ndim += ni[i] * ni[i] + ni[i] * ni[i + 1];
    }

    ndim += ni[N - 1] * ni[N - 1];


    my_malloc_init( n_begin, N, int );


/*  n_begin: starting address of each diagonal block in A_tri and G_tri
 *  the Hi,i+1 block will start at n_begin[i] + ni[i] * ni[i]
 */

    n_begin[0] = 0;
    for (i = 1; i < N; i++)
    {
        n_begin[i] = n_begin[i - 1] + ni[i - 1] * ni[i - 1] + ni[i - 1] * ni[i];
    }

    for (i = 0; i < ndim; i++)
    {
        A_tri[i] = 0.0;
    }

    ntem = 0;
    for (i = 0; i < N; i++)
    {

        for (k = 0; k < ni[i]; k++)
            for (j = 0; j < ni[i]; j++)
            {
                A_tri[n_begin[i] + j * ni[i] + k] = Aii[(j + ntem) * ntot + k + ntem];

            }

        ntem += ni[i];
    }

    ntem = 0;
    for (i = 0; i < N - 1; i++)
    {

        for (k = 0; k < ni[i + 1]; k++)
            for (j = 0; j < ni[i]; j++)
            {
                A_tri[n_begin[i] + ni[i] * ni[i] + k * ni[i] + j] =
                    Aii[j + ntem + (k + ni[i] + ntem) * ntot];

            }

        ntem += ni[i];
    }
    my_free( n_begin );
}
