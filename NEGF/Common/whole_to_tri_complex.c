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


void whole_to_tri_complex (complex double * H_tri, complex double * Hii, int N, int * ni)
{
/* Semi_tridiagonal matrix  
 *
 *    H00   H01   0    0   ...		0
 *    H10   H11   H12  0  ...		0	
 *    0     H21   H22  H23 ....		0
 *    .     .     .    .		.	
 *    .     .     .    .		.
 *    .     .     .    .     Hn-1,n-1	Hn-1,n
 *    0     0     0          Hn,n-1     Hnn
 *
 *   H_tri: output store the input matrix in the order of H00, H01, H11, H12, .... Hnn
 *   Hi,i+1  = transpose (Hi+1, i)
 *
 *   Hii: store the whole matrix in input
 *   N:   number of blocks
 *   ni:  dimension of each block
 *   for example Hii is a ni[i] x ni[i] matrix
 *               Hi,i+1 is a ni[i] x ni[i+1] matrix 
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


/*  n_begin: starting address of each diagonal block in H_tri and G_tri
 *  the Hi,i+1 block will start at n_begin[i] + ni[i] * ni[i]
 */

    n_begin[0] = 0;
    for (i = 1; i < N; i++)
    {
        n_begin[i] = n_begin[i - 1] + ni[i - 1] * ni[i - 1] + ni[i - 1] * ni[i];
    }

    for (i = 0; i < ndim; i++)
    {
        H_tri[i] = 0.0;
    }

    ntem = 0;
    for (i = 0; i < N; i++)
    {

        for (k = 0; k < ni[i]; k++)
            for (j = 0; j < ni[i]; j++)
            {
                H_tri[n_begin[i] + j * ni[i] + k] = Hii[(j + ntem) * ntot + k + ntem];

            }

        ntem += ni[i];
    }

    ntem = 0;
    for (i = 0; i < N - 1; i++)
    {

        for (k = 0; k < ni[i + 1]; k++)
            for (j = 0; j < ni[i]; j++)
            {
                H_tri[n_begin[i] + ni[i] * ni[i] + k * ni[i] + j] =
                    Hii[j + ntem + (k + ni[i] + ntem) * ntot];

            }

        ntem += ni[i];
    }
    my_free( n_begin );
}
