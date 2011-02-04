/************************** SVN Revision Information **************************
 **    $Id: tri_to_whole_complex.c 1242 2011-02-02 18:55:23Z luw $    **
******************************************************************************/
 

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "md.h"

void tri_to_whole_complex (doublecomplex * H_tri, doublecomplex * Hii, int N, int * ni)
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
 *   H_tri: store the input matrix in the order of H00, H01, H11, H12, .... Hnn
 *   Hi,i+1  = transpose (Hi+1, i)
 *
 *   Hii: store the whole matrix in output
 *   N:   number of blocks
 *   ni:  dimension of each block
 *   for example Hii is a ni[i] x ni[i] matrix
 *               Hi,i+1 is a ni[i] x ni[i+1] matrix 
 */

    int i, j, ntem, k;
    int *n_begin;
    int ione = 1;
    int ntot;


/*  find the maximum dimension of the blocks  */


    ntot = 0;
    for (i = 0; i < N; i++)
    {
        ntot += ni[i];
    }

    my_malloc_init( n_begin, N, int );


/*  n_begin: starting address of each diagonal block in H_tri and G_tri
 *  the Hi,i+1 block will start at n_begin[i] + ni[i] * ni[i]
 *  ntem_begin: starting address of one column of G_tem for each block
 */

    n_begin[0] = 0;
    for (i = 1; i < N; i++)
    {
        n_begin[i] = n_begin[i - 1] + ni[i - 1] * ni[i - 1] + ni[i - 1] * ni[i];
    }

    for (i = 0; i < ntot * ntot; i++)
    {
        Hii[i].r = 0.0, Hii[i].i = 0.0;
    }

    ntem = 0;
    for (i = 0; i < N; i++)
    {

        for (j = 0; j < ni[i]; j++)
            for (k = 0; k < ni[i]; k++)
            {
                Hii[(j + ntem) * ntot + k + ntem].r = H_tri[n_begin[i] + j * ni[i] + k].r;
                Hii[(j + ntem) * ntot + k + ntem].i = H_tri[n_begin[i] + j * ni[i] + k].i;

            }

        ntem += ni[i];
    }


    ntem = 0;
    for (i = 0; i < N - 1; i++)
    {

        for (j = 0; j < ni[i]; j++)
            for (k = 0; k < ni[i + 1]; k++)
            {
                Hii[j + ntem + (k + ntem + ni[i]) * ntot].r =
                    H_tri[n_begin[i] + ni[i] * ni[i] + k * ni[i] + j].r;
                Hii[(j + ntem) * ntot + k + ntem + ni[i]].r =
                    H_tri[n_begin[i] + ni[i] * ni[i] + k * ni[i] + j].r;

                Hii[j + ntem + (k + ntem + ni[i]) * ntot].i =
                    H_tri[n_begin[i] + ni[i] * ni[i] + k * ni[i] + j].i;
                Hii[(j + ntem) * ntot + k + ntem + ni[i]].i =
                    H_tri[n_begin[i] + ni[i] * ni[i] + k * ni[i] + j].i;


            }

        ntem += ni[i];
    }
    my_free( n_begin );
}
