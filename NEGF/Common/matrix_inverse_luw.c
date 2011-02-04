/************************** SVN Revision Information **************************
 **    $Id: matrix_inverse_luw.c 1242 2011-02-02 18:55:23Z luw $    **
******************************************************************************/
 
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "md.h"


void matrix_inverse_luw (doublecomplex * H_tri, doublecomplex * G_tri, int N, int * ni)
{
/*  Calculate the inverse of a semi-tridiagonal complex matrix
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
 *   G_tri: output inverse matrix in the same order as H_tri
 *   N:   number of blocks
 *   ni:  dimension of each block
 *   for example Hii is a ni[i] x ni[i] matrix
 *               Hi,i+1 is a ni[i] x ni[i+1] matrix 
 */

    int nmax, i, j, n1, n2, n3, n4, n5, n6, n7, n8;
    int *ipiv, *n_begin, *ntem_begin;
    doublecomplex *Hii, *Gii, *G_tem, *temp, *Hii1, *Gii0;
    doublecomplex mone, one, zero;
    int ione = 1;
    int ntot;

    mone.r = -1.0, mone.i = 0.0;
    one.r = 1.0, one.i = 0.0;
    zero.r = 0.0, zero.i = 0.0;

/*  find the maximum dimension of the blocks  */


    nmax = 0;
    ntot = 0;
    for (i = 0; i < N; i++)
    {
        ntot += ni[i];
        if (nmax < ni[i])
            nmax = ni[i];
    }

    my_malloc_init( n_begin, N, int );
    my_malloc_init( ntem_begin, N, int );
    my_malloc_init( ipiv, nmax, int );


    my_malloc_init( Hii, nmax * nmax, doublecomplex );
    my_malloc_init( Gii, nmax * nmax, doublecomplex );
    my_malloc_init( temp, nmax * nmax, doublecomplex );
    my_malloc_init( G_tem, nmax * ntot, doublecomplex );


/*  n_begin: starting address of each diagonal block in H_tri and G_tri
 *  the Hi,i+1 block will start at n_begin[i] + ni[i] * ni[i]
 *  ntem_begin: starting address of one column of G_tem for each block
 */

    n_begin[0] = 0;
    ntem_begin[0] = 0;
    for (i = 1; i < N; i++)
    {
        n_begin[i] = n_begin[i - 1] + ni[i - 1] * ni[i - 1] + ni[i - 1] * ni[i];
        ntem_begin[i] = ntem_begin[i - 1] + ni[i - 1] * nmax;
    }

/*  calculate the inverse of the first block  */

    for (i = 0; i < ni[0] * ni[0]; i++)
    {
        Hii[i].r = H_tri[i].r;
        Hii[i].i = H_tri[i].i;
    }

    get_inverse_block (Hii, Gii, ipiv, ni[0]);

    n1 = ni[0] * ni[0];
    zcopy (&n1, Gii, &ione, G_tri, &ione);

/*  iterate to get one more block  */

    for (i = 0; i < N - 1; i++)
    {
        /* get the interaction  Hi,i+1  from input H_tri 
         * Hii1 is a pointer only
         */
        Hii1 = &H_tri[n_begin[i] + ni[i] * ni[i]];
        Gii0 = &G_tri[n_begin[i]];

        /* Hii now has the matrix Hi+1,i+1  */

        for (j = 0; j < ni[i + 1] * ni[i + 1]; j++)
        {
            Hii[j].r = H_tri[j + n_begin[i + 1]].r;
            Hii[j].i = H_tri[j + n_begin[i + 1]].i;
        }


        /* calculate Hi+1,i+1 - Hi+1,i * Gii^0 * Hi,i+1  */

        n1 = ni[i + 1];
        n2 = ni[i];
        ZGEMM ("T", "N", &n1, &n2, &n2, &one, Hii1, &n2, Gii0, &n2, &zero, temp, &n1);
        ZGEMM ("N", "N", &n1, &n1, &n2, &mone, temp, &n1, Hii1, &n2, &one, Hii, &n1);


        /* now Hii store the matrix Hi+1,i+1 - Hi+1,i * Gii^0 * Hi,i+1
         * Gi+1,i+1, stored in Gii, = Hii^(-1)
         */

        get_inverse_block (Hii, Gii, ipiv, ni[i + 1]);

        n1 = ni[i + 1] * ni[i + 1];
        zcopy (&n1, Gii, &ione, &G_tri[n_begin[i + 1]], &ione);

        /*  Gi,i+1 =  Gii^0 * Hi,i+1 * Gi+1,i+1  */

        n1 = ni[i];
        n2 = ni[i + 1];
        n3 = n_begin[i] + ni[i] * ni[i];

        /* temp = Gii^0 * Hi,i+1 */
        ZGEMM ("N", "N", &n1, &n2, &n1, &mone, Gii0, &n1, Hii1, &n1, &zero, temp, &n1);

        /* G(i,i+1) = temp * G(i+1,i+1)  also == G_tem(i,i+1)  */
        ZGEMM ("N", "N", &n1, &n2, &n2, &one, temp, &n1, Gii, &n2, &zero, &G_tri[n3], &n1);

        n4 = ni[i] * ni[i + 1];
        zcopy (&n4, &G_tri[n3], &ione, &G_tem[ntem_begin[i]], &ione);

        /* update Gii  */

        ZGEMM ("N", "T", &n1, &n1, &n2, &one, temp, &n1, &G_tri[n3], &n1, &one,
               &G_tri[n_begin[i]], &n1);

        for (j = i - 1; j >= 0; j--)
        {

            n1 = ni[j];         /* dimension of block j */
            n2 = ni[i];         /* dimension of block i */
            n3 = ni[i + 1];     /* dimension of block i+1 */
            n4 = ntem_begin[j]; /* starting address of Gtem(j,*) in G_tem */
            n5 = n_begin[j];    /* starting address of G(j,j) in G_tri */
            n6 = n_begin[j] + ni[j] * ni[j];    /* starting address of G(j,j+1) in G_tri */
            n7 = ntem_begin[j + 1];     /* starting address of Gtem(j+1,*) in G_tem */
            n8 = ni[j + 1];     /* dimension of block j+1 */

            /* temp = -G0(j,i) * Hi,i+1  */
            ZGEMM ("N", "N", &n1, &n3, &n2, &mone, &G_tem[n4], &n1, Hii1, &n2, &zero, temp, &n1);

            /* G0(j, i+1) = temp * G(i+1,i+1) */
            ZGEMM ("N", "N", &n1, &n3, &n3, &one, temp, &n1, Gii, &n3, &zero, &G_tem[n4], &n1);

            /* G(j,j) = G0(j,j) + temp * G(i+1,j) */
            ZGEMM ("N", "T", &n1, &n1, &n3, &one, temp, &n1, &G_tem[n4], &n1, &one, &G_tri[n5],
                   &n1);

            /* G(j,j+1) = G0(j,j+1) + temp * G(i+1,j+1)  */
            ZGEMM ("N", "T", &n1, &n8, &n3, &one, temp, &n1, &G_tem[n7], &n8, &one, &G_tri[n6],
                   &n1);
        }                       /* end for (j--) */

    }                           /* end  for(i = 0; i < N-1; i++) */


    my_free( n_begin );
    my_free( ntem_begin );
    my_free( ipiv );
    my_free( Hii );
    my_free( Gii );
    my_free( temp );
    my_free( G_tem );
}
