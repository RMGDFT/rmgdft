/************************** SVN Revision Information **************************
 **    $Id: matrix_inverse_right_p.c 1242 2011-02-02 18:55:23Z luw $    **
******************************************************************************/
 

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "md.h"
#include "pmo.h"




void matrix_inverse_right_p (doublecomplex * H_tri, int N, int * ni, doublecomplex * Green_C)
{
/*  Calculate the inverse of a semi-tridiagonal complex matrix
 *
 *    H00   H01   0    0   ...      0
 *    H10   H11   H12  0  ...       0   
 *    0     H21   H22  H23 ....     0
 *    .     .     .    .        .   
 *    .     .     .    .        .
 *    .     .     .    .     Hn-1,n-1   Hn-1,n
 *    0     0     0          Hn,n-1     Hnn
 *
 *   H_tri: store the input matrix in the order of H00, H01, H11, H12, .... Hnn
 *   Hi,i+1  = transpose (Hi+1, i)
 *
 *   ch0: work as temperary memory
 *        output: store the inversed matrix, but only the first row of blocks are right 
 *        which will be used for non-equilibrium part calculation Gamma *G* Gamma for left lead
 *   N:   number of blocks
 *   ni:  dimension of each block
 *   for example Hii is a ni[i] x ni[i] matrix
 *               Hi,i+1 is a ni[i] x ni[i+1] matrix 
 */

    int nmax, i, j, n1, n2, n3, n4, n5, n6, n7, n8;
    int *ipiv, *n_begin1;
    doublecomplex *Hii, *Gii, *G_tem, *temp, *Hii1;
    doublecomplex mone, one, zero;
    int ione = 1;
    int ntot, k, maxrow, maxcol, *desca, *descb, *descc;
    int *descd;

    mone.r = -1.0, mone.i = 0.0;
    one.r = 1.0, one.i = 0.0;
    zero.r = 0.0, zero.i = 0.0;

/*  find the maximum dimension of the blocks  */


    ntot = pmo.ntot;
    maxrow = 0; 
    maxcol = 0;
    for (i = 0; i < N; i++)
    {
        maxrow = max(maxrow, pmo.mxllda_cond[i]);
        maxcol = max(maxcol, pmo.mxlocc_cond[i]);
    }

    my_malloc( n_begin1, N, int );
    my_malloc( ipiv, maxrow + pmo.mblock, int );

    my_malloc_init( Hii, maxrow * maxcol, doublecomplex );
    my_malloc_init( Gii, maxrow * maxcol, doublecomplex );
    my_malloc_init( temp, maxrow * maxcol, doublecomplex );
    G_tem = Green_C;



/*  n_begin: starting address of each diagonal block in H_tri and G_tri
 *  the Hi,i+1 block will start at n_begin[i] + ni[i] * ni[i]
 */

    n_begin1[0] = 0;
    for (i = 1; i < N; i++)
    {
        n_begin1[i] = n_begin1[i - 1] + pmo.mxllda_cond[i - 1] * maxcol;
    }

/*  calculate the inverse of the first block  */

    for (i = 0; i < pmo.mxllda_cond[0] * pmo.mxlocc_cond[0]; i++)
    {
        Hii[i].r = H_tri[ i].r;
        Hii[i].i = H_tri[ i].i;
    }

    desca = &pmo.desc_cond[0];
    get_inverse_block_p (Hii, Gii, ipiv, desca);

    n1 = pmo.mxllda_cond[0] * pmo.mxlocc_cond[0];
    zcopy (&n1, Gii, &ione, &G_tem[n_begin1[0]], &ione);

/*  iterate to get one more block  */

    for (i = 0; i < N - 1; i++)
    {
        /* get the interaction  Hi,i+1  from input H_tri 
         * Hii1 is a pointer only
         */
        Hii1 = &H_tri[pmo.offdiag_begin[i] ];

        /* Hii now has the matrix Hi+1,i+1  */

        for (j = 0; j < pmo.mxllda_cond[i + 1] * pmo.mxlocc_cond[i + 1]; j++)
        {
            Hii[j].r = H_tri[j + pmo.diag_begin[i + 1]].r;
            Hii[j].i = H_tri[j + pmo.diag_begin[i + 1]].i;
        }


        /* calculate Hi+1,i+1 - Hi+1,i * Gii^0 * Hi,i+1  */

        n1 = ni[i + 1];
        n2 = ni[i];

        desca = &pmo.desc_cond[ ( i   + (i+1) * ct.num_blocks) * DLEN ];
        descb = &pmo.desc_cond[ ( i   +  i    * ct.num_blocks) * DLEN ];
        descc = &pmo.desc_cond[ ( i+1 +  i    * ct.num_blocks) * DLEN ];
        descd = &pmo.desc_cond[ ( i+1 + (i+1) * ct.num_blocks) * DLEN ];

        PZGEMM ("T", "N", &n1, &n2, &n2, &one, Hii1, &ione, &ione, desca, 
               Gii, &ione, &ione, descb, &zero, temp, &ione, &ione, descc);
        PZGEMM ("N", "N", &n1, &n1, &n2, &mone, temp, &ione, &ione, descc, 
               Hii1, &ione, &ione, desca, &one, Hii, &ione, &ione, descd);


        /* now Hii store the matrix Hi+1,i+1 - Hi+1,i * Gii^0 * Hi,i+1
         * Gi+1,i+1, stored in Gii, = Hii^(-1)
         */

        get_inverse_block_p (Hii, Gii, ipiv, descd);

        n1 = pmo.mxllda_cond[i + 1] * pmo.mxlocc_cond[i + 1];
        zcopy (&n1, Gii, &ione, &G_tem[n_begin1[i + 1]], &ione);


        /*  temp = - Hi,i+1 * Gi+1,i+1  */
        n1 = ni[i];
        n2 = ni[i + 1];
        desca = &pmo.desc_cond[ ( i   + (i+1) * ct.num_blocks) * DLEN ];
        descb = &pmo.desc_cond[ ( i+1 + (i+1) * ct.num_blocks) * DLEN ];
        PZGEMM ("N", "N", &n1, &n2, &n2, &mone, Hii1, &ione, &ione, desca, 
                Gii, &ione, &ione, descb, &zero, temp, &ione, &ione, desca);


        /* Gj, i+1 = G0_j,i * temp , j = i, n-1 */
        for (j = 0; j < i + 1; j++)
        {

            n1 = ni[j];
            n2 = ni[i];
            n3 = ni[i + 1];

            desca = &pmo.desc_cond[ ( j   +  i    * ct.num_blocks) * DLEN ];
            descb = &pmo.desc_cond[ ( i   + (i+1) * ct.num_blocks) * DLEN ];
            descc = &pmo.desc_cond[ ( j   + (i+1) * ct.num_blocks) * DLEN ];

            PZGEMM ("N", "N", &n1, &n3, &n2, &one, &G_tem[n_begin1[j]], &ione, &ione, desca, 
                    temp, &ione, &ione, descb, &zero, Hii, &ione, &ione, descc);

            n1 = pmo.mxllda_cond[j] * pmo.mxlocc_cond[i + 1];
            zcopy (&n1, Hii, &ione, &G_tem[n_begin1[j]], &ione);
        }

    }                           /* end  for(i = 0; i < N-1; i++) */


    my_free(n_begin1);
    my_free(ipiv);
    my_free( Hii );
    my_free( Gii );
    my_free( temp );
}
