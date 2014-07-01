/************************** SVN Revision Information **************************
 **    $Id: matrix_inverse_p.c 1348 2011-04-28 18:32:01Z luw $    **
******************************************************************************/
 
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "main.h"
#include "init_var.h"
#include "LCR.h"
#include "pmo.h"


void matrix_inverse_cuda (complex double * H_tri, complex double * G_tri)
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

#if GPU_ENABLED
    int  i, j, n1, n2, n3, n4, n5, n6, n7, n8;
    int *ipiv;
    complex double *Hii, *Gii, *G_tem, *temp, *Imatrix;
    complex double mone, one, zero;
    int ione = 1, *ntem_begin;
    int idx;
    int ntot_row;
    int maxrow, maxcol;
    int *ni, *diag_begin, *offdiag_begin, N;
    cuDoubleComplex cuhalf, cumone, cuone,  cuzero, *Hlower, *Hupper, *Gii0;
    int size;

    cublasOperation_t transT = CUBLAS_OP_T, transN = CUBLAS_OP_N;

    int info;
    ni = ct.block_dim;
    N = ct.num_blocks;
    offdiag_begin = pmo.offdiag_begin ;
    diag_begin = pmo.diag_begin ;
    mone = -1.0;
    one = 1.0;
    zero = 0.0; 
    cuone.x = 1.0;
    cumone.x = -1.0;
    cuzero.x = 0.0;
    cuone.y = 0.0;
    cumone.y = 0.0;
    cuzero.y = 0.0;

    cuhalf.x = 0.5;
    cuhalf.y = 0.0;
    /*  find the maximum dimension of the blocks  */


    ntot_row = 0;
    maxrow=0;
    maxcol=0;
    for (i = 0; i < ct.num_blocks; i++)
    {
        ntot_row += ni[i];
        maxrow = max(maxrow, ni[i]);
        maxcol = max(maxcol, ni[i]);

    }

    my_malloc_init( ntem_begin, ct.num_blocks, int);
    size = maxrow + pmo.mblock;
    my_malloc_init( ipiv, size, int );
    size = maxrow * maxcol;
    my_malloc_init( Hii, size, complex double );
    my_malloc_init( Gii, size, complex double );
    my_malloc_init( temp, size, complex double );
    my_malloc_init( Imatrix, size, complex double );
    size = ntot_row * maxcol;
    my_malloc_init( G_tem, size, complex double );	


    /*
     *  ntem_begin: starting address of one column of G_tem for each block
     */

    ntem_begin[0] = 0;
    for (i = 1; i < ct.num_blocks; i++)
    {
        ntem_begin[i] = ntem_begin[i - 1] + pmo.mxllda_cond[i - 1] * maxcol;
    }

    /*  calculate the inverse of the first block  */


    for (i = 0; i < maxrow * maxcol; i++)
    {
        Imatrix[i] = 0.0;
    }

    for (i = 0; i < maxrow; i++)
    {
        Imatrix[i * maxrow + i] = 1.0;
    }


    n1 = ni[0];
    n2 = ni[0] * ni[0];

    cublasSetVector( maxrow * maxcol, sizeof( complex double ), Imatrix, ione, ct.gpu_Imatrix, ione );
    cublasSetVector( pmo.ntot_low, sizeof( complex double ), H_tri, ione, ct.gpu_Htri, ione );

    cublasZgemm (ct.cublas_handle, transN, transN, n1, n1, n1, &cuone, ct.gpu_Imatrix, maxrow,
            ct.gpu_Imatrix, maxrow, &cuzero, ct.gpu_Gii, n1);
    cublasZcopy (ct.cublas_handle, n2, ct.gpu_Htri, ione, ct.gpu_Hii, ione);

    magma_zgesv_gpu( n1, n1, ct.gpu_Hii, n1, ipiv, ct.gpu_Gii, n1, &info );

    cublasZcopy (ct.cublas_handle, n2, ct.gpu_Gii, ione, ct.gpu_Gtri, ione);

    //    cublasGetVector( n2, sizeof( complex double ), ct.gpu_Gtri, ione, G_tri, ione );

    //	zgetrf( &n1, &n1, Hii, &n1, ipiv, &info);
    //	zgetrs(&transn, &n1, &n1, Hii, &n1, ipiv, Gii, &n1, &info);
    //	zgesv(&n1, &n1, Hii, &n1, ipiv, Gii, &n1, &info);

    /*  iterate to get one more block  */


    for (i = 0; i < N - 1; i++)
    {
        /* get the interaction  Hi,i+1  from input H_tri 
         * Hupper is a pointer only
         */
        Hupper = &ct.gpu_Htri[offdiag_begin[i] ];
        Gii0 = &ct.gpu_Gtri[diag_begin[i] ];
        Hlower = &ct.gpu_Htri[pmo.lowoffdiag_begin[i] ];

        /* Hii now has the matrix Hi+1,i+1  */

        n2 = ni[i+1] * ni[i+1];
        cublasZcopy (ct.cublas_handle, n2, &ct.gpu_Htri[diag_begin[i+1]], ione, ct.gpu_Hii, ione);
        /*
           for (j = 0; j < ni[i+1] * ni[i + 1]; j++)
           {
           Hii[j] = H_tri[j + diag_begin[i+ 1]];
           }
         */

        /* calculate Hi+1,i+1 - Hi+1,i * Gii^0 * Hi,i+1  */

        n1 = ni[i + 1];
        n2 = ni[i];


        cublasZgemm (ct.cublas_handle, transN, transN, n1, n2, n2, &cuone, Hlower, n1,
                Gii0, n2, &cuzero, ct.gpu_temp, n1);
        cublasZgemm (ct.cublas_handle, transN, transN, n1, n1, n2, &cumone, ct.gpu_temp, n1,
                Hupper, n2, &cuone, ct.gpu_Hii, n1);


        /* now Hii store the matrix Hi+1,i+1 - Hi+1,i * Gii^0 * Hi,i+1
         * Gi+1,i+1, stored in Gii, = Hii^(-1)
         */

        cublasZgemm (ct.cublas_handle, transN, transN, n1, n1, n1, &cuone, ct.gpu_Imatrix, maxrow,
                ct.gpu_Imatrix, maxrow, &cuzero, ct.gpu_Gii, n1);
        magma_zgesv_gpu( n1, n1, ct.gpu_Hii, n1, ipiv, ct.gpu_Gii, n1, &info );

        n1 = ni[i + 1] * ni[i + 1];
        cublasZcopy (ct.cublas_handle, n1, ct.gpu_Gii, ione, &ct.gpu_Gtri[diag_begin[i+1]], ione);
        //	    zcopy (&n1, Gii, &ione, &G_tri[diag_begin[i + 1]], &ione);

        /*  Gi,i+1 =  Gii^0 * Hi,i+1 * Gi+1,i+1  */

        n1 = ni[i];
        n2 = ni[i + 1];
        n3 = offdiag_begin[i];


        /* temp = Hi+1,i * Gii^0  eq. 11 for j = I term*/
        cublasZgemm (ct.cublas_handle, transN, transN, n2, n1, n1, &cumone, Hlower, n2, Gii0, n1,
                &cuzero, ct.gpu_temp, n2);
        cublasZgemm (ct.cublas_handle, transN, transN, n2, n1, n2, &cuone, &ct.gpu_Gtri[diag_begin[i+1]], n2, ct.gpu_temp, n2,
                &cuzero, &ct.gpu_Gtri[pmo.lowoffdiag_begin[i]], n2);

        n4 = ni[i] * ni[i + 1];
        cublasZcopy (ct.cublas_handle, n4, &ct.gpu_Gtri[pmo.lowoffdiag_begin[i]], ione, &ct.gpu_Gcol[ntem_begin[i]], ione);


        /* temp = Gii^0 * Hi,i+1 */
        cublasZgemm (ct.cublas_handle, transN, transN, n1, n2, n1, &cumone, Gii0, n1,
                Hupper, n1, &cuzero, ct.gpu_temp, n1);
        //	    zgemm ("N", "N", &n1, &n2, &n1, &mone, Gii0, &n1, 
        //			    Hupper, &n1, &zero, temp, &n1);

        /* G(i,i+1) = temp * G(i+1,i+1)  also == Grow(i,i+1) eq. 8 for j = I */
        cublasZgemm (ct.cublas_handle, transN, transN, n1, n2, n2, &cuone, ct.gpu_temp, n1,
                ct.gpu_Gii, n2, &cuzero, &ct.gpu_Gtri[n3], n1);
        //	    zgemm ("N", "N", &n1, &n2, &n2, &one, temp, &n1,
        //			    Gii, &n2, &zero, &G_tri[n3], &n1);

        n4 = ni[i] * ni[i + 1];
        cublasZcopy (ct.cublas_handle, n4, &ct.gpu_Gtri[n3], ione, &ct.gpu_Grow[ntem_begin[i]], ione);
        //	    zcopy (&n4, &G_tri[n3], &ione, &G_tem[ntem_begin[i]], &ione);

        /* update Gii eq. 9 for j=I */

        cublasZgemm (ct.cublas_handle, transN, transN, n1, n1, n2, &cuone, ct.gpu_temp, n1,
                &ct.gpu_Gtri[pmo.lowoffdiag_begin[i]], n2, &cuone, &ct.gpu_Gtri[diag_begin[i]], n1);

        for (j = i - 1; j >= 0; j--)
        {

            n1 = ni[j];         /* dimension of block j */
            n2 = ni[i];         /* dimension of block i */
            n3 = ni[i + 1];     /* dimension of block i+1 */
            n4 = ntem_begin[j]; /* starting address of Gtem(j,*) in G_tem */
            n5 = diag_begin[j];    /* starting address of G(j,j) in G_tri */
            n6 = offdiag_begin[j];    /* starting address of G(j,j+1) in G_tri */
            n7 = ntem_begin[j + 1];     /* starting address of Gtem(j+1,*) in G_tem */
            n8 = ni[j + 1];     /* dimension of block j+1 */


            /* gpu_Gii = -Hi+1,i * G0(i,j)   */

            cublasZgemm (ct.cublas_handle, transN, transN, n3, n1, n2, &cumone, Hlower, n3, &ct.gpu_Gcol[n4], n2,
                    &cuzero, ct.gpu_Gii, n3);

            /*  G(I+1, j) = G(I+1, I+1) * gpu_Gii  eq. 11 */
            cublasZgemm (ct.cublas_handle, transN, transN, n3, n1, n3, &cuone, 
                    &ct.gpu_Gtri[diag_begin[i+1]], n3, ct.gpu_Gii, n3, &cuzero, &ct.gpu_Gcol[n4], n3);

            /* temp = -G0(j,i) * Hi,i+1  */

            cublasZgemm (ct.cublas_handle, transN, transN, n1, n3, n2, &cumone, &ct.gpu_Grow[n4], n1,
                    Hupper, n2, &cuzero, ct.gpu_temp, n1);

            /* G0(j, i+1) = temp * G(i+1,i+1) eq. 8 */

            cublasZgemm (ct.cublas_handle, transN, transN, n1, n3, n3, &cuone, ct.gpu_temp, n1,
                    &ct.gpu_Gtri[diag_begin[i+1]], n3, &cuzero, &ct.gpu_Grow[n4], n1);

            /* G(j,j) = G0(j,j) + temp * G(i+1,j) eq. 9*/
            cublasZgemm (ct.cublas_handle, transN, transN, n1, n1, n3, &cuone, ct.gpu_temp, n1,
                    &ct.gpu_Gcol[n4], n3, &cuone, &ct.gpu_Gtri[n5], n1);

            /* G(j,j+1) = G0(j,j+1) + temp * G(i+1,j+1)  eq. 10 */
            cublasZgemm (ct.cublas_handle, transN, transN, n1, n8, n3, &cuone, ct.gpu_temp, n1,
                    &ct.gpu_Gcol[n7], n3, &cuone, &ct.gpu_Gtri[n6], n1);

            /* G(j+1,j) = G0(j+1,j) + G(j+1,I+1) * gpu_Gii  eq. 10 */
            cublasZgemm (ct.cublas_handle, transN, transN, n8, n1, n3, &cuone, &ct.gpu_Grow[n7], n8,
                    ct.gpu_Gii, n3, &cuone, &ct.gpu_Gtri[pmo.lowoffdiag_begin[j]], n8);
        }                       /* end for (j--) */

    }                           /* end  for(i = 0; i < N-1; i++) */

// due to time-inverse symmetry, H and S are Hermitian,  G(-k) = G(k)^~
//we average the lower and upper part of G(k) and only copy the upper triangle parts back to cpu

    for(i = 0; i < ct.num_blocks - 1; i++)
    {
        n1 = ct.block_dim[i];
        n2 = ct.block_dim[i+1];
        n3 = pmo.offdiag_begin[i];
        n4 = pmo.lowoffdiag_begin[i];

        cublasZgemm (ct.cublas_handle, transT, transN, n1, n2, n2, &cuhalf, &ct.gpu_Gtri[n4], n2,
                    ct.gpu_Imatrix, maxrow, &cuhalf, &ct.gpu_Gtri[n3], n1);
    }
        


    cublasGetVector( pmo.ntot, sizeof( complex double ), ct.gpu_Gtri, ione, G_tri, ione );

    my_free( ntem_begin );
    my_free( ipiv );
    my_free( Hii );
    my_free( Gii );
    my_free( temp );
    my_free( G_tem );
    /*
       acc_free( ipiv );
       acc_free( Hii );
       acc_free( Gii );
       acc_free( temp );
       acc_free( G_tem );
     */
#endif
}
