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
#include "pmo.h"
//#include <cuComplex.h>


void matrix_inverse_anyprobe_cuda (complex double * H_tri, int N, int * ni, int iprobe, 
    complex double * Green_C_row, complex double * Green_C_col)
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

#if GPU_ENABLED
    cuDoubleComplex cumone, cuone,  cuzero, *Hii1, *Hi1i;
    cublasOperation_t transN = CUBLAS_OP_N;

    int nmax, i, j, n1, n2, n3, n4, m;
    int *ipiv, *n_begin1;
    complex double *Hii, *Gii, *temp, *Imatrix;
    complex double *Hmm1, *temp2, *temp3;
    complex double mone, one, zero;
    int ione = 1, ntot, k, maxrow, maxcol; 
    int *desca, *descb, *descc, *descd;

    mone = -1.0;
    one = 1.0;
    zero = 0.0;


    cuone.x = 1.0;
    cumone.x = -1.0;
    cuzero.x = 0.0;
    cuone.y = 0.0;
    cumone.y = 0.0;
    cuzero.y = 0.0;



    /*  find the maximum dimension of the blocks  */


    int tot_row = 0;
    ntot = pmo.ntot;
    maxrow = 0;
    maxcol = 0;
    for (i = 0; i < N; i++)
    {
        tot_row += ni[i];
        maxrow = max(maxrow, pmo.mxllda_cond[i]);
        maxcol = max(maxcol, pmo.mxlocc_cond[i]);
    }

    my_malloc( n_begin1, N, int );
    my_malloc( ipiv, maxrow + pmo.mblock, int ); 

    my_malloc_init( Imatrix, maxrow * maxcol, complex double );
    my_malloc_init( Hii, maxrow * maxcol, complex double );
    my_malloc_init( Gii, maxrow * maxcol, complex double );
    my_malloc_init( temp, maxrow * maxcol, complex double );

    for (i = 0; i < maxrow * maxcol; i++)
    {
        Imatrix[i] = 0.0;
    }

    for (i = 0; i < maxrow; i++)
    {
        Imatrix[i * maxrow + i] = 1.0;
    }


    /*  n_begin: starting address of each diagonal block in H_tri and G_tri
     *  the Hi,i+1 block will start at n_begin[i] + ni[i] * ni[i]
     */

    n_begin1[0] = 0;
    for (i = 1; i < N; i++)
    {
        n_begin1[i] = n_begin1[i - 1] + pmo.mxllda_cond[i - 1] * maxcol;
    }

    /*  find the block index (corresponding to the probe) we are interested in  */
    m = cei.probe_in_block[iprobe - 1];

    /* printf (" while matrix diag: iprobe, block =  %d %d \n", iprobe, m); */

    /* ========================== Part-I ==================================== */

    cublasSetVector( maxrow * maxcol, sizeof( complex double ), Imatrix, ione, ct.gpu_Imatrix, ione );
    cublasSetVector( pmo.ntot_low, sizeof( complex double ), H_tri, ione, ct.gpu_Htri, ione );


    int info;

    n1 = ni[0];
    n2 = n1 * n1;
    cublasZcopy (ct.cublas_handle, n2, ct.gpu_Htri, ione, ct.gpu_Hii, ione);

    cublasZgemm (ct.cublas_handle, transN, transN, n1, n1, n1, &cuone, ct.gpu_Imatrix, maxrow,
            ct.gpu_Imatrix, maxrow, &cuzero, ct.gpu_Gii, n1);
    magma_zgesv_gpu( n1, n1, ct.gpu_Hii, n1, ipiv, ct.gpu_Gii, n1, &info );



    cublasZcopy (ct.cublas_handle, n2, ct.gpu_Gii, ione, &ct.gpu_Grow[n_begin1[0]], ione);
    cublasZcopy (ct.cublas_handle, n2, ct.gpu_Gii, ione, &ct.gpu_Gcol[n_begin1[0]], ione);

    /*  calculate the inverse of the first block  */


    /*  iterate to get one more block  */
    for (i = 0; i < m; i++)
    {
        /* get the interaction  Hi,i+1  from input H_tri 
         * Hii1 is a pointer only
         */
        Hii1 = &ct.gpu_Htri[pmo.offdiag_begin[i] ];
        Hi1i = &ct.gpu_Htri[pmo.lowoffdiag_begin[i] ];

        n2 = ni[i+1] * ni[i+1];
        cublasZcopy (ct.cublas_handle, n2, &ct.gpu_Htri[pmo.diag_begin[i+1]], ione, ct.gpu_Hii, ione);


        /* calculate Hi+1,i+1 - Hi+1,i * Gii^0 * Hi,i+1  */

        n1 = ni[i + 1];
        n2 = ni[i];


        /* calculate Hi+1,i * Gii^0 = tempi+1,i */
        cublasZgemm (ct.cublas_handle, transN, transN, n1, n2, n2, &cuone, Hi1i, n1,
                ct.gpu_Gii, n2, &cuzero, ct.gpu_temp, n1);
        cublasZgemm (ct.cublas_handle, transN, transN, n1, n1, n2, &cumone, ct.gpu_temp, n1,
                Hii1, n2, &cuone, ct.gpu_Hii, n1);

        cublasZgemm (ct.cublas_handle, transN, transN, n1, n1, n1, &cuone, ct.gpu_Imatrix, maxrow,
                ct.gpu_Imatrix, maxrow, &cuzero, ct.gpu_Gii, n1);
        magma_zgesv_gpu( n1, n1, ct.gpu_Hii, n1, ipiv, ct.gpu_Gii, n1, &info );
        n2 = ni[i+1] * ni[i+1];
        cublasZcopy (ct.cublas_handle, n2, ct.gpu_Gii, ione, &ct.gpu_Grow[n_begin1[i+1]], ione);
        cublasZcopy (ct.cublas_handle, n2, ct.gpu_Gii, ione, &ct.gpu_Gcol[n_begin1[i+1]], ione);

        n1 = ni[i];
        n2 = ni[i + 1];
        cublasZgemm (ct.cublas_handle, transN, transN, n1, n2, n2, &cumone, Hii1, n1,
                ct.gpu_Gii, n2, &cuzero, ct.gpu_temp, n1);

        /* Gj, i+1 = G0_j,i * tempi,i+1 ; j = 0, i   eq.15*/
        for (j = 0; j <= i; j++)
        {

            n3 = ni[j];

            cublasZgemm (ct.cublas_handle, transN, transN, n3, n2, n1, &cuone, &ct.gpu_Grow[n_begin1[j]], n3,
                    ct.gpu_temp, n1, &cuzero, ct.gpu_Hii, n3);
            n4 = ni[j] * ni[i+1];
            cublasZcopy (ct.cublas_handle, n4, ct.gpu_Hii, ione, &ct.gpu_Grow[n_begin1[j]], ione);

        }                           /* end  for(i = 0; i < m; i++) */

        /* Gi+1, j = G_i+1,i+1 * Hi+1,i * G0_i,j; j = 0, i  eq. 16*/

        cublasZgemm (ct.cublas_handle, transN, transN, n2, n1, n2, &cumone, ct.gpu_Gii, n2,
                Hi1i, n2, &cuzero, ct.gpu_temp, n2);

        for (j = 0; j <= i; j++)
        {

            n3 = ni[j];

            cublasZgemm (ct.cublas_handle, transN, transN, n2, n3, n1, &cuone, ct.gpu_temp, n2, 
                    &ct.gpu_Gcol[n_begin1[j]], n1, &cuzero, ct.gpu_Hii, n2);
            n4 = ni[j] * ni[i+1];
            cublasZcopy (ct.cublas_handle, n4, ct.gpu_Hii, ione, &ct.gpu_Gcol[n_begin1[j]], ione);

        }                           /* end  for(i = 0; i < m; i++) */

    }


    /* ========================== Part-II starts here ============================ */

    if(m < N - 1 )  
    {

        /*  calculate the inverse of the last block  */
        n2 = ni[N-1] * ni[N-1];
        cublasZcopy (ct.cublas_handle, n2, &ct.gpu_Htri[pmo.diag_begin[N-1]], ione, ct.gpu_Hii, ione);

        n1 = ni[N-1];
        cublasZgemm (ct.cublas_handle, transN, transN, n1, n1, n1, &cuone, ct.gpu_Imatrix, maxrow,
                ct.gpu_Imatrix, maxrow, &cuzero, ct.gpu_Gii, n1);
        magma_zgesv_gpu( n1, n1, ct.gpu_Hii, n1, ipiv, ct.gpu_Gii, n1, &info );

        n2 = ni[N-1] * ni[N-1];
        cublasZcopy (ct.cublas_handle, n2, ct.gpu_Gii, ione, &ct.gpu_Grow[n_begin1[N-1]], ione);
        cublasZcopy (ct.cublas_handle, n2, ct.gpu_Gii, ione, &ct.gpu_Gcol[n_begin1[N-1]], ione);





        /*  iterate to get one more block  */

        for (i = N - 1; i > m + 1 ; i--) 
            /*    for (i = N - 1; i > m; i--)   */ 
        {
            /* get the interaction  Hi-1,i  from input H_tri 
             * Hii1 is a pointer only
             */
            Hii1 = &ct.gpu_Htri[pmo.offdiag_begin[i-1] ];
            Hi1i = &ct.gpu_Htri[pmo.lowoffdiag_begin[i-1] ];

            /* Hii now has the matrix Hi+1,i+1  */

            n2 = ni[i-1] * ni[i-1];
            cublasZcopy (ct.cublas_handle, n2, &ct.gpu_Htri[pmo.diag_begin[i-1]], ione, ct.gpu_Hii, ione);


            /* calculate Hi-1,i-1 - Hi-1,i * Gii^0 * Hi,i-1  */

            n1 = ni[i - 1];
            n2 = ni[i];

            /* calculate Hi-1,i * Gii^0 = tempi-1,i */
            cublasZgemm (ct.cublas_handle, transN, transN, n1, n2, n2, &cuone, Hii1, n1,
                    ct.gpu_Gii, n2, &cuzero, ct.gpu_temp, n1);

            /* calculate Hi-1,i-1 - tempi-1,i * Hi,i-1 = Hi-1,i-1 */
            cublasZgemm (ct.cublas_handle, transN, transN, n1, n1, n2, &cumone, ct.gpu_temp, n1,
                    Hi1i, n2, &cuone, ct.gpu_Hii, n1);

            /* now Hii store the matrix Hi-1,i-1 - Hi-1,i * Gii^0 * Hi,i-1
             * Gi-1,i-1, stored in Gii, = Hii^(-1)
             */
            cublasZgemm (ct.cublas_handle, transN, transN, n1, n1, n1, &cuone, ct.gpu_Imatrix, maxrow,
                    ct.gpu_Imatrix, maxrow, &cuzero, ct.gpu_Gii, n1);
            magma_zgesv_gpu( n1, n1, ct.gpu_Hii, n1, ipiv, ct.gpu_Gii, n1, &info );

            n2 = ni[i-1] * ni[i-1];
            cublasZcopy (ct.cublas_handle, n2, ct.gpu_Gii, ione, &ct.gpu_Grow[n_begin1[i-1]], ione);
            cublasZcopy (ct.cublas_handle, n2, ct.gpu_Gii, ione, &ct.gpu_Gcol[n_begin1[i-1]], ione);


            /*  temp[i,i-1] = - Hi,i-1 * Gi-1,i-1  */
            n1 = ni[i];
            n2 = ni[i - 1];

            cublasZgemm (ct.cublas_handle, transN, transN, n1, n2, n2, &cumone, Hi1i, n1,
                    ct.gpu_Gii, n2, &cuzero, ct.gpu_temp, n1);

            /* Gj, i-1 = G0_j,i * temp , j = i, n-1  eq. 18*/
            for (j = i; j < N; j++)
            {

                n3 = ni[j];

                cublasZgemm (ct.cublas_handle, transN, transN, n3, n2, n1, &cuone, &ct.gpu_Grow[n_begin1[j]], n3,
                        ct.gpu_temp, n1, &cuzero, ct.gpu_Hii, n3);


                n4 = ni[j] * ni[i-1];
                cublasZcopy (ct.cublas_handle, n4, ct.gpu_Hii, ione, &ct.gpu_Grow[n_begin1[j]], ione);
            }


            /* Gi-1, j, = -G_i-1,i-1 * Hi-1,i * G0_i,j , j = i, n-1  eq. 19*/
            cublasZgemm (ct.cublas_handle, transN, transN, n2, n1, n2, &cumone, 
                    ct.gpu_Gii, n2, Hii1, n2, &cuzero, ct.gpu_temp, n2);

            for (j = i; j < N; j++)
            {

                n3 = ni[j];

                cublasZgemm (ct.cublas_handle, transN, transN, n2, n3, n1, &cuone, ct.gpu_temp, n2, 
                        &ct.gpu_Gcol[n_begin1[j]], n1, &cuzero, ct.gpu_Hii, n2);


                n4 = ni[j] * ni[i-1];
                cublasZcopy (ct.cublas_handle, n4, ct.gpu_Hii, ione, &ct.gpu_Gcol[n_begin1[j]], ione);
            }

        }                           /* end  (i = N - 1; i > m + 1 ; i--) */


        /* ========================== Part-III starts here ========================== */
        /*        printf (" value of m ....  %d \n", m);        */


        /* get the interaction Hm,m+1 from input H_tri */
        Hii1 = &ct.gpu_Htri[pmo.offdiag_begin[m] ];
        Hi1i = &ct.gpu_Htri[pmo.lowoffdiag_begin[m] ];


        n1 = ni[m];
        n2 = ni[m + 1];

        /* eq. 25, the inverse part */
        /*  ct.gpu_temp = Hm+1,m *  G0mm  */
        cublasZgemm (ct.cublas_handle, transN, transN, n2, n1, n1,
                &cuone, Hi1i, n2, &ct.gpu_Grow[n_begin1[m]], n1,
                &cuzero, ct.gpu_temp, n2);
        /*  ct.gpu_Gii = Hm,m+1 * G0 m+1,m+1 */
        cublasZgemm (ct.cublas_handle, transN, transN, n1, n2, n2, &cuone, Hii1, n1,
                &ct.gpu_Grow[n_begin1[m+1]], n2, &cuzero, ct.gpu_Gii, n1);

        /*  ct.gpu_Hii = 1 - Hm,m+1 * Gm+1,m+1 * Hm+1,m * Gmm */
        cublasZgemm (ct.cublas_handle, transN, transN, n1, n1, n1, &cuone, ct.gpu_Imatrix, maxrow,
                ct.gpu_Imatrix, maxrow, &cuzero, ct.gpu_Hii, n1);
        cublasZgemm (ct.cublas_handle, transN, transN, n1, n1, n2, &cumone, ct.gpu_Gii, n1,
                ct.gpu_temp, n2, &cuone, ct.gpu_Hii, n1);

        /* set Gii to be unitary matrix */
        cublasZgemm (ct.cublas_handle, transN, transN, n1, n1, n1, &cuone, ct.gpu_Imatrix, maxrow,
                ct.gpu_Imatrix, maxrow, &cuzero, ct.gpu_Gii, n1);

        /* get the inverse of ct.gpu_Hii , see second bracket in eq. 25 */
        magma_zgesv_gpu( n1, n1, ct.gpu_Hii, n1, ipiv, ct.gpu_Gii, n1, &info );

        /* gnu_Hii = Hm+1,m * Gmm * (... )^(-1)   eq.25 */
        cublasZgemm (ct.cublas_handle, transN, transN, n2, n1, n1, &cuone, ct.gpu_temp, n2,
                ct.gpu_Gii, n1, &cuzero, ct.gpu_Hii, n2);

        /* eq. 25, first part for j <= m */
        for (j = 0; j <= m; j++)
        {
            n3 = ni[j];


            cublasZgemm (ct.cublas_handle, transN, transN, n3, n1, n1, &cuone, 
                    &ct.gpu_Grow[n_begin1[j]], n3, ct.gpu_Gii, n1, &cuzero, ct.gpu_temp, n3);

            n4 = ni[j] * ni[m];
            cublasZcopy (ct.cublas_handle, n4, ct.gpu_temp, ione, &ct.gpu_Grow[n_begin1[j]], ione);
        }

        for (j = m + 1; j < N; j++)
        {
            n3 = ni[j];

            cublasZgemm (ct.cublas_handle, transN, transN, n3, n1, n2, &cumone, &ct.gpu_Grow[n_begin1[j]], n3,
                    ct.gpu_Hii, n2, &cuzero, ct.gpu_temp, n3);

            n4 = ni[j] * ni[m];
            cublasZcopy (ct.cublas_handle, n4, ct.gpu_temp, ione, &ct.gpu_Grow[n_begin1[j]], ione);

        }

        /* eq. 24  */
        /* calculate Gmm * Hm,m+1 = temp[m,m+1] */

        cublasZgemm (ct.cublas_handle, transN, transN, n1, n2, n1, &cuone, &ct.gpu_Gcol[n_begin1[m]], n1,
                Hii1, n1, &cuzero, ct.gpu_temp, n1);



        /* calculate  temp[m,m+1] * Gm+1,m+1 = temp2[m,m+1] */
        cublasZgemm (ct.cublas_handle, transN, transN, n1, n2, n2, &cuone, ct.gpu_temp, n1,
                &ct.gpu_Gcol[n_begin1[m+1]], n2, &cuzero, ct.gpu_Gii, n1);
        cublasZgemm (ct.cublas_handle, transN, transN, n1, n1, n1, &cuone, ct.gpu_Imatrix, maxrow,
                ct.gpu_Imatrix, maxrow, &cuzero, ct.gpu_Hii, n1);
        cublasZgemm (ct.cublas_handle, transN, transN, n1, n1, n2, &cumone, ct.gpu_Gii, n1,
                Hi1i, n2, &cuone, ct.gpu_Hii, n1);
        /* Now 1 - Gmm * Hm,m+1 * Gm+1,m+1 * Hm+1,m is stored in Hii[m,m] */

        /* calculate identity[m,m] - temp2[m,m+1] * Hm+1,m = identity[m,m] */

        cublasZgemm (ct.cublas_handle, transN, transN, n1, n1, n1, &cuone, ct.gpu_Imatrix, maxrow,
                ct.gpu_Imatrix, maxrow, &cuzero, ct.gpu_Gii, n1);



        /* get the inverse of identity[m,m] and stored in temp2[m,m]           */
        magma_zgesv_gpu( n1, n1, ct.gpu_Hii, n1, ipiv, ct.gpu_Gii, n1, &info );


        /* ct.gpu_Hii store:  (1-G0mm * Hm,m+1 * G0m+1,m+1 * Hm+1,m)^-1 * G0mm * Hm,m+1 */
        cublasZgemm (ct.cublas_handle, transN, transN, n1, n2, n1, &cuone, ct.gpu_Gii, n1,
                ct.gpu_temp, n1, &cuzero, ct.gpu_Hii, n1);

        for (j = 0; j <= m; j++)
        {
            n3 = ni[j];


            cublasZgemm (ct.cublas_handle, transN, transN, n1, n3, n1, &cuone, ct.gpu_Gii, n1, 
                    &ct.gpu_Gcol[n_begin1[j]], n1, &cuzero, ct.gpu_temp, n1);

            n4 = ni[j] * ni[m];
            cublasZcopy (ct.cublas_handle, n4, ct.gpu_temp, ione, &ct.gpu_Gcol[n_begin1[j]], ione);
        }


        for (j = m + 1; j < N; j++)
        {
            n3 = ni[j];

            cublasZgemm (ct.cublas_handle, transN, transN, n1, n3, n2, &cumone, ct.gpu_Hii, n1, 
                    &ct.gpu_Gcol[n_begin1[j]], n2, &cuzero, ct.gpu_temp, n1);

            n4 = ni[j] * ni[m];
            cublasZcopy (ct.cublas_handle, n4, ct.gpu_temp, ione, &ct.gpu_Gcol[n_begin1[j]], ione);

        }



    }                           /* if statement ends here */
    /* ================================================================== */


    n4 = tot_row * maxcol;
    cublasGetVector(n4, sizeof( complex double ), ct.gpu_Grow, ione, Green_C_row, ione );
    n4 = tot_row * maxrow;
    cublasGetVector(n4, sizeof( complex double ), ct.gpu_Gcol, ione, Green_C_col, ione );



    my_free(n_begin1);
    my_free(ipiv);
    my_free( Hii );
    my_free( Imatrix );
    my_free( Gii );
    my_free( temp );
#endif
}



