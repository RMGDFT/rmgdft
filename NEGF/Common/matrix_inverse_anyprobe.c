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

void *memory_ptr_host_device(void *ptr_host, void *ptr_device);
void matrix_inverse_driver(double *, int *);

void matrix_inverse_anyprobe (complex double * H_tri_host, int N, int * ni, int iprobe, 
        complex double * G_row_host, complex double *G_col_host)
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

    int nmax, i, j, n1, n2, n3, n4, m;
    int *n_begin1;
    complex double *H_tri, *Hii, *Gii, *temp, *Hii1, *Hi1i;
    complex double *Hii_host, *Gii_host, *temp_host;
    complex double *Imatrix, *Imatrix_host;
    complex double *Grow, *Gcol;
    complex double mone, one, zero;
    int ione = 1, ntot, k, maxrow, maxcol, totrow, totcol; 
    int *desca, *descb, *descc, *descd;
    int *desce, *descf, nz;

    mone = -1.0;
    one = 1.0;
    zero = 0.0;

    /*  find the maximum dimension of the blocks  */


    ntot = pmo.ntot;
    maxrow = 0;
    maxcol = 0;
    totcol = 0;
    totrow = 0;
    for (i = 0; i < N; i++)
    {
        maxrow = rmg_max(maxrow, pmo.mxllda_cond[i]);
        maxcol = rmg_max(maxcol, pmo.mxlocc_cond[i]);
        totrow += pmo.mxllda_cond[i];
        totcol += pmo.mxlocc_cond[i];

    }

    my_malloc( n_begin1, N, int );

    my_malloc_init( Imatrix_host, maxrow * maxcol, complex double );
    my_malloc_init( Hii_host, maxrow * maxcol, complex double );
    my_malloc_init( Gii_host, maxrow * maxcol, complex double );
    my_malloc_init( temp_host, maxrow * maxcol, complex double );

    Gii = memory_ptr_host_device(Gii_host, ct.gpu_Gii);
    Imatrix = memory_ptr_host_device(Imatrix_host, ct.gpu_Imatrix);
    Hii = memory_ptr_host_device(Hii_host, ct.gpu_Hii);
    temp = memory_ptr_host_device(temp_host, ct.gpu_temp);
    Grow = memory_ptr_host_device(G_row_host, ct.gpu_Grow);
    Gcol = memory_ptr_host_device(G_col_host, ct.gpu_Gcol);
    H_tri = memory_ptr_host_device(H_tri_host, ct.gpu_Htri);


    setvector_host_device (pmo.ntot_low, sizeof(complex double), H_tri_host, ione, ct.gpu_Htri, ione);

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


    /*  calculate the inverse of the first block  */

    desca = &pmo.desc_cond[0];

    n1 = pmo.mxllda_cond[0] * pmo.mxlocc_cond[0];

    zcopy_driver (n1, H_tri, ione, Gii, ione);

    matrix_inverse_driver(Gii, desca);
    zcopy_driver (n1, Gii, ione, &Grow[n_begin1[0]], ione);
    zcopy_driver (n1, Gii, ione, &Gcol[n_begin1[0]], ione);


    /*  iterate to get one more block  */
    for (i = 0; i < m; i++)
    {
        /* get the interaction  Hi,i+1  from input H_tri 
         * Hii1 is a pointer only
         */
        Hii1 = &H_tri[pmo.offdiag_begin[i] ];
        Hi1i = &H_tri[pmo.lowoffdiag_begin[i] ];
    


        /* calculate Hi+1,i+1 - Hi+1,i * Gii^0 * Hi,i+1  */

        n1 = ni[i + 1];
        n2 = ni[i];

        desca = &pmo.desc_cond[ ( i   + (i+1) * N) * DLEN ];
        descb = &pmo.desc_cond[ ( i   +  i    * N) * DLEN ];
        descc = &pmo.desc_cond[ ( i+1 +  i    * N) * DLEN ];
        descd = &pmo.desc_cond[ ( i+1 + (i+1) * N) * DLEN ];

        /* calculate Hi+1,i * Gii^0 = tempi+1,i */
        zgemm_driver ("N", "N", n1, n2, n2, one, Hi1i, ione, ione, descc, 
                Gii, ione, ione, descb, zero, temp, ione, ione, descc);

        /* Gii now store the matrix Hi+1,i+1  */
        nz = pmo.mxllda_cond[i + 1] * pmo.mxlocc_cond[i + 1]; 
        zcopy_driver (nz, &H_tri[pmo.diag_begin[i + 1]], ione, Gii, ione);

        /* calculate Hi+1,i+1 - tempi+1,i * Hi,i+1 = Hi+1,i+1 */
        zgemm_driver ("N", "N", n1, n1, n2, mone, temp, ione, ione, descc, 
                Hii1, ione, ione, desca, one, Gii, ione, ione, descd);

        /* now Gii store the matrix Hi+1,i+1 - Hi+1,i * Gii^0 * Hi,i+1
         * Gi+1,i+1, stored in Gii, = Hii^(-1)
         */

        matrix_inverse_driver(Gii, descd);

        nz = pmo.mxllda_cond[i + 1] * pmo.mxlocc_cond[i + 1];
        zcopy_driver (nz, Gii, ione, &Grow[n_begin1[i + 1]], ione);
        zcopy_driver (nz, Gii, ione, &Gcol[n_begin1[i + 1]], ione);


        /*  tempi,i+1 = - Hi,i+1 * Gi+1,i+1  */
        n1 = ni[i];
        n2 = ni[i + 1];
        desca = &pmo.desc_cond[ ( i   + (i+1) * N) * DLEN ];
        descb = &pmo.desc_cond[ ( i+1 + (i+1) * N) * DLEN ];

        /* Gj, i+1 = -G0_j,i * Hi,i+1 * Gi+1,i+1 ; j = 0, i   eq.15*/
        zgemm_driver ("N", "N", n1, n2, n2, mone, Hii1, ione, ione, desca, 
                Gii, ione, ione, descb, zero, temp, ione, ione, desca);
        for (j = 0; j <= i; j++)
        {

            n3 = ni[j];

            descc = &pmo.desc_cond[ ( j   +  i    * N) * DLEN ];
            descd = &pmo.desc_cond[ ( j   + (i+1) * N) * DLEN ];

            zgemm_driver ("N", "N", n3, n2, n1, one, &Grow[n_begin1[j]], ione, ione, descc, 
                    temp, ione, ione, desca, zero, Hii, ione, ione, descd);

            nz = pmo.mxllda_cond[j] * pmo.mxlocc_cond[i + 1];
            zcopy_driver (nz, Hii, ione, &Grow[n_begin1[j]], ione);
        }

        /* Gi+1, j = -Gi+1,i+1 * Hi+1,i * G0_i,j ; j = 0, i   eq.15*/
        desca = &pmo.desc_cond[ ( i+1 +  i    * N) * DLEN ];
        descb = &pmo.desc_cond[ ( i+1 + (i+1) * N) * DLEN ];
        zgemm_driver ("N", "N", n2, n1, n2, mone, Gii, ione, ione, descb, 
                Hi1i, ione, ione, desca, zero, temp, ione, ione, desca);
        for (j = 0; j <= i; j++)
        {

            n3 = ni[j];

            descc = &pmo.desc_cond[ ( i   +  j   * N) * DLEN ];
            descd = &pmo.desc_cond[ ( i+1 +  j   * N) * DLEN ];

            zgemm_driver ("N", "N", n2, n3, n1, one, temp, ione, ione, desca, 
                    &Gcol[n_begin1[j]], ione, ione, descc, zero, Hii, ione, ione, descd);

            nz = pmo.mxllda_cond[i+1] * pmo.mxlocc_cond[j];
            zcopy_driver (nz, Hii, ione, &Gcol[n_begin1[j]], ione);
        }


    }                           /* end  for(i = 0; i < m; i++) */



    /* ========================== Part-II starts here ============================ */

    if(m < N - 1 )  
    {

        /*  calculate the inverse of the last block  */

        nz = pmo.mxllda_cond[N - 1] * pmo.mxlocc_cond[N - 1];
        zcopy_driver (nz, &H_tri[pmo.diag_begin[N - 1]], ione, Gii, ione);


        desca = &pmo.desc_cond[ (N-1 + (N-1) * N) * DLEN];
        matrix_inverse_driver (Gii, desca);



        n1 = pmo.mxllda_cond[N - 1] * pmo.mxlocc_cond[N - 1];
        zcopy_driver (n1, Gii, ione, &Grow[n_begin1[N - 1]], ione);
        zcopy_driver (n1, Gii, ione, &Gcol[n_begin1[N - 1]], ione);

        /*  iterate to get one more block  */

        for (i = N - 1; i > m + 1 ; i--) 
            /*    for (i = N - 1; i > m; i--)   */ 
        {
            /* get the interaction  Hi-1,i  from input H_tri 
             * Hii1 is a pointer only
             */
            Hii1 = &H_tri[pmo.offdiag_begin[i-1] ];
            Hi1i = &H_tri[pmo.lowoffdiag_begin[i-1] ];


            n1 = ni[i - 1];
            n2 = ni[i];
            desca = &pmo.desc_cond[ (i-1 +  i    * N ) * DLEN ];
            descb = &pmo.desc_cond[ (i   +  i    * N ) * DLEN ];
            descc = &pmo.desc_cond[ (i-1 + (i-1) * N ) * DLEN ];
            descc = &pmo.desc_cond[ (i   + (i-1) * N ) * DLEN ];

            /* calculate Hi-1,i * Gii^0 = tempi-1,i */
            zgemm_driver ("N", "N", n1, n2, n2, one, Hii1, ione, ione, desca, 
                    Gii, ione, ione, descb, zero, temp, ione, ione, desca);


            /* Gii now has the matrix Hi-1,i-1  */
            nz = pmo.mxllda_cond[i - 1] * pmo.mxlocc_cond[i - 1]; 
            zcopy_driver(nz, &H_tri[pmo.diag_begin[i - 1]], ione, Gii, ione);


            /* calculate Hi-1,i-1 - Hi-1,i * Gii^0 * Hi,i-1  */

            /* calculate  Hi-1,i-1 - tempi-1,i * Hi,i-1 = Hi-1,i-1 */
            zgemm_driver ("N", "N", n1, n1, n2, mone, temp, ione, ione, desca, 
                    Hi1i, ione, ione, descd, one, Gii, ione, ione, descc);


            /* now Gii store the matrix Hi-1,i-1 - Hi-1,i * Gii^0 * Hi,i-1
             * Gi-1,i-1, stored in Gii, = Hii^(-1)
             */

            matrix_inverse_driver(Gii, descc);

            nz = pmo.mxllda_cond[i - 1] * pmo.mxlocc_cond[i - 1];
            zcopy_driver (nz, Gii, ione, &Grow[n_begin1[i - 1]], ione);
            zcopy_driver (nz, Gii, ione, &Gcol[n_begin1[i - 1]], ione);


            /*  temp[i,i-1] = - Hi,i-1 * Gi-1,i-1  */
            n1 = ni[i];
            n2 = ni[i - 1];

            desca = &pmo.desc_cond[ (i   + (i-1) * N ) * DLEN ];
            descb = &pmo.desc_cond[ (i-1 + (i-1) * N ) * DLEN ];
            zgemm_driver ("N", "N", n1, n2, n2, mone, Hi1i, ione, ione, desca, 
                    Gii, ione, ione, descb, zero, temp, ione, ione, desca);


            /* Gj, i-1 = G0_j,i * temp , j = i, n-1 , eq. 18*/
            for (j = i; j < N; j++)
            {

                n3 = ni[j];

                descc = &pmo.desc_cond[ (j +  i    * N ) * DLEN ];
                descd = &pmo.desc_cond[ (j + (i-1) * N ) * DLEN ];

                zgemm_driver ("N", "N", n3, n2, n1, one, &Grow[n_begin1[j]], ione, ione, descc, 
                        temp, ione, ione, desca, zero, Hii, ione, ione, descd);

                nz = pmo.mxllda_cond[j] * pmo.mxlocc_cond[i - 1];
                zcopy_driver (nz, Hii, ione, &Grow[n_begin1[j]], ione);
            }

            /* Gi-1, j, = -G_i-1,i-1 * Hi-1,i * G0_i,j , j = i, n-1  eq. 19*/
            desca = &pmo.desc_cond[ (i-1 +  i    * N ) * DLEN ];
            descb = &pmo.desc_cond[ (i-1 + (i-1) * N ) * DLEN ];
            zgemm_driver ("N", "N", n2, n1, n2, mone, Gii, ione, ione, descb, 
                    Hii1, ione, ione, desca, zero, temp, ione, ione, desca);


            for (j = i; j < N; j++)
            {

                n3 = ni[j];

                descc = &pmo.desc_cond[ (i   +  j * N ) * DLEN ];
                descd = &pmo.desc_cond[ (i-1 +  j * N ) * DLEN ];

                zgemm_driver ("N", "N", n3, n2, n1, one, temp, ione, ione, desca, 
                        &Gcol[n_begin1[j]], ione, ione, descc, zero, Hii, ione, ione, descd);

                nz = pmo.mxllda_cond[i-1] * pmo.mxlocc_cond[j];
                zcopy_driver (nz, Hii, ione, &Gcol[n_begin1[j]], ione);
            }




        }                           /* end  (i = N - 1; i > m + 1 ; i--) */

        /* ========================== Part-III starts here ========================== */
        /*        printf (" value of m ....  %d \n", m);        */


        /* call the identity matrix and allocate its memory */

        desca = &pmo.desc_cond[ (m + m * N) * DLEN]; 
        pmo_unitary_matrix(Imatrix_host, desca); 
        nz= ni[m] * ni[m];
        setvector_host_device (nz, sizeof(complex double), Imatrix_host, ione, ct.gpu_Imatrix, ione);

        /* get the interaction Hm,m+1 from input H_tri */


        Hii1 = &ct.gpu_Htri[pmo.offdiag_begin[m] ];
        Hi1i = &ct.gpu_Htri[pmo.lowoffdiag_begin[m] ];

        /* calculate: 1 - Gmm * Hm,m+1 * Gm+1,m+1 * Hm+1,m  */

        n1 = ni[m];
        n2 = ni[m + 1];

        desca = &pmo.desc_cond[ ( m   +  m    * N) * DLEN ];
        descb = &pmo.desc_cond[ ( m   + (m+1) * N) * DLEN ];
        descc = &pmo.desc_cond[ ( m+1 +  m    * N) * DLEN ];
        descd = &pmo.desc_cond[ ( m+1 + (m+1) * N) * DLEN ];


        /* eq. 25, the inverse part */
        /*  temp = Hm+1,m *  G0mm  */
        zgemm_driver ("N", "N", n2, n1, n1, one, Hi1i, ione, ione, descc, 
                &Grow[n_begin1[m]], ione, ione, desca, zero, temp, ione, ione, descc);
        /*  Gii = Hm,m+1 * G0 m+1,m+1 */
        zgemm_driver ("N", "N", n1, n2, n2, one, Hii1, ione, ione, descb,
                &Grow[n_begin1[m+1]], ione, ione, descd, zero, Hii, ione, ione, descb);

        /*  Hii = 1 - Hm,m+1 * Gm+1,m+1 * Hm+1,m * Gmm */
        nz = pmo.mxllda_cond[m] * pmo.mxlocc_cond[m];
        zcopy_driver(nz, Imatrix, ione, Gii, ione);
        zgemm_driver ("N", "N", n1, n1, n2, mone, Hii, ione, ione, descb,
                temp, ione, ione, descc, one, Gii, ione, ione, desca);

        /* get the inverse of Hii , see second bracket in eq. 25 */
        matrix_inverse_driver(Gii, desca);

        /* gnu_Hii = Hm+1,m * Gmm * (... )^(-1)   eq.25 */
        zgemm_driver ("N", "N", n2, n1, n1, one, temp, ione, ione, descc,
                Gii, ione, ione, desca, zero, Hii, ione, ione, descc);

        /* eq. 25, first part for j <= m */
        for (j = 0; j <= m; j++)
        {
            n3 = ni[j];

            desce = &pmo.desc_cond[ ( j + m * N) * DLEN ];

            zgemm_driver ("N", "N", n3, n1, n1, one, &Grow[n_begin1[j]], ione, ione, desce, 
                    Gii, ione, ione, desca, zero, temp, ione, ione, desce);

            nz = pmo.mxllda_cond[j] * pmo.mxlocc_cond[m];
            zcopy_driver (nz, temp, ione, &Grow[n_begin1[j]], ione);
        }

        for (j = m + 1; j < N; j++)
        {
            n3 = ni[j];
            desce = &pmo.desc_cond[ ( j + (m+1) * N) * DLEN ];
            descf = &pmo.desc_cond[ ( j +  m   * N) * DLEN ];

            zgemm_driver ("N", "N", n3, n1, n2, mone, &Grow[n_begin1[j]], ione, ione, desce,
                    Hii, ione, ione, descc, zero, temp, ione, ione, descf);

            nz = pmo.mxllda_cond[j] * pmo.mxlocc_cond[m];
            zcopy_driver (nz, temp, ione, &Grow[n_begin1[j]], ione);

        }

        /* eq. 24  */
        /* calculate Gmm * Hm,m+1 = temp[m,m+1] */

        zgemm_driver ("N", "N", n1, n2, n1, one, &Gcol[n_begin1[m]], ione, ione, desca,
                Hii1, ione, ione, descb, zero, temp, ione, ione, descb);



        /* calculate  temp[m,m+1] * Gm+1,m+1 = temp2[m,m+1] */
        zgemm_driver ("N", "N", n1, n2, n2, one, temp, ione, ione, descb,
                &Gcol[n_begin1[m+1]], ione, ione, descd, zero, Hii, ione, ione, descb);

        nz = pmo.mxllda_cond[m] * pmo.mxlocc_cond[m];
        zcopy_driver (nz, Imatrix, ione, Gii, ione);

        zgemm_driver ("N", "N", n1, n1, n2, mone, Hii, ione, ione, descb,
                Hi1i, ione, ione, descc, one, Gii, ione, ione, desca);
        /* Now 1 - Gmm * Hm,m+1 * Gm+1,m+1 * Hm+1,m is stored in Gii[m,m] */

        matrix_inverse_driver(Gii, desca);


        /* Hii store:  (1-G0mm * Hm,m+1 * G0m+1,m+1 * Hm+1,m)^-1 * G0mm * Hm,m+1 */
        zgemm_driver ("N", "N", n1, n2, n1, one, Gii, ione, ione, desca,
                temp, ione, ione, descb, zero, Hii, ione, ione, descb);

        for (j = 0; j <= m; j++)
        {
            n3 = ni[j];

            desce = &pmo.desc_cond[ ( m + j * N) * DLEN ];

            zgemm_driver ("N", "N", n1, n3, n1, one, Gii, ione, ione, desca, 
                    &Gcol[n_begin1[j]], ione, ione, desce, zero, temp, ione, ione, desce);

            nz = pmo.mxllda_cond[m] * pmo.mxlocc_cond[j];
            zcopy_driver (nz, temp, ione, &Gcol[n_begin1[j]], ione);
        }


        for (j = m + 1; j < N; j++)
        {
            n3 = ni[j];

            desce = &pmo.desc_cond[ ( m   + j * N) * DLEN ];
            descf = &pmo.desc_cond[ ( m+1 + j * N) * DLEN ];
            zgemm_driver ("N", "N", n1, n3, n2, mone, Hii, ione, ione, descb, 
                    &Gcol[n_begin1[j]], ione, ione, descf, zero, temp, ione, ione, desce);

            nz = pmo.mxllda_cond[m] * pmo.mxlocc_cond[j];
            zcopy_driver (nz, temp, ione, &Gcol[n_begin1[j]], ione);

        }



    }                           /* if statement ends here */
    /* ================================================================== */


    n4 = totrow * maxcol;
    getvector_device_host(n4, sizeof(complex double), ct.gpu_Grow, ione, G_row_host, ione);
    n4 = totcol * maxrow;
    getvector_device_host(n4, sizeof(complex double), ct.gpu_Gcol, ione, G_col_host, ione);



    my_free(n_begin1);
    my_free( Hii_host );
    my_free( Gii_host );
    my_free( temp_host );
    my_free( Imatrix_host );
}

