/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "main.h"
#include "pmo.h"


void matrix_inverse_blocknm (complex double * H_tri, int N_blocks, int * ni,
        int m, int n, complex double * Green_C)
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
 *   Green_C: block (m2,m1) will be returned which will be used to
 *   calculate the transmission.
 *   N:   number of blocks
 *   ni:  dimension of each block
 *   for example Hii is a ni[i] x ni[i] matrix
 *               Hi,i+1 is a ni[i] x ni[i+1] matrix 
 */

    int nmax, i, j, n1, n2, n3, n4;
    int *ipiv, *n_begin1;
    complex double *Hii, *temp, *Hii1;
    complex double *Hmm1, *temp2, *temp3, *identity;
    complex double *Gmm, *Gnm, *Gnn;
    complex double mone, one, zero;
    int ione = 1, ntot, k, maxrow, maxcol; 
	int *desca, *descb, *descc, *descd;

    mone = -1.0;
    one = 1.0;
    zero = 0.0;

/*  find the maximum dimension of the blocks  */
    assert(m <= n);
    assert(m != N_blocks);

    ntot = pmo.ntot;
    maxrow = 0;
    maxcol = 0;
    for (i = 0; i < N_blocks; i++)
    {
        maxrow = max(maxrow, pmo.mxllda_cond[i]);
        maxcol = max(maxcol, pmo.mxlocc_cond[i]);
    }

    my_malloc( n_begin1, N_blocks, int );
    my_malloc( ipiv, maxrow + pmo.mblock, int ); 

    my_malloc_init( Hii, maxrow * maxcol, complex double );
    my_malloc_init( Gnn, maxrow * maxcol, complex double );
    my_malloc_init( Gmm, maxrow * maxcol, complex double );
    my_malloc_init( Gnm, maxrow * maxcol, complex double );
    my_malloc_init( temp, maxrow * maxcol, complex double );


/*  n_begin: starting address of each diagonal block in H_tri and G_tri
 *  the Hi,i+1 block will start at n_begin[i] + ni[i] * ni[i]
 */

    n_begin1[0] = 0;
    for (i = 1; i < N_blocks; i++)
    {
        n_begin1[i] = n_begin1[i - 1] + pmo.mxllda_cond[i - 1] * maxcol;
    }

    /* printf (" while matrix diag: iprobe, block =  %d %d \n", iprobe, m); */

/* ========================== Part-I ==================================== */

   
/*  calculate the inverse of the first block  */

    for (i = 0; i < pmo.mxllda_cond[0] * pmo.mxlocc_cond[0]; i++)
    {
        Hii[i] = H_tri[ i];
    }

    desca = &pmo.desc_cond[0];
    get_inverse_block_p (Hii, Gmm, ipiv, desca);


/*----------------------------
    for(idx =0; idx < n1; idx++)
    {
        if(pct.gridpe ==0) 
        printf (" Green_C %d %d %f %f \n", i, j, &Green_C[[0] + idx].r, &Green_C[[0] + idx].i);  
    }
----------------------------*/
	
/*  iterate to get one more block  */
    for (i = 0; i < m; i++)
    {
        /* get the interaction  Hi,i+1  from input H_tri 
         * Hii1 is a pointer only
         */
        Hii1 = &H_tri[pmo.offdiag_begin[i] ];

        /* Hii now has the matrix Hi+1,i+1  */
        for (j = 0; j < pmo.mxllda_cond[i + 1] * pmo.mxlocc_cond[i + 1]; j++)
        {
            Hii[j] = H_tri[j + pmo.diag_begin[i + 1]];
        }


        /* calculate Hi+1,i+1 - Hi+1,i * Gii^0 * Hi,i+1  */
		
        n1 = ni[i + 1];
        n2 = ni[i];

        desca = &pmo.desc_cond[ ( i   + (i+1) * N_blocks) * DLEN ];
        descb = &pmo.desc_cond[ ( i   +  i    * N_blocks) * DLEN ];
        descc = &pmo.desc_cond[ ( i+1 +  i    * N_blocks) * DLEN ];
        descd = &pmo.desc_cond[ ( i+1 + (i+1) * N_blocks) * DLEN ];

        /* calculate Hi+1,i * Gmm^0 = tempi+1,i */
        PZGEMM ("C", "N", &n1, &n2, &n2, &one, Hii1, &ione, &ione, desca, 
               Gmm, &ione, &ione, descb, &zero, temp, &ione, &ione, descc);
        /* calculate Hi+1,i+1 - tempi+1,i * Hi,i+1 = Hi+1,i+1 */
        PZGEMM ("N", "N", &n1, &n1, &n2, &mone, temp, &ione, &ione, descc, 
               Hii1, &ione, &ione, desca, &one, Hii, &ione, &ione, descd);


        /* now Hii store the matrix Hi+1,i+1 - Hi+1,i * Gii^0 * Hi,i+1
         * Gi+1,i+1, stored in Gii, = Hii^(-1)
         */

        get_inverse_block_p (Hii, Gmm, ipiv, descd);


    }                           /* end  for(i = 0; i < m; i++) */


    
/* ========================== Part-II starts here ============================ */


/*  calculate the inverse of the last block  */

    for (i = 0; i < pmo.mxllda_cond[N_blocks - 1] * pmo.mxlocc_cond[N_blocks - 1]; i++)
    {
        Hii[i] = H_tri[pmo.diag_begin[N_blocks - 1] + i];
    }

    desca = &pmo.desc_cond[ (N_blocks-1 + (N_blocks-1) * N_blocks) * DLEN];
    get_inverse_block_p (Hii, Gnn, ipiv, desca);


/*  iterate to get one more block  */


    for (i = N_blocks - 1; i > n ; i--) 
    {
        /* get the interaction  Hi-1,i  from input H_tri 
         * Hii1 is a pointer only
         */
        Hii1 = &H_tri[pmo.offdiag_begin[i-1] ];

        /* Hii now has the matrix Hi+1,i+1  */

        for (j = 0; j < pmo.mxllda_cond[i - 1] * pmo.mxlocc_cond[i - 1]; j++)
        {
            Hii[j] = H_tri[j + pmo.diag_begin[i - 1]];
        }


        /* calculate Hi-1,i-1 - Hi-1,i * Gii^0 * Hi,i-1  */

        n1 = ni[i - 1];
        n2 = ni[i];
        desca = &pmo.desc_cond[ (i-1 +  i    * N_blocks ) * DLEN ];
        descb = &pmo.desc_cond[ (i   +  i    * N_blocks ) * DLEN ];
        descc = &pmo.desc_cond[ (i-1 + (i-1) * N_blocks ) * DLEN ];

        /* calculate Hi-1,i * Gii^0 = tempi-1,i */
        PZGEMM ("N", "N", &n1, &n2, &n2, &one, Hii1, &ione, &ione, desca, 
                Gnn, &ione, &ione, descb, &zero, temp, &ione, &ione, desca);

        /* calculate Hi-1,i-1 - tempi-1,i * Hi,i-1 = Hi-1,i-1 */
        PZGEMM ("N", "C", &n1, &n1, &n2, &mone, temp, &ione, &ione, desca, 
                Hii1, &ione, &ione, desca, &one, Hii, &ione, &ione, descc);


        /* now Hii store the matrix Hi-1,i-1 - Hi-1,i * Gii^0 * Hi,i-1
         * Gi-1,i-1, stored in Gii, = Hii^(-1)
         */

        get_inverse_block_p (Hii, Gnn, ipiv, descc);
    }

/* ========================== Part-III starts here ============================ */

    /*  copy Gnn to Gnm  */

    n1 = pmo.mxllda_cond[n] * pmo.mxlocc_cond[n];
    zcopy (&n1, Gnn, &ione, Gnm, &ione);

    for (i = n; i > m+1 ; i--) 
    {
        /* get the interaction  Hi-1,i  from input H_tri 
         * Hii1 is a pointer only
         */
        Hii1 = &H_tri[pmo.offdiag_begin[i-1] ];

        /* Hii now has the matrix Hi+1,i+1  */

        for (j = 0; j < pmo.mxllda_cond[i - 1] * pmo.mxlocc_cond[i - 1]; j++)
        {
            Hii[j] = H_tri[j + pmo.diag_begin[i - 1]];
        }


        /* calculate Hi-1,i-1 - Hi-1,i * Gii^0 * Hi,i-1  */

        n1 = ni[i - 1];
        n2 = ni[i];
        desca = &pmo.desc_cond[ (i-1 +  i    * N_blocks ) * DLEN ];
        descb = &pmo.desc_cond[ (i   +  i    * N_blocks ) * DLEN ];
        descc = &pmo.desc_cond[ (i-1 + (i-1) * N_blocks ) * DLEN ];

        /* calculate Hi-1,i * Gii^0 = tempi-1,i */
        PZGEMM ("N", "N", &n1, &n2, &n2, &one, Hii1, &ione, &ione, desca, 
                Gnn, &ione, &ione, descb, &zero, temp, &ione, &ione, desca);

        /* calculate Hi-1,i-1 - tempi-1,i * Hi,i-1 = Hi-1,i-1 */
        PZGEMM ("N", "C", &n1, &n1, &n2, &mone, temp, &ione, &ione, desca, 
                Hii1, &ione, &ione, desca, &one, Hii, &ione, &ione, descc);


        /* now Hii store the matrix Hi-1,i-1 - Hi-1,i * Gii^0 * Hi,i-1
         * Gi-1,i-1, stored in Gnn, = Hii^(-1)
         */

        get_inverse_block_p (Hii, Gnn, ipiv, descc);


        /*  temp[n,i-1] = - G[n,i) Hi,i-1  */
        n1 = ni[n];
        n2 = ni[i];
        n3 = ni[i - 1];

        desca = &pmo.desc_cond[ (n   +  i    * N_blocks ) * DLEN ];
        descb = &pmo.desc_cond[ (i-1 +  i    * N_blocks ) * DLEN ];
        descc = &pmo.desc_cond[ (n   + (i-1) * N_blocks ) * DLEN ];
        PZGEMM ("N", "C", &n1, &n3, &n2, &mone, Gnm, &ione, &ione, desca, 
                Hii1, &ione, &ione, descb, &zero, temp, &ione, &ione, descc);


        /* Gn, i-1 =  temp[n,i-1] * G[i-1,i-1] , j = i, n-1 */


        n1 = ni[n];
        n2 = ni[i - 1];

        desca = &pmo.desc_cond[ (n   + (i-1) * N_blocks ) * DLEN ];
        descc = &pmo.desc_cond[ (i-1 + (i-1) * N_blocks ) * DLEN ];

        PZGEMM ("N", "N", &n1, &n2, &n2, &one, temp, &ione, &ione, desca, 
                Gnn, &ione, &ione, descc, &zero, Gnm, &ione, &ione, desca);


    }

    /* now, Gnn store the Green function G^0_{m+1, m+1) 
     *      Gnm store the Green function G^0_{n, m+1}
     *      Gmm store the Green function G^0_{m, m}
     */


/* ========================== Part-IV starts here ========================== */
/*        printf (" value of m ....  %d \n", m);        */


    /* call the identity matrix and allocate its memory */

    my_malloc_init( identity, maxrow * maxcol, complex double );  
    desca = &pmo.desc_cond[ (m + m * N_blocks) * DLEN]; 
    pmo_unitary_matrix(identity, desca); 


    /* allocate memory for temp2 and temp3 */
    my_malloc_init( temp2, maxrow * maxcol, complex double );
    my_malloc_init( temp3, maxrow * maxcol, complex double );

    /* get the interaction Hm,m+1 from input H_tri */
    Hmm1 = &H_tri[pmo.offdiag_begin[m] ];

    /* calculate: 1 - Gmm * Hm,m+1 * Gm+1,m+1 * Hm+1,m  */

    n1 = ni[m];
    n2 = ni[m + 1];

    desca = &pmo.desc_cond[ ( m   +  m    * N_blocks) * DLEN ];
    descb = &pmo.desc_cond[ ( m   + (m+1) * N_blocks) * DLEN ];
    descc = &pmo.desc_cond[ ( m+1 + (m+1) * N_blocks) * DLEN ];


    /* calculate Gmm * Hm,m+1 = temp[m,m+1] */
    PZGEMM ("N", "N", &n1, &n2, &n1, &one, Gmm, &ione, &ione, desca, 
            Hmm1, &ione, &ione, descb, &zero, temp, &ione, &ione, descb);


    /* calculate  temp[m,m+1] * Gm+1,m+1 = temp2[m,m+1] */
    PZGEMM ("N", "N", &n1, &n2, &n2, &one, temp, &ione, &ione, descb, 
            Gnn, &ione, &ione, descc, &zero, temp2, &ione, &ione, descb);


    /* calculate identity[m,m] - temp2[m,m+1] * Hm+1,m = identity[m,m] */
    PZGEMM ("N", "C", &n1, &n1, &n2, &mone, temp2, &ione, &ione, descb, 
            Hmm1, &ione, &ione, descb, &one, identity, &ione, &ione, desca);

    /* Now 1 - Gmm * Hm,m+1 * Gm+1,m+1 * Hm+1,m is stored in identity[m,m] */

    /* get the inverse of identity[m,m] and stored in temp2[m,m]           */
    get_inverse_block_p (identity, temp2, ipiv, desca); /* check ipiv for iprobe */


    /* calculate G[m,m] = temp2[m,m] * G^0[m,m] */
    PZGEMM ("N", "N", &n1, &n1, &n1, &one, temp2, &ione, &ione, desca, 
            Gmm, &ione, &ione, desca, &zero, temp3, &ione, &ione, desca);

    n4 = pmo.mxllda_cond[m] * pmo.mxlocc_cond[m];
    zcopy (&n4, temp3, &ione, Gmm, &ione);

    if(n == m)
    {

        n4 = pmo.mxllda_cond[m] * pmo.mxlocc_cond[m];
        zcopy (&n4, Gmm, &ione, Green_C, &ione);
    }
    else
    {


        n1 = ni[n];
        n2 = ni[m];
        n3 = ni[m+1];

        desca = &pmo.desc_cond[ ( n   + (m+1) * N_blocks) * DLEN ];
        descb = &pmo.desc_cond[ ( n   +  m    * N_blocks) * DLEN ];
        descc = &pmo.desc_cond[ ( m   + (m+1) * N_blocks) * DLEN ];
        descd = &pmo.desc_cond[ ( m   +  m    * N_blocks) * DLEN ];


        /* temp[n,m] = -G^0[n, m+1] * Hm+1,m */
        PZGEMM ("N", "C", &n1, &n2, &n3, &mone,  Gnm, &ione, &ione, desca, 
                Hmm1, &ione, &ione, descc, &zero, temp, &ione, &ione, descb);

        /* Gnm = temp[n,m] * G[m,m] */
        PZGEMM ("N", "N", &n1, &n2, &n2, &one,  temp, &ione, &ione, descb, 
                Gmm, &ione, &ione, descd, &zero, Green_C, &ione, &ione, descb);
    }

    my_free(n_begin1);
    my_free(ipiv);
    my_free( Hii );
    my_free( Gmm );
    my_free( Gnn );
    my_free( Gnm );
    my_free( temp );

    my_free( temp2 );
    my_free( temp3 );
    my_free( identity );


}   
/* ================================================================== */

