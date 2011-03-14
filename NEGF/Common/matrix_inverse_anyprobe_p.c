/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "md.h"
#include "pmo.h"


void matrix_inverse_anyprobe_p (doublecomplex * H_tri, int N, int * ni, int iprobe, doublecomplex * Green_C)
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
    int *ipiv, *n_begin1;
    doublecomplex *Hii, *Gii, *temp, *Hii1;
    doublecomplex *Hmm1, *temp2, *temp3, *identity;
    doublecomplex mone, one, zero;
    int ione = 1, ntot, k, maxrow, maxcol; 
	int *desca, *descb, *descc, *descd;

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

    for (i = 0; i < pmo.mxllda_cond[0] * pmo.mxlocc_cond[0]; i++)
    {
        Hii[i].r = H_tri[ i].r;
        Hii[i].i = H_tri[ i].i;
    }

    desca = &pmo.desc_cond[0];
    get_inverse_block_p (Hii, Gii, ipiv, desca);

    n1 = pmo.mxllda_cond[0] * pmo.mxlocc_cond[0];
    zcopy (&n1, Gii, &ione, &Green_C[n_begin1[0]], &ione);

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
            Hii[j].r = H_tri[j + pmo.diag_begin[i + 1]].r;
            Hii[j].i = H_tri[j + pmo.diag_begin[i + 1]].i;
        }


        /* calculate Hi+1,i+1 - Hi+1,i * Gii^0 * Hi,i+1  */
		
        n1 = ni[i + 1];
        n2 = ni[i];

        desca = &pmo.desc_cond[ ( i   + (i+1) * N) * DLEN ];
        descb = &pmo.desc_cond[ ( i   +  i    * N) * DLEN ];
        descc = &pmo.desc_cond[ ( i+1 +  i    * N) * DLEN ];
        descd = &pmo.desc_cond[ ( i+1 + (i+1) * N) * DLEN ];

        /* calculate Hi+1,i * Gii^0 = tempi+1,i */
        PZGEMM ("T", "N", &n1, &n2, &n2, &one, Hii1, &ione, &ione, desca, 
               Gii, &ione, &ione, descb, &zero, temp, &ione, &ione, descc);
        /* calculate Hi+1,i+1 - tempi+1,i * Hi,i+1 = Hi+1,i+1 */
        PZGEMM ("N", "N", &n1, &n1, &n2, &mone, temp, &ione, &ione, descc, 
               Hii1, &ione, &ione, desca, &one, Hii, &ione, &ione, descd);


        /* now Hii store the matrix Hi+1,i+1 - Hi+1,i * Gii^0 * Hi,i+1
         * Gi+1,i+1, stored in Gii, = Hii^(-1)
         */

        get_inverse_block_p (Hii, Gii, ipiv, descd);

        n1 = pmo.mxllda_cond[i + 1] * pmo.mxlocc_cond[i + 1];
        zcopy (&n1, Gii, &ione, &Green_C[n_begin1[i + 1]], &ione);


        /*  tempi,i+1 = - Hi,i+1 * Gi+1,i+1  */
        n1 = ni[i];
        n2 = ni[i + 1];
        desca = &pmo.desc_cond[ ( i   + (i+1) * N) * DLEN ];
        descb = &pmo.desc_cond[ ( i+1 + (i+1) * N) * DLEN ];
        PZGEMM ("N", "N", &n1, &n2, &n2, &mone, Hii1, &ione, &ione, desca, 
                Gii, &ione, &ione, descb, &zero, temp, &ione, &ione, desca);


        /* Gj, i+1 = G0_j,i * tempi,i+1 ; j = 0, i */
        for (j = 0; j <= i; j++)
        {

            n1 = ni[j];
            n2 = ni[i];
            n3 = ni[i + 1];

            desca = &pmo.desc_cond[ ( j   +  i    * N) * DLEN ];
            descb = &pmo.desc_cond[ ( i   + (i+1) * N) * DLEN ];
            descc = &pmo.desc_cond[ ( j   + (i+1) * N) * DLEN ];

            PZGEMM ("N", "N", &n1, &n3, &n2, &one, &Green_C[n_begin1[j]], &ione, &ione, desca, 
                    temp, &ione, &ione, descb, &zero, Hii, &ione, &ione, descc);

            n1 = pmo.mxllda_cond[j] * pmo.mxlocc_cond[i + 1];
            zcopy (&n1, Hii, &ione, &Green_C[n_begin1[j]], &ione);
        }
    }                           /* end  for(i = 0; i < m; i++) */


    
/* ========================== Part-II starts here ============================ */

	if(m < N - 1 )  
	{

/*  calculate the inverse of the last block  */

    for (i = 0; i < pmo.mxllda_cond[N - 1] * pmo.mxlocc_cond[N - 1]; i++)
    {
        Hii[i].r = H_tri[pmo.diag_begin[N - 1] + i].r;
        Hii[i].i = H_tri[pmo.diag_begin[N - 1] + i].i;
    }

    desca = &pmo.desc_cond[ (N-1 + (N-1) * N) * DLEN];
    get_inverse_block_p (Hii, Gii, ipiv, desca);


    n1 = pmo.mxllda_cond[N - 1] * pmo.mxlocc_cond[N - 1];
    zcopy (&n1, Gii, &ione, &Green_C[n_begin1[N - 1]], &ione);

/*  iterate to get one more block  */

    for (i = N - 1; i > m + 1 ; i--) 
/*    for (i = N - 1; i > m; i--)   */ 
    {
        /* get the interaction  Hi-1,i  from input H_tri 
         * Hii1 is a pointer only
         */
        Hii1 = &H_tri[pmo.offdiag_begin[i-1] ];

        /* Hii now has the matrix Hi+1,i+1  */

        for (j = 0; j < pmo.mxllda_cond[i - 1] * pmo.mxlocc_cond[i - 1]; j++)
        {
            Hii[j].r = H_tri[j + pmo.diag_begin[i - 1]].r;
            Hii[j].i = H_tri[j + pmo.diag_begin[i - 1]].i;
        }


        /* calculate Hi-1,i-1 - Hi-1,i * Gii^0 * Hi,i-1  */

        n1 = ni[i - 1];
        n2 = ni[i];
        desca = &pmo.desc_cond[ (i-1 +  i    * N ) * DLEN ];
        descb = &pmo.desc_cond[ (i   +  i    * N ) * DLEN ];
        descc = &pmo.desc_cond[ (i-1 + (i-1) * N ) * DLEN ];

        /* calculate Hi-1,i * Gii^0 = tempi-1,i */
        PZGEMM ("N", "N", &n1, &n2, &n2, &one, Hii1, &ione, &ione, desca, 
                Gii, &ione, &ione, descb, &zero, temp, &ione, &ione, desca);

        /* calculate Hi-1,i-1 - tempi-1,i * Hi,i-1 = Hi-1,i-1 */
        PZGEMM ("N", "T", &n1, &n1, &n2, &mone, temp, &ione, &ione, desca, 
                Hii1, &ione, &ione, desca, &one, Hii, &ione, &ione, descc);


        /* now Hii store the matrix Hi-1,i-1 - Hi-1,i * Gii^0 * Hi,i-1
         * Gi-1,i-1, stored in Gii, = Hii^(-1)
         */

        get_inverse_block_p (Hii, Gii, ipiv, descc);


        n1 = pmo.mxllda_cond[i - 1] * pmo.mxlocc_cond[i - 1];
        zcopy (&n1, Gii, &ione, &Green_C[n_begin1[i - 1]], &ione);


        /*  temp[i,i-1] = - Hi,i-1 * Gi-1,i-1  */
        n1 = ni[i];
        n2 = ni[i - 1];

        desca = &pmo.desc_cond[ (i-1 +  i    * N ) * DLEN ];
        descb = &pmo.desc_cond[ (i-1 + (i-1) * N ) * DLEN ];
        descc = &pmo.desc_cond[ (i   + (i-1) * N ) * DLEN ];
        PZGEMM ("T", "N", &n1, &n2, &n2, &mone, Hii1, &ione, &ione, desca, 
                Gii, &ione, &ione, descb, &zero, temp, &ione, &ione, descc);


        /* Gj, i-1 = G0_j,i * temp , j = i, n-1 */
        for (j = i; j < N; j++)
        {

            n1 = ni[j];
            n2 = ni[i];
            n3 = ni[i - 1];

            desca = &pmo.desc_cond[ (j +  i    * N ) * DLEN ];
            descb = &pmo.desc_cond[ (i + (i-1) * N ) * DLEN ];
            descc = &pmo.desc_cond[ (j + (i-1) * N ) * DLEN ];

            PZGEMM ("N", "N", &n1, &n3, &n2, &one, &Green_C[n_begin1[j]], &ione, &ione, desca, 
                    temp, &ione, &ione, descb, &zero, Hii, &ione, &ione, descc);


            n1 = pmo.mxllda_cond[j] * pmo.mxlocc_cond[i - 1];
            zcopy (&n1, Hii, &ione, &Green_C[n_begin1[j]], &ione);
        }

    }                           /* end  (i = N - 1; i > m + 1 ; i--) */

/* ========================== Part-III starts here ========================== */
/*        printf (" value of m ....  %d \n", m);        */


		/* call the identity matrix and allocate its memory */

		my_malloc_init( identity, maxrow * maxcol, doublecomplex );  
		desca = &pmo.desc_cond[ (m + m * N) * DLEN]; 
		pmo_unitary_matrix(identity, desca); 


		/* allocate memory for temp2 and temp3 */
		my_malloc_init( temp2, maxrow * maxcol, doublecomplex );
		my_malloc_init( temp3, maxrow * maxcol, doublecomplex );

		/* get the interaction Hm,m+1 from input H_tri */
		Hmm1 = &H_tri[pmo.offdiag_begin[m] ];

		/* calculate: 1 - Gmm * Hm,m+1 * Gm+1,m+1 * Hm+1,m  */

		n1 = ni[m];
		n2 = ni[m + 1];

		desca = &pmo.desc_cond[ ( m   +  m    * N) * DLEN ];
		descb = &pmo.desc_cond[ ( m   + (m+1) * N) * DLEN ];
		descc = &pmo.desc_cond[ ( m+1 + (m+1) * N) * DLEN ];

		/*
		   if(pct.gridpe ==0)
		   {
		   printf (" desca 0....  %d %d %d \n", desca[0], descb[0], descc[0]);
		   printf (" desca 1....  %d %d %d \n", desca[1], descb[1], descc[1]);
		   printf (" desca 2....  %d %d %d \n", desca[2], descb[2], descc[2]);
		   printf (" desca 3....  %d %d %d \n", desca[3], descb[3], descc[3]);
		   printf (" desca 4....  %d %d %d \n", desca[4], descb[4], descc[4]);
		   printf (" desca 5....  %d %d %d \n", desca[5], descb[5], descc[5]);
		   printf (" desca 6....  %d %d %d \n", desca[6], descb[6], descc[6]);
		   printf (" desca 7....  %d %d %d \n", desca[7], descb[7], descc[7]);
		   printf (" desca 8....  %d %d %d \n", desca[8], descb[8], descc[8]);
		   }
		 */

		/* calculate Gmm * Hm,m+1 = temp[m,m+1] */
		PZGEMM ("N", "N", &n1, &n2, &n1, &one, &Green_C[n_begin1[m]], &ione, &ione, desca, 
				Hmm1, &ione, &ione, descb, &zero, temp, &ione, &ione, descb);


		/* calculate  temp[m,m+1] * Gm+1,m+1 = temp2[m,m+1] */
		PZGEMM ("N", "N", &n1, &n2, &n2, &one, temp, &ione, &ione, descb, 
				&Green_C[n_begin1[m + 1]], &ione, &ione, descc, &zero, temp2, &ione, &ione, descb);


		/* calculate identity[m,m] - temp2[m,m+1] * Hm+1,m = identity[m,m] */
		PZGEMM ("N", "T", &n1, &n1, &n2, &mone, temp2, &ione, &ione, descb, 
				Hmm1, &ione, &ione, descb, &one, identity, &ione, &ione, desca);

		/* Now 1 - Gmm * Hm,m+1 * Gm+1,m+1 * Hm+1,m is stored in identity[m,m] */

		/* get the inverse of identity[m,m] and stored in temp2[m,m]           */
		get_inverse_block_p (identity, temp2, ipiv, desca); /* check ipiv for iprobe */


		/* calculate temp2[m,m] * temp[m,m+1] = temp3[m,m+1] */
		PZGEMM ("N", "N", &n1, &n2, &n1, &one, temp2, &ione, &ione, desca, 
				temp, &ione, &ione, descb, &zero, temp3, &ione, &ione, descb);


		/* Gmj = temp2[m,m] * Gmj - temp3[m,m+1] * Gm+1,j; j = 0, N-1 */
		/* this will be calculated in two parts                      */


		/* Gmj = temp2[m,m] * Gmj; j = 0, m */
		/* which is: Gjm = (Gmj)^T = Gjm * temp2[m,m]; j = 0, m */
		for (j = 0; j <= m; j++)
		{
			n3 = ni[j];

			descc = &pmo.desc_cond[ (j + m * N ) * DLEN ];

			PZGEMM ("N", "T", &n3, &n1, &n1, &one,  &Green_C[n_begin1[j]], &ione, &ione, descc, 
					temp2, &ione, &ione, desca, &zero, temp, &ione, &ione, descc);

			n4 = pmo.mxllda_cond[j] * pmo.mxlocc_cond[m];
			zcopy (&n4, temp, &ione, &Green_C[n_begin1[j]], &ione);
		}


		/* Gmj = - temp3[m,m+1] * Gm+1,j; j = m+1, N-1 */
		/* which is: Gjm = (Gmj)^T = - Gj,m+1 * temp3[m+1,m]; j = m+1, N-1 */
		for (j = m + 1; j < N; j++)
		{
			n3 = ni[j];

			descc = &pmo.desc_cond[ (j +  m    * N ) * DLEN ];
			descd = &pmo.desc_cond[ (j + (m+1) * N ) * DLEN ];

			PZGEMM ("N", "T", &n3, &n1, &n2, &mone, &Green_C[n_begin1[j]], &ione, &ione, descd, 
					temp3, &ione, &ione, descb, &zero, temp, &ione, &ione, descc);

			n4 = pmo.mxllda_cond[j] * pmo.mxlocc_cond[m];
			zcopy (&n4, temp, &ione, &Green_C[n_begin1[j]], &ione);

		}


		my_free( temp2 );
		my_free( temp3 );
		my_free( identity );

	}                           /* if statement ends here */
	/* ================================================================== */

	my_free(n_begin1);
	my_free(ipiv);
	my_free( Hii );
	my_free( Gii );
	my_free( temp );
}
