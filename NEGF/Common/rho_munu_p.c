/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include "md.h"
#include "pmo.h"



void rho_munu_p (complex double * rho_mn, complex double * green_C, complex double * sigma_L, int iprobe)
{

    complex double *temp;
    complex double one, zero;
    int i, j, n_green, n1, n2;
    char fcd_n = 'N', fcd_c = 'C';

    int N, *ni, nrow, ncol, nL, N1;
    int *desca, *descb, *descc, *desc_col, *desc_lead;
    int *descd, maxrow, maxcol, ione =1;

    N = ct.num_blocks;
    ni = ct.block_dim;


	N1 = cei.probe_in_block[iprobe - 1];
    nL = ct.block_dim[N1];
    nrow = pmo.mxllda_cond[N1];
    ncol = pmo.mxlocc_cond[N1];

    desc_lead = &pmo.desc_cond[ (N1 + N1 * ct.num_blocks) * DLEN ];
    desc_col= &pmo.desc_cond[ (N1 * ct.num_blocks) * DLEN ];


	one = 1.0;
    zero = 0.0;

    maxrow = 0;
    maxcol = 0;
    for (i = 0; i < N; i++)
    {
        maxrow = max(maxrow, pmo.mxllda_cond[i]);
        maxcol = max(maxcol, pmo.mxlocc_cond[i]);
    }

    /*  Gamma = i*(sigma - simga^+)  */
    for (i = 0; i < nrow * ncol; i++)
    {
        sigma_L[i] = -cimag(sigma_L[i]);
    }

    n_green = 0;

    n1 = maxrow * maxcol;
    my_malloc_init( temp, n1, complex double );

    for (i = 0; i < N - 1; i++)
    {
        n1 = ni[i];
        n2 = ni[i + 1];

        /*  temp = G_i0 * Gamma  */
        desca = &desc_col[i * DLEN];
        descb = &pmo.desc_cond[ (i + i * ct.num_blocks ) * DLEN ];
        descc = &desc_col[ (i+1) * DLEN];
        descd = &pmo.desc_cond[ (i + (i+1) * ct.num_blocks ) * DLEN ];
        PZGEMM (&fcd_n, &fcd_n, &n1, &nL, &nL, &one, &green_C[n_green], &ione, &ione, desca,
                sigma_L, &ione, &ione, desc_lead, &zero, temp, &ione, &ione, desca);

        /* rho_mn (i,i) = temp * G_i0^, the block (i,i) */

        PZGEMM (&fcd_n, &fcd_c, &n1, &n1, &nL, &one, temp, &ione, &ione, desca, 
                &green_C[n_green], &ione, &ione, desca, &zero,
                &rho_mn[pmo.diag_begin[i] ], &ione, &ione, descb);

        /* rho_mn (i,i+1) = temp * G_i+10^, the block (i,i) */
        n_green += pmo.mxllda_cond[i] * maxcol;
        PZGEMM (&fcd_n, &fcd_c, &n1, &n2, &nL, &one, temp, &ione, &ione, desca, 
               &green_C[n_green], &ione, &ione, descc, &zero,
               &rho_mn[pmo.offdiag_begin[i]], &ione, &ione, descd);

    }

    /* calculate the last block  */
    n1 = ni[N - 1];

    desca = &desc_col[ (N-1) * DLEN ];
    descb = &pmo.desc_cond[ (N-1 + (N-1) * ct.num_blocks ) * DLEN ];
    PZGEMM (&fcd_n, &fcd_n, &n1, &nL, &nL, &one, &green_C[n_green], &ione, &ione, desca,
           sigma_L, &ione, &ione, desc_lead, &zero, temp, &ione, &ione, desca);

    PZGEMM (&fcd_n, &fcd_c, &n1, &n1, &nL, &one, temp, &ione, &ione, desca, 
           &green_C[n_green], &ione, &ione, desca, &zero,
           &rho_mn[pmo.diag_begin[N-1] ], &ione, &ione, descb);




    my_free( temp );


}
