/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
/*
 *  Calculate the Sigma matrix,  eq. 20 of PRB 65, 165401   
 *  V = HLC - e * SLC
 *  HLC, SLC:  nL * nC real matrix
 *  Green:  nL * nL complex matrix
 *  sigma: nC * nC complex matrix as output
 *  sigma = V^+ Green * V
 */

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include "main.h"
#include "pmo.h"



void Sigma_cond_p (complex double *sigma, REAL *HLC, REAL *SLC, REAL eneR, REAL eneI, complex double *green, int iprobe)
{

    complex double *ch, *ch1;
    int i;
    char fcd_n = 'N', fcd_t = 'T';
    complex double alpha, beta;
    int nmax, ndim, ione =1, nrow, ncol, n0, n1; 
	int *desca, *descb, *descc, *descd;


    alpha.r = 1.0;
    alpha.i = 0.0;
    beta.r = 0.0;
    beta.i = 0.0;

    n0 = cei.probe_in_block[iprobe - 1];  /* n0 => block index */

    ndim = ct.block_dim[n0];
    nmax = lcr[iprobe].num_states;

	
    nrow = pmo.mxllda_cond[n0];
    ncol = pmo.mxlocc_lead[iprobe-1];
    n1 = nrow * ncol;

    desca = &pmo.desc_cond_lead[ (n0 + (iprobe - 1) * ct.num_blocks) * DLEN ]; /* (C,L) */
    descb = &pmo.desc_lead[ ( iprobe-1)       * DLEN ];                        /* (L,L) */
    descc = &pmo.desc_cond[ ( n0 + n0 * ct.num_blocks) * DLEN ];               /* (C,C) */
    

    my_malloc_init( ch, n1, complex double);
    my_malloc_init( ch1, n1, complex double);

    for (i = 0; i < n1; i++)
    {
        ch[i].r = eneR * SLC[i] - HLC[i] * Ha_eV;
        ch[i].i = eneI * SLC[i];
    }

	/*    ch1(C,L) = ch(C,L) * green(L,L)  */
    /*    sigma(C,C) = ch1(C,L) * ch(L,C)  */ 
		  

    PZGEMM (&fcd_n, &fcd_n, &ndim, &nmax, &nmax, &alpha, ch, &ione, &ione, desca,
            green, &ione, &ione, descb, &beta, ch1, &ione, &ione, desca);

    PZGEMM (&fcd_n, &fcd_t, &ndim, &ndim, &nmax, &alpha, ch1, &ione, &ione, desca,
            ch, &ione, &ione, desca, &beta, sigma, &ione, &ione, descc);



    my_free(ch);
    my_free(ch1);

}
