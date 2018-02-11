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
#include "init_var.h"
#include "LCR.h"
#include "pmo.h"



void Sigma_p (complex double *sigma, complex double *ch, complex double
        *ch01, complex double *ch10, complex double *green, int iprobe)
{

    int i;
    char fcd_n = 'N', fcd_t = 'T';
    complex double alpha, beta;
    int nmax, ndim, ione =1, nrow, ncol, n0, n1; 
	int *desca, *descb, *descc, *descd;


    alpha = 1.0;
    beta = 0.0;

    n0 = cei.probe_in_block[iprobe - 1];  /* n0 => block index */

    ndim = ct.block_dim[n0];
    nmax = lcr[iprobe].num_states;

	
    nrow = pmo.mxllda_cond[n0];
    ncol = pmo.mxlocc_lead[iprobe-1];
    n1 = nrow * ncol;

    desca = &pmo.desc_cond_lead[ (n0 + (iprobe - 1) * ct.num_blocks) * DLEN ]; /* (C,L) */
    descb = &pmo.desc_lead[ ( iprobe-1)       * DLEN ];                        /* (L,L) */
    descc = &pmo.desc_cond[ ( n0 + n0 * ct.num_blocks) * DLEN ];               /* (C,C) */
    
    descd = &pmo.desc_lead_cond[ (n0 + (iprobe - 1) * ct.num_blocks) * DLEN ]; /* (C,L) */


	/*    ch(C,L) = ch01(C,L) * green(L,L)  */
    /*    sigma(C,C) = ch1(C,L) * ch(L,C)  */ 
		  

    pzgemm (&fcd_n, &fcd_n, &ndim, &nmax, &nmax, &alpha, ch01, &ione, &ione, desca,
            green, &ione, &ione, descb, &beta, ch, &ione, &ione, desca);

    pzgemm (&fcd_n, &fcd_n, &ndim, &ndim, &nmax, &alpha, ch, &ione, &ione, desca,
            ch10, &ione, &ione, descd, &beta, sigma, &ione, &ione, descc);



}
