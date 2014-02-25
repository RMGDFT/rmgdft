/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include "main.h"
#include "init_var_negf.h"
#include "LCR.h"
#include "pmo.h"



void rho_munu_cuda (complex double * rho_mn, complex double * green_C, complex double * gamma, int iprobe)
{

#if GPU_ENABLED
    cuDoubleComplex cuone,  cuzero;
    cublasOperation_t transC = CUBLAS_OP_C, transN = CUBLAS_OP_N;

    int i, j, n_green, n1, n2;

    int N, *ni, nrow, ncol, nL, N1;
    int maxrow, maxcol, ione =1;

    N = ct.num_blocks;
    ni = ct.block_dim;


    N1 = cei.probe_in_block[iprobe - 1];
    nL = ct.block_dim[N1];
    nrow = pmo.mxllda_cond[N1];
    ncol = pmo.mxlocc_cond[N1];


    cuone.x = 1.0;
    cuone.y = 0.0;

    cuzero.x = 0.0;
    cuzero.y = 0.0;

    maxrow = 0;
    maxcol = 0;
    for (i = 0; i < N; i++)
    {
        maxrow = max(maxrow, pmo.mxllda_cond[i]);
        maxcol = max(maxcol, pmo.mxlocc_cond[i]);
    }

    n_green = 0;

    n1 = maxrow * maxcol;
    cublasSetVector( nL * nL, sizeof( complex double ), gamma, ione, ct.gpu_Gii, ione );

    for (i = 0; i < N - 1; i++)
    {
        n1 = ni[i];
        n2 = ni[i + 1];

        /*  temp = G_i0 * Gamma  */
        cublasZgemm(ct.cublas_handle, transN, transN, n1, nL, nL, &cuone, &ct.gpu_Gtem[n_green], n1,
                ct.gpu_Gii, nL, &cuzero, ct.gpu_temp, n1);

        /* rho_mn (i,i) = temp * G_i0^, the block (i,i) */

        cublasZgemm(ct.cublas_handle, transN, transC, n1, n1, nL, &cuone, ct.gpu_temp, n1,
                &ct.gpu_Gtem[n_green], n1, &cuzero, &ct.gpu_Gtri[pmo.diag_begin[i]], n1);

        /* rho_mn (i,i+1) = temp * G_i+10^, the block (i,i) */
        n_green += pmo.mxllda_cond[i] * maxcol;


        cublasZgemm(ct.cublas_handle, transN, transC, n1, n2, nL, &cuone, ct.gpu_temp, n1,
                &ct.gpu_Gtem[n_green], n2, &cuzero, &ct.gpu_Gtri[pmo.offdiag_begin[i]], n1);


    }

    /* calculate the last block  */
    n1 = ni[N - 1];

    cublasZgemm(ct.cublas_handle, transN, transN, n1, nL, nL, &cuone, &ct.gpu_Gtem[n_green], n1,
            ct.gpu_Gii, nL, &cuzero, ct.gpu_temp, n1);
    cublasZgemm(ct.cublas_handle, transN, transC, n1, n1, nL, &cuone, ct.gpu_temp, n1,
            &ct.gpu_Gtem[n_green], n1, &cuzero, &ct.gpu_Gtri[pmo.diag_begin[N-1]], n1);


    cublasGetVector(pmo.ntot, sizeof( complex double ), ct.gpu_Gtri, ione, rho_mn, ione );



#endif
}
