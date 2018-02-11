/************************** SVN Revision Information **************************
 **    $Id$    **
 ******************************************************************************/


#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include "make_conf.h"

#if GPU_ENABLED
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cublas_v2.h>
#endif

#include "main.h"
#include "Scalapack.h"
#include "blas.h"



void zgemm_driver (char *transa, char *transb, int m, int n, int k, 
double complex alpha, double complex *A, int ia, int ja, int *desca,
double complex *B, int ib, int jb, int *descb, double complex beta, 
double complex *C, int ic, int jc, int *descc)
{

    int nprow, npcol, myrow, mycol;
    int lda=desca[8], ldb=descb[8], ldc = descc[8];
    int ictxt = desca[1];

    Cblacs_gridinfo (ictxt, &nprow, &npcol, &myrow, &mycol);
    if(nprow*npcol <1) 
    {
        printf ("error in zgemmdriver nprow= %d npcol=%d \n", nprow, npcol);
        fflush (NULL);
        exit (0);
    }

#if GPU_ENABLED

    if(nprow*npcol != 1)
    {
        printf ("GPU ENALBED but nprow*npcol !=1  nprow= %d npcol=%d \n", nprow, npcol);
        fflush (NULL);
        exit (0);
    }


    cublasStatus_t custat;
    cublasOperation_t cu_transA = CUBLAS_OP_N, cu_transB = CUBLAS_OP_N;

    if(!strcmp(transa, "t")) cu_transA = CUBLAS_OP_T;
    if(!strcmp(transa, "T")) cu_transA = CUBLAS_OP_T;
    if(!strcmp(transa, "c")) cu_transA = CUBLAS_OP_C;
    if(!strcmp(transa, "C")) cu_transA = CUBLAS_OP_C;

    if(!strcmp(transb, "t")) cu_transB = CUBLAS_OP_T;
    if(!strcmp(transb, "T")) cu_transB = CUBLAS_OP_T;
    if(!strcmp(transb, "c")) cu_transB = CUBLAS_OP_C;
    if(!strcmp(transb, "C")) cu_transB = CUBLAS_OP_C;

    custat = cublasZgemm (ct.cublas_handle, cu_transA, cu_transB, m, n, k, 
            (cuDoubleComplex *)&alpha,
            (cuDoubleComplex*)A, lda,
            (cuDoubleComplex*)B, ldb,
            (cuDoubleComplex*)&beta, (cuDoubleComplex*)C, ldc );

#else
    //  use scalapack if nprow * npcol > 1
    if(nprow*npcol > 1)  
    {
        pzgemm (transa, transb, &m, &n, &k, &alpha, A, &ia, &ja, desca,
                B, &ib, &jb, descb, &beta, C, &ic, &jc, descc);
    }
    else
    {

        zgemm(transa, transb, &m, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc);
    }
#endif


}
