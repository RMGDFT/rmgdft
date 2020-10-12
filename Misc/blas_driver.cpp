#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <complex>


#if CUDA_ENABLED
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cublas_v2.h>
#endif

#include "blas.h"
#include "blacs.h"
#include "Scalapack.h"

#include "typedefs.h"
#include "rmg_control.h"

void my_sync_device()
{
#if CUDA_ENABLED
    DeviceSynchronize();
#endif
}

void zcopy_driver (int n, std::complex<double> *A, int ia, std::complex<double> *B, int ib) 
{

#if CUDA_ENABLED
    cublasZcopy (ct.cublas_handle, n, (cuDoubleComplex *)A, ia, (cuDoubleComplex *)B, ib);
#else
    zcopy (&n, A, &ia, B, &ib);
#endif
}


void zaxpy_driver (int n, std::complex<double> alpha, std::complex<double> *A, int ia, std::complex<double> *B, int ib) 
{

#if CUDA_ENABLED
    cublasZaxpy (ct.cublas_handle, n, (cuDoubleComplex *)&alpha, (cuDoubleComplex *)A, ia, (cuDoubleComplex *)B, ib);
#else
    zaxpy (&n, &alpha, A, &ia, B, &ib);
#endif
}

void dzasum_driver(int n, std::complex<double> *A, int ia, double *sum)
{
#if CUDA_ENABLED
    cublasDzasum (ct.cublas_handle, n, (cuDoubleComplex *)A, ia, sum);
#else
    *sum = dzasum(&n, (double *)A, &ia);
#endif
}


void dcopy_driver (int n, double *A, int ia, double *B, int ib) 
{

#if CUDA_ENABLED
    cublasDcopy (ct.cublas_handle, n, A, ia, B, ib);
#else
    dcopy (&n, A, &ia, B, &ib);
#endif
}


void daxpy_driver (int n, double alpha, double *A, int ia, double *B, int ib) 
{

#if CUDA_ENABLED
    cublasDaxpy (ct.cublas_handle, n, &alpha, A, ia, B, ib);
#else
    daxpy (&n, &alpha, A, &ia, B, &ib);
#endif
}

void dscal_driver(int n, double beta, double *A, int ione)
{
#if CUDA_ENABLED
    cublasDscal (ct.cublas_handle, n, &beta, A, ione);
#else
    dscal(&n, &beta, A, &ione);
#endif

}

void dgemm_driver (char *transa, char *transb, int m, int n, int k, 
double alpha, double *A, int ia, int ja, int *desca,
double *B, int ib, int jb, int *descb, double beta, 
double *C, int ic, int jc, int *descc)
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

#if CUDA_ENABLED

    if(nprow*npcol != 1)
    {
        printf ("GPU ENALBED but nprow*npcol !=1  nprow= %d npcol=%d \n", nprow, npcol);
        fflush (NULL);
        exit (0);
    }


    cublasOperation_t cu_transA = CUBLAS_OP_N, cu_transB = CUBLAS_OP_N;

    if(!strcmp(transa, "t")) cu_transA = CUBLAS_OP_T;
    if(!strcmp(transa, "T")) cu_transA = CUBLAS_OP_T;

    if(!strcmp(transb, "t")) cu_transB = CUBLAS_OP_T;
    if(!strcmp(transb, "T")) cu_transB = CUBLAS_OP_T;

    cublasDgemm (ct.cublas_handle, cu_transA, cu_transB, m, n, k, 
            &alpha, A, lda, B, ldb, &beta, C, ldc );

#else
    //  use scalapack if nprow * npcol > 1
    if(nprow*npcol > 1)  
    {
        pdgemm (transa, transb, &m, &n, &k, &alpha, A, &ia, &ja, desca,
                B, &ib, &jb, descb, &beta, C, &ic, &jc, descc);
    }
    else
    {

        dgemm(transa, transb, &m, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc);
    }
#endif


}


void zgemm_driver (char *transa, char *transb, int m, int n, int k, 
std::complex<double> alpha, std::complex<double> *A, int ia, int ja, int *desca,
std::complex<double> *B, int ib, int jb, int *descb, std::complex<double> beta, 
std::complex<double> *C, int ic, int jc, int *descc)
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

#if CUDA_ENABLED

    if(nprow*npcol != 1)
    {
        printf ("GPU ENALBED but nprow*npcol !=1  nprow= %d npcol=%d \n", nprow, npcol);
        fflush (NULL);
        exit (0);
    }


    cublasOperation_t cu_transA = CUBLAS_OP_N, cu_transB = CUBLAS_OP_N;

    if(!strcmp(transa, "t")) cu_transA = CUBLAS_OP_T;
    if(!strcmp(transa, "T")) cu_transA = CUBLAS_OP_T;
    if(!strcmp(transa, "c")) cu_transA = CUBLAS_OP_C;
    if(!strcmp(transa, "C")) cu_transA = CUBLAS_OP_C;

    if(!strcmp(transb, "t")) cu_transB = CUBLAS_OP_T;
    if(!strcmp(transb, "T")) cu_transB = CUBLAS_OP_T;
    if(!strcmp(transb, "c")) cu_transB = CUBLAS_OP_C;
    if(!strcmp(transb, "C")) cu_transB = CUBLAS_OP_C;

    cublasZgemm (ct.cublas_handle, cu_transA, cu_transB, m, n, k, 
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
