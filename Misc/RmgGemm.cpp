#include <complex>
#include <typeinfo>
#include <string.h>

#include "const.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "RmgGemm.h"
#include "GpuAlloc.h"
#include "ErrorFuncs.h"
#include "RmgTimer.h"
#include "transition.h"


#if CUDA_ENABLED
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cublas_v2.h>
#endif

#define         dgemm           RMG_FC_GLOBAL(dgemm, DGEMM)
#define         zgemm           RMG_FC_GLOBAL(zgemm, ZGEMM)

extern "C" {
void dgemm(const char *, const char *, int *, int *, int *, double *, double *, int *, double *, int *, double *, double *, int *);
void zgemm(const char *, const char *, int *, int *, int *, std::complex<double> *, std::complex<double> *, int *, std::complex<double> *, int *, std::complex<double> *, std::complex<double> *, int *);
}


/*
  These functions are used to hide the details of the matrix multiplication data types and GPU 
  utilization from the higher level routines.

  The first 13 arguments are the same as the standard dgemm args but with scalar quantities passed
  by value instead of by reference.

*/


template void RmgGemm<double>(char *, char *, int, int, int, double, double *, int, double *, int, 
                                  double, double *, int);

template void RmgGemm<std::complex<double> >(char *, char *, int, int, int, std::complex<double>, 
                      std::complex<double> *, int, std::complex<double> *, int, 
                      std::complex<double>, std::complex<double> *, int);


template <typename DataType> void RmgGemm(char *transa, char *transb, int m, int n, int k, 
                             DataType alpha, DataType *A, int lda, DataType *B, int ldb, DataType beta, 
                             DataType *C, int ldc)
{

#if BLAS_PROFILE
    if(typeid(DataType) == typeid(std::complex<double>))
    {
        if(pct.gridpe==0) printf("ZGEMM CALL m=%d n=%d k=%d\n",m,n,k);
    }
    else
    {
        if(pct.gridpe==0) printf("DGEMM CALL m=%d n=%d k=%d\n",m,n,k);
    }
#endif

#if CUDA_ENABLED

    cublasStatus_t custat;
    cublasOperation_t cu_transA = CUBLAS_OP_N, cu_transB = CUBLAS_OP_N;
    DataType ZERO_t(0.0);

    if(!strcmp(transa, "t")) cu_transA = CUBLAS_OP_T;
    if(!strcmp(transa, "T")) cu_transA = CUBLAS_OP_T;
    if(!strcmp(transa, "c")) cu_transA = CUBLAS_OP_C;
    if(!strcmp(transa, "C")) cu_transA = CUBLAS_OP_C;

    if(!strcmp(transb, "t")) cu_transB = CUBLAS_OP_T;
    if(!strcmp(transb, "T")) cu_transB = CUBLAS_OP_T;
    if(!strcmp(transb, "c")) cu_transB = CUBLAS_OP_C;
    if(!strcmp(transb, "C")) cu_transB = CUBLAS_OP_C;

    int ka = m;
    if(!strcmp("n", transa)) ka = k;
    if(!strcmp("N", transa)) ka = k;

    int kb = k;
    if(!strcmp("n", transb)) kb = n;
    if(!strcmp("N", transb)) kb = n;

    DeviceSynchronize();
    if(typeid(DataType) == typeid(std::complex<double>)) {
        custat = cublasZgemm(ct.cublas_handle, cu_transA, cu_transB, m, n, k,
                            (cuDoubleComplex *)&alpha,
                            (cuDoubleComplex*)A, lda,
                            (cuDoubleComplex*)B, ldb,
                            (cuDoubleComplex*)&beta, (cuDoubleComplex*)C, ldc );
        ProcessCublasError(custat);
        RmgGpuError(__FILE__, __LINE__, custat, "Problem executing cublasZgemm");
    }
    else {
        custat = cublasDgemm(ct.cublas_handle, cu_transA, cu_transB, m, n, k,
                            (double*)&alpha,
                            (double*)A, lda,
                            (double*)B, ldb,
                            (double*)&beta, (double*)C, ldc );
        ProcessCublasError(custat);
        RmgGpuError(__FILE__, __LINE__, custat, "Problem executing cublasDgemm");
    }
    DeviceSynchronize();
    return;

#else

    
//    RmgTimer *RT = new RmgTimer("gemmmmm ");
    if(typeid(DataType) == typeid(std::complex<double>))
    {
        if(ct.use_alt_zgemm)
            MyZgemm(transa, transb, m, n, k, (std::complex<double> *)(&alpha), (std::complex<double> *)A, lda, 
             (std::complex<double> *)B, ldb, (std::complex<double> *)(&beta), (std::complex<double> *)C, ldc);
        else
            zgemm(transa, transb, &m, &n, &k, (std::complex<double> *)&alpha, (std::complex<double> *)A, &lda,
            (std::complex<double> *)B, &ldb, (std::complex<double> *)&beta, (std::complex<double> *)C, &ldc);
    }
    else {
        dgemm(transa, transb, &m, &n, &k, (double *)&alpha, (double *)A, &lda, 
        (double *)B, &ldb, (double *)&beta, (double *)C, &ldc);
    }
//    delete RT;

#endif
}
