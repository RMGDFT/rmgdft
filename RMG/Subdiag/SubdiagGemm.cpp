#include <complex>
#include "const.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "rmg_alloc.h"
#include "rmg_error.h"
#include "Subdiag.h"
#include "GpuAlloc.h"


#if GPU_ENABLED
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cublas_v2.h>
#endif

#define         dgemm   dgemm_
#define         zgemm   zgemm_


extern "C" {
void dgemm(const char *, const char *, int *, int *, int *, double *, double *, int *, double *, int *, double *, double *, int *);
void zgemm(const char *, const char *, int *, int *, int *, std::complex<double> *, std::complex<double> *, int *, std::complex<double> *, int *, std::complex<double> *, std::complex<double> *, int *);
}




/*
  These functions are used to simplify the subspace diagonalization code by hiding the
  details of the matrix multiplication data type and GPU utilization from the higher level routines.

  The first 13 arguments are the same as the standard dgemm args but with scalar quantities passed
  by value instead of by reference. The last three arguments are used only when GPU_ENABLED is true.
  In that case if

  [ABC]gpu == NULL   transfer [ABC] to gpu and perform matrix multiplication there. 

  If [ABC]gpu != NULL then data does not need to be transferred to the GPU.

  Copy C gpu data back to C.

*/


template void SubdiagGemm<double>(char *, char *, int, int, int, double, double *, int, double *, int, 
                                  double, double *, int, double *, double *, double *);
template void SubdiagGemm<std::complex<double> >(char *, char *, int, int, int, std::complex<double>, std::complex<double> *, int, std::complex<double> *, int, 
                                  std::complex<double>, std::complex<double> *, int, std::complex<double> *, std::complex<double> *, std::complex<double> *);

template <typename DataType> void SubdiagGemm(char *transa, char *transb, int m, int n, int k, 
                             DataType alpha, DataType *A, int lda, DataType *B, int ldb, DataType beta, DataType *C, int ldc,
                             DataType *Agpu, DataType *Bgpu, DataType *Cgpu )
{

#if GPU_ENABLED
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

    int ka = m;
    if(!strcmp("n", transa)) ka = k;
    if(!strcmp("N", transa)) ka = k;

    int kb = k;
    if(!strcmp("n", transb)) kb = n;
    if(!strcmp("N", transb)) kb = n;


    DataType *Agpu1;
    DataType *Bgpu1;
    DataType *Cgpu1;

    if(Agpu == NULL) {

        Agpu1 = (DataType *)GpuMalloc(ka * lda * sizeof( DataType ));
        custat = cublasSetVector(ka * lda , sizeof( DataType ), A, 1, Agpu1, 1 );

    }
    else {
        Agpu1 = Agpu;
    }
    if(Bgpu == NULL) {

        Bgpu1 = (DataType *)GpuMalloc(kb * ldb * sizeof( DataType ));
        custat = cublasSetVector(kb * ldb , sizeof( DataType ), B, 1, Bgpu1, 1 );

    }
    else {
        Bgpu1 = Bgpu;
    }
    if(Cgpu == NULL) {

        Cgpu1 = (DataType *)GpuMalloc(n * ldc * sizeof( DataType ));
        custat = cublasSetVector(n * ldc , sizeof( DataType ), C, 1, Cgpu1, 1 );

    }
    else {
        Cgpu1 = Cgpu;
    }

    if(typeid(DataType) == typeid(std::complex<double>)) {
        custat = cublasZgemm(ct.cublas_handle, cu_transA, cu_transB, m, n, k,
                            (cuDoubleComplex *)&alpha,
                            (cuDoubleComplex*)Agpu1, lda,
                            (cuDoubleComplex*)Bgpu1, ldb,
                            (cuDoubleComplex*)&beta, (cuDoubleComplex*)Cgpu1, ldc );
    }
    else {
        custat = cublasDgemm(ct.cublas_handle, cu_transA, cu_transB, m, n, k,
                            (double*)&alpha,
                            (double*)Agpu1, lda,
                            (double*)Bgpu1, ldb,
                            (double*)&beta, (double*)Cgpu1, ldc );
    }

    // Retreive data from the GPU.
    custat = cublasGetVector(n * ldc, sizeof( DataType ), Cgpu1, 1, C, 1 );

    if(Cgpu == NULL) GpuFree(Cgpu1);
    if(Bgpu == NULL) GpuFree(Bgpu1);
    if(Agpu == NULL) GpuFree(Agpu1);

#else

    if(typeid(DataType) == typeid(std::complex<double>)) {
        zgemm(transa, transb, &m, &n, &k, (std::complex<double> *)&alpha, (std::complex<double> *)A, &lda, 
             (std::complex<double> *)B, &ldb, (std::complex<double> *)&beta, (std::complex<double> *)C, &ldc);
    }
    else {
        dgemm(transa, transb, &m, &n, &k, (double *)&alpha, (double *)A, &lda, 
        (double *)B, &ldb, (double *)&beta, (double *)C, &ldc);
    }

#endif
}
