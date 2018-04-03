#include <complex>
#include <typeinfo>
#include <string.h>

#include "const.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "RmgGemm.h"
#include "GpuAlloc.h"
#include "ErrorFuncs.h"


#if GPU_ENABLED
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cublas_v2.h>
#endif

#define         dsymm   dsymm_
#define         zsymm   zsymm_


extern "C" {
void dsymm(const char *, const char *, int *, int *, double *, double *, int *, double *, int *, double *, double *, int *);
void zsymm(const char *, const char *, int *, int *, std::complex<double> *, std::complex<double> *, int *, std::complex<double> *, int *, std::complex<double> *, std::complex<double> *, int *);
}




/*
  These functions are used to hide the details of the matrix multiplication data types and GPU 
  utilization from the higher level routines.

  The first 13 arguments are the same as the standard dgemm args but with scalar quantities passed
  by value instead of by reference. The last three arguments are used only when GPU_ENABLED is true.
  In that case if

  [ABC]gpu == NULL   transfer [ABC] to gpu and perform matrix multiplication there. 

  If [ABC]gpu != NULL and Copyto[ABC]gpu is true then data needs to be transferred to the GPU.

  If CopyfromCgpu flag is true and Cgpu!=NULL then copy from Cgpu to C. Otherwise leave
  data in Cgpu for reuse.

*/


template void RmgSymm<double>(char *, char *, int, int, double, double *, int, double *, int, 
                                  double, double *, int);
template void RmgSymm<std::complex<double> >(char *, char *, int, int, std::complex<double>, 
                      std::complex<double> *, int, std::complex<double> *, int, 
                      std::complex<double>, std::complex<double> *, int);

template <typename DataType> void RmgSymm(char *side, char *uplo, int m, int n, DataType alpha,
                             DataType *A, int lda, DataType *B, int ldb, DataType beta, DataType *C, int ldc)
{

#if GPU_ENABLED
    cublasStatus_t custat;
    cublasSideMode_t cu_side = CUBLAS_SIDE_LEFT;
    cublasFillMode_t cu_uplo = CUBLAS_FILL_MODE_LOWER;
    DataType ZERO_t(0.0);

    cudaDeviceSynchronize();

    if(!strcmp(side, "l")) cu_side = CUBLAS_SIDE_LEFT;
    if(!strcmp(side, "L")) cu_side = CUBLAS_SIDE_LEFT;
    if(!strcmp(side, "r")) cu_side = CUBLAS_SIDE_RIGHT;
    if(!strcmp(side, "R")) cu_side = CUBLAS_SIDE_RIGHT;

    if(!strcmp(uplo, "l")) cu_uplo = CUBLAS_FILL_MODE_LOWER;
    if(!strcmp(uplo, "L")) cu_uplo = CUBLAS_FILL_MODE_LOWER;
    if(!strcmp(uplo, "u")) cu_uplo = CUBLAS_FILL_MODE_UPPER;
    if(!strcmp(uplo, "U")) cu_uplo = CUBLAS_FILL_MODE_UPPER;

    int ka = m;
    if(!strcmp("r", side)) ka = n;
    if(!strcmp("R", side)) ka = n;

   if(typeid(DataType) == typeid(std::complex<double>)) {
        custat = cublasZsymm(ct.cublas_handle, cu_side, cu_uplo, m, n,
                            (cuDoubleComplex *)&alpha,
                            (cuDoubleComplex*)A, lda,
                            (cuDoubleComplex*)B, ldb,
                            (cuDoubleComplex*)&beta, (cuDoubleComplex*)C, ldc );
        ProcessCublasError(custat);
        RmgCudaError(__FILE__, __LINE__, custat, "Problem executing cublasZgemm");
    }
    else {
        custat = cublasDsymm(ct.cublas_handle, cu_side, cu_uplo, m, n,
                            (double*)&alpha,
                            (double*)A, lda,
                            (double*)B, ldb,
                            (double*)&beta, (double*)C, ldc );
        ProcessCublasError(custat);
        RmgCudaError(__FILE__, __LINE__, custat, "Problem executing cublasDgemm");
    }

#else

    if(typeid(DataType) == typeid(std::complex<double>)) {
        zsymm(side, uplo, &m, &n, (std::complex<double> *)&alpha, (std::complex<double> *)A, &lda, 
             (std::complex<double> *)B, &ldb, (std::complex<double> *)&beta, (std::complex<double> *)C, &ldc);
    }
    else {
        dsymm(side, uplo, &m, &n, (double *)&alpha, (double *)A, &lda, 
        (double *)B, &ldb, (double *)&beta, (double *)C, &ldc);
    }

#endif
}

