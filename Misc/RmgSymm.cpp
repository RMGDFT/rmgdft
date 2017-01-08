#include <complex>
#include <typeinfo>
#include <string.h>
#include "make_conf.h"
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
void ProcessCublasError(cublasStatus_t custat);
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
                                  double, double *, int, double *, double *, double *, bool, bool, bool, bool);
template void RmgSymm<std::complex<double> >(char *, char *, int, int, std::complex<double>, 
                      std::complex<double> *, int, std::complex<double> *, int, 
                      std::complex<double>, std::complex<double> *, int, std::complex<double> *, 
                      std::complex<double> *, std::complex<double> *, bool, bool, bool, bool);

template <typename DataType> void RmgSymm(char *side, char *uplo, int m, int n,
                             DataType alpha, DataType *A, int lda, DataType *B, int ldb, DataType beta, 
                             DataType *C, int ldc, DataType *Agpu, DataType *Bgpu, 
                             DataType *Cgpu, bool CopytoAgpu, bool CopytoBgpu, 
                             bool CopytoCgpu, bool CopyfromCgpu )
{

#if GPU_ENABLED
    cublasStatus_t custat;
    cublasSideMode_t cu_side = CUBLAS_SIDE_LEFT;
    cublasFillMode_t cu_uplo = CUBLAS_FILL_MODE_LOWER;
    DataType ZERO_t(0.0);

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

    DataType *Agpu1;
    DataType *Bgpu1;
    DataType *Cgpu1;

    if(Agpu == NULL) {

        Agpu1 = (DataType *)GpuMallocDevice(ka * lda * sizeof( DataType ));
        custat = cublasSetVector(ka * lda , sizeof( DataType ), A, 1, Agpu1, 1 );
        RmgCudaError(__FILE__, __LINE__, custat, "Problem transferring A matrix to GPU.");

    }
    else {
        Agpu1 = Agpu;
        if(CopytoAgpu) {
            custat = cublasSetVector(ka * lda , sizeof( DataType ), A, 1, Agpu1, 1 );
            RmgCudaError(__FILE__, __LINE__, custat, "Problem transferring A matrix to GPU.");
        }
    }
    if(Bgpu == NULL) {

        Bgpu1 = (DataType *)GpuMallocDevice(n * ldb * sizeof( DataType ));
        custat = cublasSetVector(n * ldb , sizeof( DataType ), B, 1, Bgpu1, 1 );
        RmgCudaError(__FILE__, __LINE__, custat, "Problem transferring B matrix to GPU.");

    }
    else {
        Bgpu1 = Bgpu;
        if(CopytoBgpu) {
            custat = cublasSetVector(n * ldb , sizeof( DataType ), B, 1, Bgpu1, 1 );
            RmgCudaError(__FILE__, __LINE__, custat, "Problem transferring A matrix to GPU.");
        }
    }
    if(Cgpu == NULL) {

        Cgpu1 = (DataType *)GpuMallocDevice(n * ldc * sizeof( DataType ));
        // No need to copy if beta is zero
        if(beta != ZERO_t) {
            custat = cublasSetVector(n * ldc , sizeof( DataType ), C, 1, Cgpu1, 1 );
            RmgCudaError(__FILE__, __LINE__, custat, "Problem transferring C matrix to GPU.");
        }

    }
    else {
        Cgpu1 = Cgpu;
        if(CopytoCgpu) {
            custat = cublasSetVector(n * ldc , sizeof( DataType ), C, 1, Cgpu1, 1 );
            RmgCudaError(__FILE__, __LINE__, custat, "Problem transferring C matrix to GPU.");
        }
    }

    if(typeid(DataType) == typeid(std::complex<double>)) {
        custat = cublasZsymm(ct.cublas_handle, cu_side, cu_uplo, m, n,
                            (cuDoubleComplex *)&alpha,
                            (cuDoubleComplex*)Agpu1, lda,
                            (cuDoubleComplex*)Bgpu1, ldb,
                            (cuDoubleComplex*)&beta, (cuDoubleComplex*)Cgpu1, ldc );
        ProcessCublasError(custat);
        RmgCudaError(__FILE__, __LINE__, custat, "Problem executing cublasZgemm");
    }
    else {
        custat = cublasDsymm(ct.cublas_handle, cu_side, cu_uplo, m, n,
                            (double*)&alpha,
                            (double*)Agpu1, lda,
                            (double*)Bgpu1, ldb,
                            (double*)&beta, (double*)Cgpu1, ldc );
        ProcessCublasError(custat);
        RmgCudaError(__FILE__, __LINE__, custat, "Problem executing cublasDgemm");
    }

    // Retreive data from the GPU.
    if((Cgpu == NULL) || ((Cgpu != NULL) && CopyfromCgpu)) {
        custat = cublasGetVector(n * ldc, sizeof( DataType ), Cgpu1, 1, C, 1 );
        RmgCudaError(__FILE__, __LINE__, custat, "Problem transferring C matrix from GPU to system memory.");
    }

    if(Cgpu == NULL) GpuFreeDevice(Cgpu1);
    if(Bgpu == NULL) GpuFreeDevice(Bgpu1);
    if(Agpu == NULL) GpuFreeDevice(Agpu1);

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

