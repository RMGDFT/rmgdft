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
#endif

#define         dgemm   dgemm_
#define         zgemm   zgemm_


extern "C" {
void dgemm(const char *, const char *, int *, int *, int *, double *, double *, int *, double *, int *, double *, double *, int *);
void zgemm(const char *, const char *, int *, int *, int *, std::complex<double> *, std::complex<double> *, int *, std::complex<double> *, int *, std::complex<double> *, std::complex<double> *, int *);
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


template void RmgGemm<double>(char *, char *, int, int, int, double, double *, int, double *, int, 
                                  double, double *, int, double *, double *, double *, bool, bool, bool, bool);
template void RmgGemm<std::complex<double> >(char *, char *, int, int, int, std::complex<double>, 
                      std::complex<double> *, int, std::complex<double> *, int, 
                      std::complex<double>, std::complex<double> *, int, std::complex<double> *, 
                      std::complex<double> *, std::complex<double> *, bool, bool, bool, bool);

template <typename DataType> void RmgGemm(char *transa, char *transb, int m, int n, int k, 
                             DataType alpha, DataType *A, int lda, DataType *B, int ldb, DataType beta, 
                             DataType *C, int ldc, DataType *Agpu, DataType *Bgpu, 
                             DataType *Cgpu, bool CopytoAgpu, bool CopytoBgpu, 
                             bool CopytoCgpu, bool CopyfromCgpu )
{

#if GPU_ENABLED

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

    if((m <= 512) || (k <= 512)) {
//        cublasXtSetCpuRatio(ct.cublasXt_handle, CUBLASXT_GEMM, CUBLASXT_DOUBLE, 1.0);
    }
    else {
//        cublasXtSetCpuRatio(ct.cublasXt_handle, CUBLASXT_GEMM, CUBLASXT_DOUBLE, 0.0);
    }
    cudaPointerAttributes attrib_A, attrib_B, attrib_C;
    cudaError_t cudaErrA, cudaErrB, cudaErrC;

    cudaErrA = cudaPointerGetAttributes(&attrib_A, A);
    cudaErrB = cudaPointerGetAttributes(&attrib_B, B);
    cudaErrC = cudaPointerGetAttributes(&attrib_C, C);
    bool UsingCudaMemory = ((cudaErrA == cudaSuccess) && (cudaErrB == cudaSuccess) && (cudaErrC == cudaSuccess));
//if(!UsingCudaMemory)printf("NOT USING\n");
    // Check attributes and if all matrices are located on the host then
    // we can use the XT interface.
    if((attrib_A.memoryType == cudaMemoryTypeHost) &&
       (attrib_B.memoryType == cudaMemoryTypeHost) &&
       (attrib_C.memoryType == cudaMemoryTypeHost) &&
       (Agpu == NULL) &&
       (Bgpu == NULL) &&
       (Cgpu == NULL) &&
        UsingCudaMemory) {
//printf("HOST GEMM\n");

            if(typeid(DataType) == typeid(std::complex<double>)) {
                custat = cublasXtZgemm(ct.cublasXt_handle, cu_transA, cu_transB, m, n, k,
                                    (cuDoubleComplex *)&alpha,
                                    (cuDoubleComplex*)A, lda,
                                    (cuDoubleComplex*)B, ldb,
                                    (cuDoubleComplex*)&beta, (cuDoubleComplex*)C, ldc );
                ProcessCublasError(custat);
                RmgCudaError(__FILE__, __LINE__, custat, "Problem executing cublasXtZgemm");
            }
            else {
                custat = cublasXtDgemm(ct.cublasXt_handle, cu_transA, cu_transB, m, n, k,
                                    (double*)&alpha,
                                    (double*)A, lda,
                                    (double*)B, ldb,
                                    (double*)&beta, (double*)C, ldc );
                ProcessCublasError(custat);
                RmgCudaError(__FILE__, __LINE__, custat, "Problem executing cublasXtDgemm");
            }
            return;

    }

    // If all matrices are located on the device then we can use the regular interface
    if((attrib_A.memoryType == cudaMemoryTypeDevice) &&
       (attrib_B.memoryType == cudaMemoryTypeDevice) &&
       (attrib_C.memoryType == cudaMemoryTypeDevice) &&
        UsingCudaMemory) {
//printf("DEVICE GEMM\n");
            if(typeid(DataType) == typeid(std::complex<double>)) {
                custat = cublasZgemm(ct.cublas_handle, cu_transA, cu_transB, m, n, k,
                                    (cuDoubleComplex *)&alpha,
                                    (cuDoubleComplex*)A, lda,
                                    (cuDoubleComplex*)B, ldb,
                                    (cuDoubleComplex*)&beta, (cuDoubleComplex*)C, ldc );
                ProcessCublasError(custat);
                RmgCudaError(__FILE__, __LINE__, custat, "Problem executing cublasZgemm");
            }
            else {
                custat = cublasDgemm(ct.cublas_handle, cu_transA, cu_transB, m, n, k,
                                    (double*)&alpha,
                                    (double*)A, lda,
                                    (double*)B, ldb,
                                    (double*)&beta, (double*)C, ldc );
                ProcessCublasError(custat);
                RmgCudaError(__FILE__, __LINE__, custat, "Problem executing cublasDgemm");
            }
            return;
    }


//printf("MIXED GEMM\n");


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

        Bgpu1 = (DataType *)GpuMallocDevice(kb * ldb * sizeof( DataType ));
        custat = cublasSetVector(kb * ldb , sizeof( DataType ), B, 1, Bgpu1, 1 );
        RmgCudaError(__FILE__, __LINE__, custat, "Problem transferring B matrix to GPU.");

    }
    else {
        Bgpu1 = Bgpu;
        if(CopytoBgpu) {
            custat = cublasSetVector(kb * ldb , sizeof( DataType ), B, 1, Bgpu1, 1 );
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
        custat = cublasZgemm(ct.cublas_handle, cu_transA, cu_transB, m, n, k,
                            (cuDoubleComplex *)&alpha,
                            (cuDoubleComplex*)Agpu1, lda,
                            (cuDoubleComplex*)Bgpu1, ldb,
                            (cuDoubleComplex*)&beta, (cuDoubleComplex*)Cgpu1, ldc );
        ProcessCublasError(custat);
        RmgCudaError(__FILE__, __LINE__, custat, "Problem executing cublasZgemm");
    }
    else {
        custat = cublasDgemm(ct.cublas_handle, cu_transA, cu_transB, m, n, k,
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
        zgemm(transa, transb, &m, &n, &k, (std::complex<double> *)&alpha, (std::complex<double> *)A, &lda, 
             (std::complex<double> *)B, &ldb, (std::complex<double> *)&beta, (std::complex<double> *)C, &ldc);
    }
    else {
        dgemm(transa, transb, &m, &n, &k, (double *)&alpha, (double *)A, &lda, 
        (double *)B, &ldb, (double *)&beta, (double *)C, &ldc);
    }

#endif
}

#if GPU_ENABLED
void ProcessCublasError(cublasStatus_t custat)
{
    if(custat==CUBLAS_STATUS_SUCCESS)
        return;

    if(custat==CUBLAS_STATUS_NOT_INITIALIZED)
    {
        printf("'CUBLAS_STATUS_NOT_INITIALIZED'");
    }
    else if(custat==CUBLAS_STATUS_ALLOC_FAILED)
    {
        printf("CUBLAS_STATUS_ALLOC_FAILED");
    }
    else if(custat==CUBLAS_STATUS_INVALID_VALUE)
    {
        printf("CUBLAS_STATUS_INVALID_VALUE");
    }
    else if(custat==CUBLAS_STATUS_ARCH_MISMATCH)
    {
        printf("CUBLAS_STATUS_ARCH_MISMATCH");
    }
    else if(custat==CUBLAS_STATUS_MAPPING_ERROR)
    {
        printf("CUBLAS_STATUS_MAPPING_ERROR");
    }
    else if(custat==CUBLAS_STATUS_EXECUTION_FAILED)
    {
        printf("CUBLAS_STATUS_EXECUTION_FAILED");
    }
    else if(custat==CUBLAS_STATUS_INTERNAL_ERROR)
    {
        printf("CUBLAS_STATUS_INTERNAL_ERROR");
    }

    printf("UNKNOWN CUBLAS ERROR");

}
#endif
