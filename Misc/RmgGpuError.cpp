
#include "rmg_error.h"
#include "ErrorFuncs.h"


#if CUDA_ENABLED
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cublas_v2.h>

void RmgCudaError(const char *file, int line, const cudaError_t cudaStatus, const char * errorMessage) {

    if (cudaStatus != cudaSuccess) {
        rmg_error_handler(file, line, errorMessage);
    }

}


void RmgCudaError(const char *file, int line, const cublasStatus_t status, const char * errorMessage) {

    if (status != CUBLAS_STATUS_SUCCESS) {
        rmg_error_handler(file, line, errorMessage);
    }

}


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
    else
    {
        printf("UNKNOWN CUBLAS ERROR");
    }

}

#endif




#if HIP_ENABLED

#include <hipblas.h>

void RmgHipError(const char *file, int line, const hipError_t hipStatus, const char * errorMessage)
{
    if (hipStatus != hipSuccess) {
        rmg_error_handler(file, line, errorMessage);
    }
}


void RmgHipError(const char *file, int line, const hipblasStatus_t status, const char * errorMessage)
{
    if (status != HIPBLAS_STATUS_SUCCESS) {
        rmg_error_handler(file, line, errorMessage);
    }
}


void ProcessHipblasError(hipblasStatus_t hipstat)
{
    if(hipstat==HIPBLAS_STATUS_SUCCESS)
        return;

    if(hipstat==HIPBLAS_STATUS_NOT_INITIALIZED)
    {
        printf("'HIPBLAS_STATUS_NOT_INITIALIZED'");
    }
    else if(hipstat==HIPBLAS_STATUS_ALLOC_FAILED)
    {
        printf("HIPBLAS_STATUS_ALLOC_FAILED");
    }
    else if(hipstat==HIPBLAS_STATUS_INVALID_VALUE)
    {
        printf("HIPBLAS_STATUS_INVALID_VALUE");
    }
    else if(hipstat==HIPBLAS_STATUS_ARCH_MISMATCH)
    {
        printf("HIPBLAS_STATUS_ARCH_MISMATCH");
    }
    else if(hipstat==HIPBLAS_STATUS_MAPPING_ERROR)
    {
        printf("HIPBLAS_STATUS_MAPPING_ERROR");
    }
    else if(hipstat==HIPBLAS_STATUS_EXECUTION_FAILED)
    {
        printf("HIPBLAS_STATUS_EXECUTION_FAILED");
    }
    else if(hipstat==HIPBLAS_STATUS_INTERNAL_ERROR)
    {
        printf("HIPBLAS_STATUS_INTERNAL_ERROR");
    }
    else
    {
        printf("UNKNOWN HIPBLAS ERROR");
    }

}

#endif

