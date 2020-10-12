#include <stdio.h>
#include "ErrorFuncs.h"


#if CUDA_ENABLED
#include <cublas_v2.h>
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

