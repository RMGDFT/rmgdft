#include <stdio.h>
#include "ErrorFuncs.h"


#if HIP_ENABLED
#include <hipblas.h>
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

