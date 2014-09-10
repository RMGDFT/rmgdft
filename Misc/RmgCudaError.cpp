#include "make_conf.h"
#include "rmg_error.h"
#include "ErrorFuncs.h"


#if GPU_ENABLED
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

#endif
