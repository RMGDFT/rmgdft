#ifndef RMG_ErrorFuncs_H
#define RMG_ErrorFuncs_H 1

#include "make_conf.h"

#if __cplusplus


#if GPU_ENABLED
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cublas_v2.h>

void RmgCudaError(const char *file, int line, const cudaError_t cudaStatus, const char * errorMessage);
void RmgCudaError(const char *file, int line, const cublasStatus_t status, const char * errorMessage);
void ProcessCublasError(cublasStatus_t custat);
#endif

#endif
#endif
