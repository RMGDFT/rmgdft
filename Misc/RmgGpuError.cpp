
#include "rmg_error.h"
#include "ErrorFuncs.h"


#if CUDA_ENABLED
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cublas_v2.h>

void RmgGpuError(const char *file, int line, const cudaError_t cudaStatus, const char * errorMessage) {

    if (cudaStatus != cudaSuccess) {
        rmg::error(errorMessage);
    }

}


void RmgGpuError(const char *file, int line, const cublasStatus_t status, const char * errorMessage) {

    if (status != CUBLAS_STATUS_SUCCESS) {
        rmg::error(errorMessage);
    }

}


#endif




#if HIP_ENABLED

#include <hipblas/hipblas.h>

void RmgGpuError(const char *file, int line, const hipError_t hipStatus, const char * errorMessage)
{
    if (hipStatus != hipSuccess) {
        rmg::error(errorMessage);
    }
}


void RmgGpuError(const char *file, int line, const hipblasStatus_t status, const char * errorMessage)
{
    if (status != HIPBLAS_STATUS_SUCCESS) {
        rmg::error(errorMessage);
    }
}


#elif SYCL_ENABLED
#include "rmgtypedefs.h"
void RmgGpuError(const char *file, int line, const gpuError_t hipStatus, const char * errorMessage)
{
//    if (hipStatus != hipSuccess) {
//        rmg::error(errorMessage);
//    }
}


//void RmgGpuError(const char *file, int line, const hipblasStatus_t status, const char * errorMessage)
//{
//    if (status != HIPBLAS_STATUS_SUCCESS) {
//        rmg::error(errorMessage);
//    }
//}

#endif

