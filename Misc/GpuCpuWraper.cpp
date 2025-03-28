#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <complex>
#include "GpuAlloc.h"
#include "rmgtypedefs.h"
#include "transition.h"


#if CUDA_ENABLED
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cublas_v2.h>
#endif


void RmgMemcpy(void *A_dest, void *A_source, size_t a_size)
{
    if(A_dest == A_source) return;
#if CUDA_ENABLED
        cudaMemcpy(A_dest, A_source, a_size, cudaMemcpyDefault);
#elif SYCL_ENABLED
        gpuMemcpy(A_dest, A_source, a_size, gpuMemcpyHostToDevice);
#elif HIP_ENABLED
        hipMemcpy(A_dest, A_source, a_size, hipMemcpyDefault);
#else
       memcpy(A_dest, A_source, a_size); 
#endif 

} 

void MemcpyHostDevice (size_t a_size, void *A_host, void *A_device)
{
#if !CUDA_ENABLED && !SYCL_ENABLED && !HIP_ENABLED
    return;
#endif
    if(ct.tddft_gpu)
    {
#if CUDA_ENABLED
        cudaMemcpy(A_device, A_host, a_size, cudaMemcpyDefault);
#elif SYCL_ENABLED
        gpuMemcpy(A_device, A_host, a_size, gpuMemcpyHostToDevice);
#elif HIP_ENABLED
        hipMemcpyHtoD(A_device, A_host, a_size );
#endif 
    }
    else
    {
        memcpy(A_device, A_host, a_size);
    }
}

void MemcpyDeviceHost (size_t a_size, void *A_device, void *A_host)
{
#if !CUDA_ENABLED && !SYCL_ENABLED && !HIP_ENABLED
    return;
#endif

    if(ct.tddft_gpu)
    {
#if CUDA_ENABLED
    cudaMemcpy(A_host, A_device, a_size, cudaMemcpyDefault);
#elif SYCL_ENABLED
    gpuMemcpy(A_host, A_device, a_size, gpuMemcpyDeviceToHost);
#elif HIP_ENABLED
    hipMemcpyDtoH(A_host, A_device, a_size );
#endif 
    }
    else
    {
        memcpy(A_host, A_device, a_size);
    }
}


template double *MemoryPtrHostDevice(double *ptr_host, double *ptr_device);
template std::complex<double> *MemoryPtrHostDevice(std::complex<double> *ptr_host, std::complex<double> *ptr_device);
template <typename T> T *MemoryPtrHostDevice(T *ptr_host, T *ptr_device)
{
    if(ct.tddft_gpu)
    {
#if CUDA_ENABLED || HIP_ENABLED || SYCL_ENABLED
    return ptr_device;
#endif 
    }
    return ptr_host;
}


