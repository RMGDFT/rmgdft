#include "GpuAlloc.h"

/************************** SVN Revision Information **************************
 **    $Id$    **
 ******************************************************************************/


#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <complex>


#if CUDA_ENABLED
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cublas_v2.h>
#endif



void MemcpyHostDevice (size_t a_size, void *A_host, void *A_device)
{

#if CUDA_ENABLED
  cudaMemcpy(A_device, A_host, a_size, cudaMemcpyDefault);
#elif HIP_ENABLED
    hipMemcpyHtoD(A_device, A_host, a_size );
#endif 
}

void MemcpyDeviceHost (size_t a_size, void *A_device, void *A_host)
{

#if CUDA_ENABLED
  cudaMemcpy(A_host, A_device, a_size, cudaMemcpyDefault);
#elif HIP_ENABLED
    hipMemcpyDtoH(A_host, A_device, a_size );
#endif 
}


template double *MemoryPtrHostDevice(double *ptr_host, double *ptr_device);
template std::complex<double> *MemoryPtrHostDevice(std::complex<double> *ptr_host, std::complex<double> *ptr_device);
template <typename T> T *MemoryPtrHostDevice(T *ptr_host, T *ptr_device)
{
#if CUDA_ENABLED || HIP_ENABLED
    return ptr_device;
#else
    return ptr_host;
#endif 
}


