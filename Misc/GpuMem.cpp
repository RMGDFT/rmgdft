/*
 *
 * Copyright 2020 The RMG Project Developers. See the COPYRIGHT file 
 * at the top-level directory of this distribution or in the current
 * directory.
 * 
 * This file is part of RMG. 
 * RMG is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * any later version.
 *
 * RMG is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
*/

/*
   Wrapper functions that simplify using various cuda, hip and sycl runtime calls.

   These functions use the same return semantics as their cuda or hip counterparts
   but for now RMG does not try to recover from an error condition so if the functions
   are not successful rmg::error is called and execution terminated so the caller
   does not need to check return codes.

*/

#include "GpuAlloc.h"
#include "rmg_error.h"
#include "transition.h"
#include "rmg_error.h"
#include "ErrorFuncs.h"


void MallocHostOrDevice(void **ptr, size_t size)
{
    if(ct.tddft_gpu)
    {
        gpuMalloc(ptr, size);
    }
    else
    {
        *ptr = malloc(size);
    }
}
void FreeHostOrDevice(void *ptr)
{
    if(ct.tddft_gpu)
    {
        gpuFree(ptr);
    }
    else
    {
        free(ptr);
    }
}

#if HIP_ENABLED

#ifdef __HIP_PLATFORM_NVCC__
#include <cuda.h>
#include <cuda_runtime_api.h>
#endif

hipError_t gpuMalloc(void **ptr, size_t size, const std::source_location loc)
{
    rmg::error(hipMalloc(ptr, size), loc);
    return hipSuccess;
}

hipError_t gpuMallocManaged(void **ptr, size_t size, const std::source_location loc)
{
    rmg::error(hipMallocManaged(ptr, size), loc);
    return hipSuccess;
}

hipError_t gpuMallocHost(void **ptr, size_t size, const std::source_location loc)
{
    if(ct.gpu_managed_memory)
    {
        rmg::error(hipMallocManaged(ptr, size), loc);
    }
    else
    {
        rmg::error(hipHostMalloc(ptr, size, hipHostMallocNumaUser), loc);
    }
    return hipSuccess;
}

hipError_t gpuFree(void *ptr, const std::source_location loc)
{
    rmg::error(hipFree(ptr), loc);
    return hipSuccess;
}

hipError_t gpuFreeHost(void *ptr, const std::source_location loc)
{
    if(ct.gpu_managed_memory)
    {
        rmg::error(hipFree(ptr), loc);
    }
    else
    {
        rmg::error(hipFreeHost(ptr), loc);
    }
    return hipSuccess;
}

hipError_t gpuMemcpy(void *dst, const void *src, size_t sizeBytes, hipMemcpyKind kind, const std::source_location loc)
{
    rmg::error(hipMemcpy(dst, src, sizeBytes, kind), loc);
    return hipSuccess;
}

hipError_t gpuMemcpyAsync (void *dst, const void *src, size_t sizeBytes, hipMemcpyKind kind, hipStream_t stream, const std::source_location loc)
{
    rmg::error(hipMemcpyAsync (dst, src, sizeBytes, kind, stream), loc);
    return hipSuccess;
}

hipError_t gpuMemPrefetchAsync ( const void* devPtr, size_t count, int  dstDevice, hipStream_t stream, const std::source_location loc)
{
    // Bit of a hack until HIP implements this
#ifdef __HIP_PLATFORM_NVCC__
    rmg::error(cudaMemPrefetchAsync (devPtr, count, dstDevice, stream), loc);
    return hipSuccess;
#else
    return hipSuccess;
#endif
}
hipError_t gpuStreamCreateWithFlags (hipStream_t *stream, unsigned int flags, const std::source_location loc)
{
    rmg::error(hipStreamCreateWithFlags (stream, flags), loc);
    return hipSuccess;
}
hipError_t gpuStreamDestroy (hipStream_t stream, const std::source_location loc)
{
    rmg::error(hipStreamDestroy (stream), loc);
    return hipSuccess;
}
hipError_t gpuMemcpy2D (void *dst, size_t dpitch, const void *src, size_t spitch, size_t width, size_t height, hipMemcpyKind kind, const std::source_location loc)
{
    rmg::error(hipMemcpy2D (dst, dpitch, src, spitch, width, height, kind), loc);
    return hipSuccess;
}
hipError_t gpuDeviceReset (const std::source_location loc)
{
    rmg::error(hipDeviceReset(), loc);
    return hipSuccess;
}
hipError_t gpuSetDevice (int deviceId, const std::source_location loc)
{
    rmg::error(hipSetDevice (deviceId), loc);
    return hipSuccess;
}
hipError_t gpuGetDevice (int *deviceId, const std::source_location loc)
{
    rmg::error(hipGetDevice (deviceId), loc);
    return hipSuccess;
}
hipError_t gpuSetDeviceFlags (unsigned flags, const std::source_location loc)
{
    rmg::error(hipSetDeviceFlags (flags), loc);
    return hipSuccess;
}
hipError_t gpuHostRegister(void *hostPtr, size_t sizeBytes, unsigned int flags, const std::source_location loc)
{
    rmg::error(hipHostRegister(hostPtr, sizeBytes, flags), loc);
    return hipSuccess;
}
hipError_t gpuHostUnregister(void *hostPtr, const std::source_location loc)	
{
    rmg::error(hipHostUnregister(hostPtr), loc);
    return hipSuccess;
}
hipError_t gpuGetDeviceCount(int *count, const std::source_location loc)
{
    rmg::error(hipGetDeviceCount(count), loc);
    return hipSuccess;
}
hipError_t gpuStreamSynchronize (hipStream_t stream, const std::source_location loc)
{
    rmg::error(hipStreamSynchronize (stream), loc);
    return hipSuccess;
}
#elif CUDA_ENABLED

cudaError_t gpuMalloc(void **ptr, size_t size, const std::source_location loc)
{
    rmg::error(cudaMalloc(ptr, size), loc);
    return cudaSuccess;
}

cudaError_t gpuMallocManaged(void **ptr, size_t size, const std::source_location loc)
{
    rmg::error(cudaMallocManaged(ptr, size), loc);
    return cudaSuccess;
}

cudaError_t gpuMallocHost(void **ptr, size_t size, const std::source_location loc)
{
    if(ct.gpu_managed_memory)
    {
        rmg::error(gpuMallocManaged(ptr, size), loc);
    }
    else
    {
        rmg::error(cudaMallocHost(ptr, size), loc);
    }
    return cudaSuccess;
}

cudaError_t gpuFree(void *ptr, const std::source_location loc)
{
    rmg::error(cudaFree(ptr), loc);
    return cudaSuccess;
}

cudaError_t gpuFreeHost(void *ptr, const std::source_location loc)
{
    if(ct.gpu_managed_memory)
    {
        rmg::error(cudaFree(ptr), loc);
    }
    else
    {
        rmg::error(cudaFreeHost(ptr), loc);
    }
    return cudaSuccess;
}

cudaError_t gpuMemcpy(void *dst, const void *src, size_t sizeBytes, cudaMemcpyKind kind, const std::source_location loc)
{
    rmg::error(cudaMemcpy(dst, src, sizeBytes, kind), loc);
    return cudaSuccess;
}

cudaError_t gpuMemcpyAsync (void *dst, const void *src, size_t sizeBytes, cudaMemcpyKind kind, cudaStream_t stream, const std::source_location loc)
{
    rmg::error(cudaMemcpyAsync (dst, src, sizeBytes, kind, stream), loc);
    return cudaSuccess;
}

cudaError_t gpuMemPrefetchAsync ( const void* devPtr, size_t count, int _dstDevice, cudaStream_t stream, const std::source_location loc)
{
#if CUDA_VERSION_MAJOR >= 13
    struct cudaMemLocation dstDevice = {
        .type = cudaMemLocationTypeDevice,
        .id = _dstDevice,
    };
    rmg::error(cudaMemPrefetchAsync (devPtr, count, dstDevice, 0, stream), loc);
#else
    rmg::error(cudaMemPrefetchAsync (devPtr, count, _dstDevice, stream), loc);
#endif
    return cudaSuccess;
}
cudaError_t gpuStreamCreateWithFlags (cudaStream_t *stream, unsigned int flags, const std::source_location loc)
{
    rmg::error(cudaStreamCreateWithFlags (stream, flags), loc);
    return cudaSuccess;
}
cudaError_t gpuStreamDestroy (cudaStream_t stream, const std::source_location loc)
{
    rmg::error(cudaStreamDestroy (stream), loc);
    return cudaSuccess;
}
cudaError_t gpuMemcpy2D (void *dst, size_t dpitch, const void *src, size_t spitch, size_t width, size_t height, cudaMemcpyKind kind, const std::source_location loc)
{
    rmg::error(cudaMemcpy2D (dst, dpitch, src, spitch, width, height, kind), loc);
    return cudaSuccess;
}
cudaError_t gpuDeviceReset (const std::source_location loc)
{
    rmg::error(cudaDeviceReset(), loc);
    return cudaSuccess;
}
cudaError_t gpuSetDevice (int deviceId, const std::source_location loc)
{
    rmg::error(cudaSetDevice (deviceId), loc);
    return cudaSuccess;
}
cudaError_t gpuGetDevice (int *deviceId, const std::source_location loc)
{
    rmg::error(cudaGetDevice (deviceId), loc);
    return cudaSuccess;
}
cudaError_t gpuSetDeviceFlags (unsigned flags, const std::source_location loc)
{
    rmg::error(cudaSetDeviceFlags (flags), loc);
    return cudaSuccess;
}
cudaError_t gpuHostRegister(void *hostPtr, size_t sizeBytes, unsigned int flags, const std::source_location loc)
{   
    rmg::error(cudaHostRegister(hostPtr, sizeBytes, flags), loc);
    return cudaSuccess;
}
cudaError_t gpuHostUnregister(void *hostPtr, const std::source_location loc)	
{
    rmg::error(cudaHostUnregister(hostPtr), loc);
    return cudaSuccess;
}
cudaError_t gpuGetDeviceCount(int *count, const std::source_location loc)
{
    rmg::error(cudaGetDeviceCount(count), loc);
    return cudaSuccess;
}
cudaError_t gpuStreamSynchronize (cudaStream_t stream, const std::source_location loc)
{
    rmg::error(cudaStreamSynchronize (stream), loc);
    return cudaSuccess;
}
#elif SYCL_ENABLED
#include <CL/sycl.hpp>
#include "GpuAlloc.h"
#include <omp.h>

//void FreeHostOrDevice(void *ptr)
//{
//    int device_id = omp_get_default_device();
//    omp_target_free(ptr, device_id);
//}
int gpuMalloc(void **ptr, size_t size)
{
    *ptr = cl::sycl::malloc_device(size, ct.sycl_Q);

    if(!ptr)
    {
        rmg::error("Error allocating host memory. Terminating.");
        return -1;
    }
    return 0;
}

//cudaError_t gpuMallocManaged(void **ptr, size_t size)
//{
//    cudaError_t cuerr = cudaMallocManaged(ptr, size);
//    if(cuerr != cudaSuccess)
//    rmg::error("Error allocating gpu memory. Terminating.");
//    return cuerr;
//}

int gpuMallocHost(void **ptr, size_t size)
{
    *ptr = cl::sycl::malloc_host(size, ct.sycl_Q);
    if(!ptr)
    {
        rmg::error("Error allocating host memory. Terminating.");
        return -1;
    }
    return 0;
}

int gpuFree(void *ptr)
{
    cl::sycl::free(ptr, ct.sycl_Q);
    return 0;
}

int gpuFreeHost(void *ptr)
{
    cl::sycl::free(ptr, ct.sycl_Q);
    return 0;
}
int gpuMemcpy(void *dst, const void *src, size_t sizeBytes, int kind)
{
    ct.sycl_Q.memcpy( dst, src, sizeBytes);
    ct.sycl_Q.wait();
    return 0;
}

#if 0
cudaError_t gpuMemcpyAsync (void *dst, const void *src, size_t sizeBytes, cudaMemcpyKind kind, cudaStream_t stream)
{
    return cudaMemcpyAsync (dst, src, sizeBytes, kind, stream);
}

cudaError_t gpuMemPrefetchAsync ( const void* devPtr, size_t count, int  dstDevice, cudaStream_t stream)
{
    return cudaMemPrefetchAsync (devPtr, count, dstDevice, stream);
}
cudaError_t gpuStreamCreateWithFlags (cudaStream_t *stream, unsigned int flags)
{
    return cudaStreamCreateWithFlags (stream, flags);
}
cudaError_t gpuStreamDestroy (cudaStream_t stream)
{
    return cudaStreamDestroy (stream);
}
cudaError_t gpuMemcpy2D (void *dst, size_t dpitch, const void *src, size_t spitch, size_t width, size_t height, cudaMemcpyKind kind)
{
    return cudaMemcpy2D (dst, dpitch, src, spitch, width, height, kind);
}
cudaError_t gpuDeviceReset (void)
{
    return cudaDeviceReset();
}
cudaError_t gpuSetDevice (int deviceId)
{
    return cudaSetDevice (deviceId);
}
cudaError_t gpuGetDevice (int *deviceId)
{
    return cudaGetDevice (deviceId);
}
cudaError_t gpuSetDeviceFlags (unsigned flags)
{
    return cudaSetDeviceFlags (flags);
}
cudaError_t gpuHostRegister(void *hostPtr, size_t sizeBytes, unsigned int flags)
{   
    return cudaHostRegister(hostPtr, sizeBytes, flags);
}
cudaError_t gpuHostUnregister(void *hostPtr)	
{
    return cudaHostUnregister(hostPtr);
}
cudaError_t gpuGetDeviceCount(int *count)
{
    return cudaGetDeviceCount(count);
}
#endif
#else
void gpuMalloc(void **ptr, size_t size)
{
    ptr = NULL;
}
void gpuFree(void *ptr)
{}
#endif

