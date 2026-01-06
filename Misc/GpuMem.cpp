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



#include "GpuAlloc.h"
#include "rmg_error.h"
#include "transition.h"
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
    hipError_t hiperr = hipMalloc(ptr, size);
    if(hiperr != hipSuccess)
        rmg::error("Error allocating gpu memory. Terminating.");
    return hiperr;
}

hipError_t gpuMallocManaged(void **ptr, size_t size, const std::source_location loc)
{
    hipError_t hiperr = hipMallocManaged(ptr, size);
    //hipError_t hiperr = hipHostMalloc(ptr, size, hipHostMallocNumaUser);

    if(hiperr != hipSuccess)
        rmg::error("Error allocating managed memory. Terminating.");
    return hiperr;
}

hipError_t gpuMallocHost(void **ptr, size_t size, const std::source_location loc)
{
    if(ct.gpu_managed_memory)
    {
        return gpuMallocManaged(ptr, size);
    }
    else
    {
        hipError_t hiperr = hipHostMalloc(ptr, size, hipHostMallocNumaUser);
        if(hiperr != hipSuccess)
        {
            std::cout << "memory size: " << size/1024 << " MB" << std::endl; 
            rmg::error("Error allocating pinned host memory . Terminating.");
        }
        return hiperr;
    }
}

hipError_t gpuFree(void *ptr, const std::source_location loc)
{
    return hipFree(ptr);
}

hipError_t gpuFreeHost(void *ptr, const std::source_location loc)
{
    if(ct.gpu_managed_memory)
    {
        return hipFree(ptr);
    }
    else
    {
        return hipFreeHost(ptr);
    }
}

hipError_t gpuMemcpy(void *dst, const void *src, size_t sizeBytes, hipMemcpyKind kind, const std::source_location loc)
{
    return hipMemcpy(dst, src, sizeBytes, kind);
}

hipError_t gpuMemcpyAsync (void *dst, const void *src, size_t sizeBytes, hipMemcpyKind kind, hipStream_t stream, const std::source_location loc)
{
    return hipMemcpyAsync (dst, src, sizeBytes, kind, stream);
}

hipError_t gpuMemPrefetchAsync ( const void* devPtr, size_t count, int  dstDevice, hipStream_t stream, const std::source_location loc)
{
    // Bit of a hack until HIP implements this
#ifdef __HIP_PLATFORM_NVCC__
    cudaMemPrefetchAsync (devPtr, count, dstDevice, stream);
    return hipSuccess;
#else
    return hipSuccess;
#endif
}
hipError_t gpuStreamCreateWithFlags (hipStream_t *stream, unsigned int flags, const std::source_location loc)
{
    return hipStreamCreateWithFlags (stream, flags);
}
hipError_t gpuStreamDestroy (hipStream_t stream, const std::source_location loc)
{
    return hipStreamDestroy (stream);
}
hipError_t gpuMemcpy2D (void *dst, size_t dpitch, const void *src, size_t spitch, size_t width, size_t height, hipMemcpyKind kind, const std::source_location loc)
{
    return hipMemcpy2D (dst, dpitch, src, spitch, width, height, kind);
}
hipError_t gpuDeviceReset (const std::source_location loc)
{
    return hipDeviceReset();
}
hipError_t gpuSetDevice (int deviceId, const std::source_location loc)
{
    return hipSetDevice (deviceId);
}
hipError_t gpuGetDevice (int *deviceId, const std::source_location loc)
{
    return hipGetDevice (deviceId);
}
hipError_t gpuSetDeviceFlags (unsigned flags, const std::source_location loc)
{
    return hipSetDeviceFlags (flags);
}
hipError_t gpuHostRegister(void *hostPtr, size_t sizeBytes, unsigned int flags, const std::source_location loc)
{
    return hipHostRegister(hostPtr, sizeBytes, flags);
}
hipError_t gpuHostUnregister(void *hostPtr, const std::source_location loc)	
{
    return hipHostUnregister(hostPtr);
}
hipError_t gpuGetDeviceCount(int *count, const std::source_location loc)
{
    return hipGetDeviceCount(count);
}
#elif CUDA_ENABLED

cudaError_t gpuMalloc(void **ptr, size_t size, const std::source_location loc)
{
    cudaError_t cuerr = cudaMalloc(ptr, size);
    if(cuerr != cudaSuccess)
    {
        std::cout << "size to be allocated " << size/1024.0/1024.0 <<" MB"<< std::endl;
        rmg::error("Error allocating gpu memory. Terminating.");
    }
    return cuerr;
}

cudaError_t gpuMallocManaged(void **ptr, size_t size, const std::source_location loc)
{
    cudaError_t cuerr = cudaMallocManaged(ptr, size);
    if(cuerr != cudaSuccess)
        rmg::error("Error allocating gpu memory. Terminating.");
    return cuerr;
}

cudaError_t gpuMallocHost(void **ptr, size_t size, const std::source_location loc)
{
    if(ct.gpu_managed_memory)
    {
        return gpuMallocManaged(ptr, size);
    }
    else
    {
        cudaError_t cuerr = cudaMallocHost(ptr, size);
        if(cuerr != cudaSuccess)
            rmg::error("Error allocating pinned host memory. Terminating.");
        return cuerr;
    }
}

cudaError_t gpuFree(void *ptr, const std::source_location loc)
{
    return cudaFree(ptr);
}

cudaError_t gpuFreeHost(void *ptr, const std::source_location loc)
{
    if(ct.gpu_managed_memory)
    {
        return cudaFree(ptr);
    }
    else
    {
        return cudaFreeHost(ptr);
    }
}

cudaError_t gpuMemcpy(void *dst, const void *src, size_t sizeBytes, cudaMemcpyKind kind, const std::source_location loc)
{
    return cudaMemcpy(dst, src, sizeBytes, kind);
}

cudaError_t gpuMemcpyAsync (void *dst, const void *src, size_t sizeBytes, cudaMemcpyKind kind, cudaStream_t stream, const std::source_location loc)
{
    return cudaMemcpyAsync (dst, src, sizeBytes, kind, stream);
}

cudaError_t gpuMemPrefetchAsync ( const void* devPtr, size_t count, int _dstDevice, cudaStream_t stream, const std::source_location loc)
{
#if CUDA_VERSION_MAJOR >= 13
    struct cudaMemLocation dstDevice = {
        .type = cudaMemLocationTypeDevice,
        .id = _dstDevice,
    };
    return cudaMemPrefetchAsync (devPtr, count, dstDevice, 0, stream);
#else
    return cudaMemPrefetchAsync (devPtr, count, _dstDevice, stream);
#endif

}
cudaError_t gpuStreamCreateWithFlags (cudaStream_t *stream, unsigned int flags, const std::source_location loc)
{
    return cudaStreamCreateWithFlags (stream, flags);
}
cudaError_t gpuStreamDestroy (cudaStream_t stream, const std::source_location loc)
{
    return cudaStreamDestroy (stream);
}
cudaError_t gpuMemcpy2D (void *dst, size_t dpitch, const void *src, size_t spitch, size_t width, size_t height, cudaMemcpyKind kind, const std::source_location loc)
{
    return cudaMemcpy2D (dst, dpitch, src, spitch, width, height, kind);
}
cudaError_t gpuDeviceReset (const std::source_location loc)
{
    return cudaDeviceReset();
}
cudaError_t gpuSetDevice (int deviceId, const std::source_location loc)
{
    return cudaSetDevice (deviceId);
}
cudaError_t gpuGetDevice (int *deviceId, const std::source_location loc)
{
    return cudaGetDevice (deviceId);
}
cudaError_t gpuSetDeviceFlags (unsigned flags, const std::source_location loc)
{
    return cudaSetDeviceFlags (flags);
}
cudaError_t gpuHostRegister(void *hostPtr, size_t sizeBytes, unsigned int flags, const std::source_location loc)
{   
    return cudaHostRegister(hostPtr, sizeBytes, flags);
}
cudaError_t gpuHostUnregister(void *hostPtr, const std::source_location loc)	
{
    return cudaHostUnregister(hostPtr);
}
cudaError_t gpuGetDeviceCount(int *count, const std::source_location loc)
{
    return cudaGetDeviceCount(count);
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

