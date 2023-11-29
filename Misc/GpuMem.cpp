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



#if HIP_ENABLED

#ifdef __HIP_PLATFORM_NVCC__
    #include <cuda.h>
    #include <cuda_runtime_api.h>
#endif
void MallocHostOrDevice(void **ptr, size_t size)
{
    gpuMalloc(ptr, size);
}
void FreeHostOrDevice(void *ptr)
{
    gpuFree(ptr);
}

hipError_t gpuMalloc(void **ptr, size_t size)
{
    hipError_t hiperr = hipMalloc(ptr, size);
    if(hiperr != hipSuccess)
    rmg_error_handler(__FILE__, __LINE__, "Error allocating gpu memory. Terminating.");
    return hiperr;
}

hipError_t gpuMallocManaged(void **ptr, size_t size)
{
    hipError_t hiperr = hipMallocManaged(ptr, size);
    //hipError_t hiperr = hipHostMalloc(ptr, size, hipHostMallocNumaUser);

    if(hiperr != hipSuccess)
    rmg_error_handler(__FILE__, __LINE__, "Error allocating managed memory. Terminating.");
    return hiperr;
}

hipError_t gpuMallocHost(void **ptr, size_t size)
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
            rmg_error_handler(__FILE__, __LINE__, "Error allocating pinned host memory . Terminating.");
        }
        return hiperr;
    }
}

hipError_t gpuFree(void *ptr)
{
    return hipFree(ptr);
}

hipError_t gpuFreeHost(void *ptr)
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

hipError_t gpuMemcpy(void *dst, const void *src, size_t sizeBytes, hipMemcpyKind kind)
{
    return hipMemcpy(dst, src, sizeBytes, kind);
}

hipError_t gpuMemcpyAsync (void *dst, const void *src, size_t sizeBytes, hipMemcpyKind kind, hipStream_t stream)
{
    return hipMemcpyAsync (dst, src, sizeBytes, kind, stream);
}

hipError_t gpuMemPrefetchAsync ( const void* devPtr, size_t count, int  dstDevice, hipStream_t stream)
{
// Bit of a hack until HIP implements this
#ifdef __HIP_PLATFORM_NVCC__
    cudaMemPrefetchAsync (devPtr, count, dstDevice, stream);
    return hipSuccess;
#else
    return hipSuccess;
#endif
}
hipError_t gpuStreamCreateWithFlags (hipStream_t *stream, unsigned int flags)
{
    return hipStreamCreateWithFlags (stream, flags);
}
hipError_t gpuStreamDestroy (hipStream_t stream)
{
    return hipStreamDestroy (stream);
}
hipError_t gpuMemcpy2D (void *dst, size_t dpitch, const void *src, size_t spitch, size_t width, size_t height, hipMemcpyKind kind)
{
    return hipMemcpy2D (dst, dpitch, src, spitch, width, height, kind);
}
hipError_t gpuDeviceReset (void)
{
    return hipDeviceReset();
}
hipError_t gpuSetDevice (int deviceId)
{
    return hipSetDevice (deviceId);
}
hipError_t gpuGetDevice (int *deviceId)
{
    return hipGetDevice (deviceId);
}
hipError_t gpuSetDeviceFlags (unsigned flags)
{
    return hipSetDeviceFlags (flags);
}
hipError_t gpuHostRegister(void *hostPtr, size_t sizeBytes, unsigned int flags)
{
    return hipHostRegister(hostPtr, sizeBytes, flags);
}
hipError_t gpuHostUnregister(void *hostPtr)	
{
    return hipHostUnregister(hostPtr);
}
hipError_t gpuGetDeviceCount(int *count)
{
    return hipGetDeviceCount(count);
}
#elif CUDA_ENABLED

void MallocHostOrDevice(void **ptr, size_t size)
{
   gpuMalloc(ptr, size);
}
void FreeHostOrDevice(void *ptr)
{
    gpuFree(ptr);
}
cudaError_t gpuMalloc(void **ptr, size_t size)
{
    cudaError_t cuerr = cudaMalloc(ptr, size);
    if(cuerr != cudaSuccess)
    {
        std::cout << "size to be allocated " << size/1024.0/1024.0 <<" MB"<< std::endl;
        rmg_error_handler(__FILE__, __LINE__, "Error allocating gpu memory. Terminating.");
    }
    return cuerr;
}

cudaError_t gpuMallocManaged(void **ptr, size_t size)
{
    cudaError_t cuerr = cudaMallocManaged(ptr, size);
    if(cuerr != cudaSuccess)
    rmg_error_handler(__FILE__, __LINE__, "Error allocating gpu memory. Terminating.");
    return cuerr;
}

cudaError_t gpuMallocHost(void **ptr, size_t size)
{
    if(ct.gpu_managed_memory)
    {
        return gpuMallocManaged(ptr, size);
    }
    else
    {
        cudaError_t cuerr = cudaMallocHost(ptr, size);
        if(cuerr != cudaSuccess)
        rmg_error_handler(__FILE__, __LINE__, "Error allocating pinned host memory. Terminating.");
        return cuerr;
    }
}

cudaError_t gpuFree(void *ptr)
{
    return cudaFree(ptr);
}

cudaError_t gpuFreeHost(void *ptr)
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

cudaError_t gpuMemcpy(void *dst, const void *src, size_t sizeBytes, cudaMemcpyKind kind)
{
    return cudaMemcpy(dst, src, sizeBytes, kind);
}

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
#elif SYCL_ENABLED
#include <omp.h>

void MallocHostOrDevice(void **ptr, size_t size)
{
   gpuMalloc(ptr, size);
}
void FreeHostOrDevice(void *ptr)
{
    int device_id = omp_get_default_device();
    omp_target_free(ptr, device_id);
}
int gpuMalloc(void **ptr, size_t size)
{
    int device_id = omp_get_default_device();
    *ptr = (void *) omp_target_alloc_device(size, device_id);

    if(!ptr)
    {
        rmg_error_handler(__FILE__, __LINE__, "Error allocating host memory. Terminating.");
        return -1;
    }
    return 0;
}

//cudaError_t gpuMallocManaged(void **ptr, size_t size)
//{
//    cudaError_t cuerr = cudaMallocManaged(ptr, size);
//    if(cuerr != cudaSuccess)
//    rmg_error_handler(__FILE__, __LINE__, "Error allocating gpu memory. Terminating.");
//    return cuerr;
//}

int gpuMallocHost(void **ptr, size_t size)
{
    int device_id = omp_get_initial_device();
    *ptr = (void *) omp_target_alloc_shared(size, device_id);

    if(!ptr)
    {
        rmg_error_handler(__FILE__, __LINE__, "Error allocating host memory. Terminating.");
        return -1;
    }
    return 0;
}

int gpuFree(void *ptr)
{
    int device_id = omp_get_default_device();
    omp_target_free(ptr, device_id);
    return 0;
}

int gpuFreeHost(void *ptr)
{
    int device_id = omp_get_initial_device();
    omp_target_free(ptr, device_id);
    return 0;
}
#if 0
int gpuMemcpy(void *dst, const void *src, size_t sizeBytes, cudaMemcpyKind kind)
{
    return cudaMemcpy(dst, src, sizeBytes, kind);
}

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
void MallocHostOrDevice(void **ptr, size_t size)
{
    *ptr = malloc(size);
}
void FreeHostOrDevice(void *ptr)
{
    free(ptr);
}
void gpuMalloc(void **ptr, size_t size)
{
   ptr = NULL;
}
void gpuFree(void *ptr)
{}
#endif

