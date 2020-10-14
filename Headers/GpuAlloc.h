#ifndef GPU_ALLOC_H
#define GPU_ALLOC_H 1

#include "stddef.h"

void *DGpuMallocDevice(size_t size, const char *fname, size_t line);
void InitGpuMalloc(size_t size);
void DGpuFreeDevice( void *ptr, const char *fname, size_t line );

void *DGpuMallocHost(size_t size, const char *fname, size_t line);
void InitGpuMallocHost(size_t size);
void DGpuFreeHost( void *ptr, const char *fname, size_t line);

void *DGpuMallocManaged(size_t size, const char *fname, size_t line);
void DGpuFreeManaged( void *ptr, const char *fname, size_t line);

#define  GpuMallocHost(x) DGpuMallocHost (x,__FILE__,__LINE__)
#define  GpuFreeHost(x) DGpuFreeHost (x,__FILE__,__LINE__)
#define  GpuMallocDevice(x) DGpuMallocDevice (x,__FILE__,__LINE__)
#define  GpuFreeDevice(x) DGpuFreeDevice (x,__FILE__,__LINE__)
#define  GpuMallocManaged(x) DGpuMallocManaged (x,__FILE__,__LINE__)
#define  GpuFreeManaged(x) DGpuFreeManaged (x,__FILE__,__LINE__)


// The functions above manage blocks allocated from a pool. These are direct calls to
// the underlying functions.
#if HIP_ENABLED
#include <hip/hip_runtime.h>
hipError_t gpuMalloc(void **ptr, size_t size);
hipError_t gpuMallocManaged(void **ptr, size_t size);
hipError_t gpuMallocHost(void **ptr, size_t size);
hipError_t gpuFree(void *ptr);
hipError_t gpuFreeHost(void *ptr);
hipError_t gpuMemcpy(void *dst, const void *src, size_t sizeBytes, hipMemcpyKind kind);
hipError_t gpuMemcpyAsync (void *dst, const void *src, size_t sizeBytes, hipMemcpyKind kind, hipStream_t stream);
hipError_t gpuStreamSynchronize (hipStream_t stream);
hipError_t gpuMemPrefetchAsync ( const void* devPtr, size_t count, int  dstDevice, hipStream_t stream);
hipError_t gpuStreamCreateWithFlags (hipStream_t *stream, unsigned int flags);
hipError_t gpuStreamDestroy (hipStream_t stream);
hipError_t gpuMemcpy2D (void *dst, size_t dpitch, const void *src, size_t spitch, size_t width, size_t height, hipMemcpyKind kind);
hipError_t gpuDeviceReset (void);
hipError_t gpuSetDevice (int deviceId);
hipError_t gpuGetDevice (int *deviceId);
hipError_t gpuSetDeviceFlags (unsigned flags);

#elif CUDA_ENABLED
#include <cuda.h>
#include <cuda_runtime_api.h>
cudaError_t gpuMalloc(void **ptr, size_t size);
cudaError_t gpuMallocManaged(void **ptr, size_t size);
cudaError_t gpuMallocHost(void **ptr, size_t size);
cudaError_t gpuFree(void *ptr);
cudaError_t gpuFreeHost(void *ptr);
cudaError_t gpuMemcpy(void *dst, const void *src, size_t sizeBytes, cudaMemcpyKind kind);
cudaError_t gpuMemcpyAsync (void *dst, const void *src, size_t sizeBytes, cudaMemcpyKind kind, cudaStream_t stream);
cudaError_t gpuStreamSynchronize (cudaStream_t stream);
cudaError_t gpuMemPrefetchAsync ( const void* devPtr, size_t count, int  dstDevice, cudaStream_t stream);
cudaError_t gpuStreamCreateWithFlags (cudaStream_t *stream, unsigned int flags);
cudaError_t gpuStreamDestroy (cudaStream_t stream);
cudaError_t gpuMemcpy2D (void *dst, size_t dpitch, const void *src, size_t spitch, size_t width, size_t height, cudaMemcpyKind kind);
cudaError_t gpuDeviceReset (void);
cudaError_t gpuSetDevice (int deviceId);
cudaError_t gpuGetDevice (int *deviceId);
cudaError_t gpuSetDeviceFlags (unsigned flags);

#endif

#endif

