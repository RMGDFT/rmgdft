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
#elif CUDA_ENABLED
#include <cuda.h>
#include <cuda_runtime_api.h>
cudaError_t gpuMalloc(void **ptr, size_t size);
cudaError_t gpuMallocManaged(void **ptr, size_t size);
cudaError_t gpuMallocHost(void **ptr, size_t size);
cudaError_t gpuFree(void *ptr);
cudaError_t gpuFreeHost(void *ptr);
#endif

#endif

