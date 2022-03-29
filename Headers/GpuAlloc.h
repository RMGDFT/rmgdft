#ifndef GPU_ALLOC_H
#define GPU_ALLOC_H 1

#include "stddef.h"

void *DGpuMallocDevice(size_t size, const char *fname, size_t line);
void InitGpuMalloc(size_t size);
void DGpuFreeDevice( void *ptr, const char *fname, size_t line );

void *DGpuMallocHost(size_t size, const char *fname, size_t line);
void InitGpuMallocHost(size_t size);
void DGpuFreeHost( void *ptr, const char *fname, size_t line);

void *DRmgMallocHost(size_t size, const char *fname, size_t line);
void DRmgFreeHost( void *ptr, const char *fname, size_t line);

#define  GpuMallocHost(x) DGpuMallocHost (x,__FILE__,__LINE__)
#define  GpuFreeHost(x) DGpuFreeHost (x,__FILE__,__LINE__)
#define  GpuMallocDevice(x) DGpuMallocDevice (x,__FILE__,__LINE__)
#define  GpuFreeDevice(x) DGpuFreeDevice (x,__FILE__,__LINE__)
#define  RmgMallocHost(x) DRmgMallocHost (x,__FILE__,__LINE__)
#define  RmgFreeHost(x) DRmgFreeHost (x,__FILE__,__LINE__)


void MallocHostOrDevice(void **ptr, size_t size);
void FreeHostOrDevice(void *ptr);
// The functions above manage blocks allocated from a pool. These are direct calls to
// the underlying functions.
#if HIP_ENABLED
#include <hip/hip_runtime.h>

#define gpublasDcopy hipblasDcopy
#define gpublasDdgmm hipblasDdgmm
#define gpublasDgeam hipblasDgeam
#define gpublasDscal hipblasDscal
#define gpublasDaxpy hipblasDaxpy
#define GPUBLAS_SIDE_LEFT HIPBLAS_SIDE_LEFT
#define GPUBLAS_SIDE_RIGHT HIPBLAS_SIDE_RIGHT
#define GPUBLAS_OP_N HIPBLAS_OP_N
#define gpuCpuDeviceId hipCpuDeviceId
#define gpuStream_t hipStream_t

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
hipError_t gpuHostRegister(void *hostPtr, size_t sizeBytes, unsigned int flags);
hipError_t gpuHostUnregister(void *hostPtr);
hipError_t gpuGetDeviceCount(int *count);

#elif CUDA_ENABLED
#include <cuda.h>
#include <cuda_runtime_api.h>

#define gpublasDcopy cublasDcopy
#define gpublasDdgmm cublasDdgmm
#define gpublasDgeam cublasDgeam
#define gpublasDscal cublasDscal
#define gpublasDaxpy cublasDaxpy
#define GPUBLAS_SIDE_LEFT CUBLAS_SIDE_LEFT
#define GPUBLAS_SIDE_RIGHT CUBLAS_SIDE_RIGHT
#define GPUBLAS_OP_N CUBLAS_OP_N
#define gpuCpuDeviceId cudaCpuDeviceId
#define gpuStream_t hipStream_t

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
cudaError_t gpuHostRegister(void *hostPtr, size_t sizeBytes, unsigned int flags);
cudaError_t gpuHostUnregister(void *hostPtr);
cudaError_t gpuGetDeviceCount(int *count); 

#define Cuda_error(err)                                                                     \
{                                                                                           \
    cudaError_t err_ = (err);                                                               \
    if (err_ != cudaSuccess) {                                                              \
        printf("CUDA error %d at %s: line %d\n", err_, __FILE__, __LINE__);                       \
        throw std::runtime_error("CUDA error");                                             \
    }                                                                                       \
}
#define Cusolver_status(custat)                                                             \
{                                                                                           \
    cusolverStatus_t custat_ = (custat);                                                    \
    if (custat_ != CUSOLVER_STATUS_SUCCESS) {                                               \
        printf("cusolver error %d at %s: line %d\n", custat_, __FILE__, __LINE__);                \
        throw std::runtime_error("cusolver error");                                         \
    }                                                                                       \
}



#else
void gpuMalloc(void **ptr, size_t size);
void gpuFree(void *ptr);
#endif


void MemcpyHostDevice (size_t a_size, void *A_host, void *A_device);
void MemcpyDeviceHost (size_t a_size, void *A_device, void *A_host);
template <typename T> T *MemoryPtrHostDevice(T *ptr_host, T *ptr_device);

#endif

