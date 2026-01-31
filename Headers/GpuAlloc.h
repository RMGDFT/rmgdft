#pragma once
#ifndef GPU_ALLOC_H
#define GPU_ALLOC_H 1

#include "stddef.h"
#include <source_location>

void *DRmgMallocHost(size_t size, const char *fname, size_t line);
void DRmgFreeHost( void *ptr, const char *fname, size_t line);

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
#define gpublasScopy hipblasScopy
#define gpublasSdgmm hipblasSdgmm
#define gpublasSgeam hipblasSgeam
#define gpublasSscal hipblasSscal
#define gpublasSaxpy hipblasSaxpy
#define gpublasCcopy hipblasCcopy
#define gpublasCdgmm hipblasCdgmm
#define gpublasCgeam hipblasCgeam
#define gpublasCscal hipblasCscal
#define gpublasCaxpy hipblasCaxpy
#define gpublasZcopy hipblasZcopy
#define gpublasZdgmm hipblasZdgmm
#define gpublasZgeam hipblasZgeam
#define gpublasZscal hipblasZscal
#define gpublasZaxpy hipblasZaxpy

#define GPUBLAS_SIDE_LEFT HIPBLAS_SIDE_LEFT
#define GPUBLAS_SIDE_RIGHT HIPBLAS_SIDE_RIGHT
#define GPUBLAS_OP_N HIPBLAS_OP_N
#define GPUBLAS_OP_T HIPBLAS_OP_T
#define gpuCpuDeviceId hipCpuDeviceId
#define gpuStream_t hipStream_t

void gpuMalloc(void **ptr, size_t size, std::source_location loc = std::source_location::current());
void gpuMallocManaged(void **ptr, size_t size, std::source_location loc = std::source_location::current());
void gpuMallocHost(void **ptr, size_t size, std::source_location loc = std::source_location::current());
void gpuFree(void *ptr, std::source_location loc = std::source_location::current());
void gpuFreeHost(void *ptr, std::source_location loc = std::source_location::current());
void gpuMemcpy(void *dst, const void *src, size_t sizeBytes, hipMemcpyKind kind, std::source_location loc = std::source_location::current());
void gpuMemcpyAsync (void *dst, const void *src, size_t sizeBytes, hipMemcpyKind kind, hipStream_t stream, std::source_location loc = std::source_location::current());
void gpuStreamSynchronize (hipStream_t stream, std::source_location loc = std::source_location::current());
void gpuMemPrefetchAsync ( const void* devPtr, size_t count, int  dstDevice, hipStream_t stream, std::source_location loc = std::source_location::current());
void gpuStreamCreateWithFlags (hipStream_t *stream, unsigned int flags, std::source_location loc = std::source_location::current());
void gpuStreamDestroy (hipStream_t stream, std::source_location loc = std::source_location::current());
void gpuMemcpy2D (void *dst, size_t dpitch, const void *src, size_t spitch, size_t width, size_t height, hipMemcpyKind kind, std::source_location loc = std::source_location::current());
void gpuDeviceReset (std::source_location loc = std::source_location::current());
void gpuSetDevice (int deviceId, std::source_location loc = std::source_location::current());
void gpuGetDevice (int *deviceId, std::source_location loc = std::source_location::current());
void gpuSetDeviceFlags (unsigned flags, std::source_location loc = std::source_location::current());
void gpuHostRegister(void *hostPtr, size_t sizeBytes, unsigned int flags, std::source_location loc = std::source_location::current());
void gpuHostUnregister(void *hostPtr, std::source_location loc = std::source_location::current());
void gpuGetDeviceCount(int *count, std::source_location loc = std::source_location::current());

#elif CUDA_ENABLED
#include <cuda.h>
#include <cuda_runtime_api.h>

#define gpublasDcopy cublasDcopy
#define gpublasDdgmm cublasDdgmm
#define gpublasDgeam cublasDgeam
#define gpublasDscal cublasDscal
#define gpublasDaxpy cublasDaxpy
#define gpublasScopy cublasScopy
#define gpublasSdgmm cublasSdgmm
#define gpublasSgeam cublasSgeam
#define gpublasSscal cublasSscal
#define gpublasSaxpy cublasSaxpy
#define gpublasCcopy cublasCcopy
#define gpublasCdgmm cublasCdgmm
#define gpublasCgeam cublasCgeam
#define gpublasCscal cublasCscal
#define gpublasCaxpy cublasCaxpy
#define gpublasZcopy cublasZcopy
#define gpublasZdgmm cublasZdgmm
#define gpublasZgeam cublasZgeam
#define gpublasZscal cublasZscal
#define gpublasZaxpy cublasZaxpy

#define GPUBLAS_SIDE_LEFT CUBLAS_SIDE_LEFT
#define GPUBLAS_SIDE_RIGHT CUBLAS_SIDE_RIGHT
#define GPUBLAS_OP_N CUBLAS_OP_N
#define GPUBLAS_OP_T CUBLAS_OP_T
#define gpuCpuDeviceId cudaCpuDeviceId
#define gpuStream_t cudaStream_t

void gpuMalloc(void **ptr, size_t size, std::source_location loc = std::source_location::current());
void gpuMallocManaged(void **ptr, size_t size, std::source_location loc = std::source_location::current());
void gpuMallocHost(void **ptr, size_t size, std::source_location loc = std::source_location::current());
void gpuFree(void *ptr, std::source_location loc = std::source_location::current());
void gpuFreeHost(void *ptr, std::source_location loc = std::source_location::current());
void gpuMemcpy(void *dst, const void *src, size_t sizeBytes, cudaMemcpyKind kind, std::source_location loc = std::source_location::current());
void gpuMemcpyAsync (void *dst, const void *src, size_t sizeBytes, cudaMemcpyKind kind, cudaStream_t stream, std::source_location loc = std::source_location::current());
void gpuStreamSynchronize (cudaStream_t stream, std::source_location loc = std::source_location::current());
void gpuMemPrefetchAsync ( const void* devPtr, size_t count, int  dstDevice, cudaStream_t stream, std::source_location loc = std::source_location::current());
void gpuStreamCreateWithFlags (cudaStream_t *stream, unsigned int flags, std::source_location loc = std::source_location::current());
void gpuStreamDestroy (cudaStream_t stream, std::source_location loc = std::source_location::current());
void gpuMemcpy2D (void *dst, size_t dpitch, const void *src, size_t spitch, size_t width, size_t height, cudaMemcpyKind kind, std::source_location loc = std::source_location::current());
void gpuDeviceReset (std::source_location loc = std::source_location::current());
void gpuSetDevice (int deviceId, std::source_location loc = std::source_location::current());
void gpuGetDevice (int *deviceId, std::source_location loc = std::source_location::current());
void gpuSetDeviceFlags (unsigned flags, std::source_location loc = std::source_location::current());
void gpuHostRegister(void *hostPtr, size_t sizeBytes, unsigned int flags, std::source_location loc = std::source_location::current());
void gpuHostUnregister(void *hostPtr, std::source_location loc = std::source_location::current());
void gpuGetDeviceCount(int *count, std::source_location loc = std::source_location::current()); 

#elif SYCL_ENABLED

int gpuMalloc(void **ptr, size_t size);
int gpuMallocManaged(void **ptr, size_t size);
int gpuMallocHost(void **ptr, size_t size);
int gpuFree(void *ptr);
int gpuFreeHost(void *ptr);
int gpuMemcpy(void *dst, const void *src, size_t sizeBytes, int kind);
//int gpuMemcpyAsync (void *dst, const void *src, size_t sizeBytes, cudaMemcpyKind kind, cudaStream_t stream);
//int gpuStreamSynchronize (cudaStream_t stream);
//int gpuMemPrefetchAsync ( const void* devPtr, size_t count, int  dstDevice, cudaStream_t stream);
//int gpuStreamCreateWithFlags (cudaStream_t *stream, unsigned int flags);
//int gpuStreamDestroy (cudaStream_t stream);
//int gpuMemcpy2D (void *dst, size_t dpitch, const void *src, size_t spitch, size_t width, size_t height, cudaMemcpyKind kind);
int gpuDeviceReset (void);
int gpuSetDevice (int deviceId);
int gpuGetDevice (int *deviceId);
int gpuSetDeviceFlags (unsigned flags);
int gpuHostRegister(void *hostPtr, size_t sizeBytes, unsigned int flags);
int gpuHostUnregister(void *hostPtr);
int gpuGetDeviceCount(int *count); 

#else
void gpuMalloc(void **ptr, size_t size);
void gpuFree(void *ptr);
#endif


void RmgMemcpy (void *A_dest, void *A_src, size_t a_size);
void MemcpyHostDevice (size_t a_size, void *A_host, void *A_device);
void MemcpyDeviceHost (size_t a_size, void *A_device, void *A_host);
template <typename T> T *MemoryPtrHostDevice(T *ptr_host, T *ptr_device);

#endif

