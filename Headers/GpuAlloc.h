#ifndef GPU_ALLOC_H
#define GPU_ALLOC_H 1

#include "stddef.h"

void *DGpuMallocDevice(size_t size, const char *fname, size_t line);
void InitGpuMalloc(size_t size);
void DGpuFreeDevice( void *ptr, const char *fname, size_t line );

void *DGpuMallocHost(size_t size, const char *fname, size_t line);
void InitGpuMallocHost(size_t size);
void DGpuFreeHost( void *ptr, const char *fname, size_t line);

#define  GpuMallocHost(x) DGpuMallocHost (x,__FILE__,__LINE__)
#define  GpuFreeHost(x) DGpuFreeHost (x,__FILE__,__LINE__)
#define  GpuMallocDevice(x) DGpuMallocDevice (x,__FILE__,__LINE__)
#define  GpuFreeDevice(x) DGpuFreeDevice (x,__FILE__,__LINE__)

#endif
