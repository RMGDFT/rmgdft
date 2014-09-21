#ifndef GPU_ALLOC_H
#define GPU_ALLOC_H 1

#include "stddef.h"

void *GpuMalloc(size_t size);
void InitGpuMalloc(size_t size);
void GpuFree( void *ptr );

void *GpuMallocHost(size_t size);
void InitGpuMallocHost(size_t size);
void GpuFreeHost( void *ptr );


#endif
