#ifndef GPU_ALLOC_H
#define GPU_ALLOC_H 1

#include "stddef.h"

void *GpuMalloc(size_t size);
void InitGpuMalloc(size_t size);
void GpuFree( void *ptr );



#endif
