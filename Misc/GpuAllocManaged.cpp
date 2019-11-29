/*
 *
 * Copyright 2014 The RMG Project Developers. See the COPYRIGHT file 
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

#include <complex>
#include <sys/mman.h>
#include "RmgException.h"

#if GPU_ENABLED

#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cublas_v2.h>


void *DGpuMallocManaged(size_t size, const char *fname, size_t line)
{
    void *ptr;
    cudaError_t custat;
    custat = cudaMallocManaged ( &ptr, size+16, cudaMemAttachGlobal );
    RmgCudaError(fname, line, custat, "Error: cudaMallocManaged failed.\n");
    return ptr;
}

void DGpuFreeManaged(void *ptr, const char *fname, size_t line)
{
    cudaFree(ptr);
}

#else

void *DGpuMallocManaged(size_t size, const char *fname, size_t line)
{
    void *ptr;
    if(NULL == (ptr = malloc(size + 16))) {
        printf("\n memory size required %lu ", size);
        throw RmgFatalException() << "memory cannot be allocated: " << fname << " at line " << line << ".\n";
    }
    return ptr;
}

void DGpuFreeManaged(void *ptr, const char *fname, size_t line)
{
    free(ptr);
}

#endif
