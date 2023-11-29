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

#if HIP_ENABLED

#include <hip/hip_runtime.h>
void *DRmgMallocHost(size_t size, const char *fname, size_t line)
{
    void *ptr;
    hipError_t hipstat;
    if(ct.gpu_managed_memory)
    {
        hipstat = hipMallocManaged(&ptr, size);
        RmgGpuError(fname, line, hipstat, "Error: hipMallocManaaged failed.\n");
    }
    else
    {
        hipstat = hipHostMalloc( &ptr, size+16, hipHostMallocNumaUser);
        RmgGpuError(fname, line, hipstat, "Error: hipHostMalloc failed.\n");
    }
    return ptr;
}

void DRmgFreeHost(void *ptr, const char *fname, size_t line)
{
    if(ct.gpu_managed_memory)
    {
        hipFree(ptr);
    }
    else
    {
        hipFreeHost(ptr);
    }
}

#elif CUDA_ENABLED

#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cublas_v2.h>


void *DRmgMallocHost(size_t size, const char *fname, size_t line)
{
    void *ptr;
    cudaError_t custat;
    if(ct.gpu_managed_memory)
    {
        custat = cudaMallocManaged ( &ptr, size+16 );
        RmgGpuError(fname, line, custat, "Error: cudaMallocManaged failed.\n");
    }
    else
    {
        custat = cudaMallocHost ( &ptr, size+16 );
        RmgGpuError(fname, line, custat, "Error: cudaMallocHost failed.\n");
    }
    return ptr;
}

void DRmgFreeHost(void *ptr, const char *fname, size_t line)
{
    if(ct.gpu_managed_memory)
    {
        cudaFree(ptr);
    }
    else
    {
        cudaFreeHost(ptr);
    }
}

#elif SYCL_ENABLED

void *DRmgMallocHost(size_t size, const char *fname, size_t line)
{
    void *ptr;
    int retval = gpuMallocHost ( &ptr, size+16 );
    return ptr;
}

void DRmgFreeHost(void *ptr, const char *fname, size_t line)
{
    gpuFreeHost(ptr);
}

#else

void *DRmgMallocHost(size_t size, const char *fname, size_t line)
{
    void *ptr;
    if(NULL == (ptr = malloc(size + 16))) {
        printf("\n memory size required %lu ", size);
        throw RmgFatalException() << "memory cannot be allocated: " << fname << " at line " << line << ".\n";
    }
    return ptr;
}
void DRmgFreeHost(void *ptr, const char *fname, size_t line)
{
    free(ptr);
}

#endif

