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
#include <hip/hip_runtime.h>

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
    if(hiperr != hipSuccess)
    rmg_error_handler(__FILE__, __LINE__, "Error allocating managed memory. Terminating.");
    return hiperr;
}

hipError_t gpuMallocHost(void **ptr, size_t size)
{
    hipError_t hiperr = hipMallocHost(ptr, size);
    if(hiperr != hipSuccess)
    rmg_error_handler(__FILE__, __LINE__, "Error allocating pinned host memory. Terminating.");
    return hiperr;
}

void gpuFree(void *ptr)
{
    hipFree(ptr);
}

#elif CUDA_ENABLED
#include <cuda.h>
#include <cuda_runtime_api.h>

cudaError_t gpuMalloc(void **ptr, size_t size)
{
    cudaError_t cuerr = cudaMalloc(ptr, size);
    if(cuerr != cudaSuccess)
    rmg_error_handler(__FILE__, __LINE__, "Error allocating gpu memory. Terminating.");
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
    cudaError_t cuerr = cudaMallocHost(ptr, size);
    if(cuerr != cudaSuccess)
    rmg_error_handler(__FILE__, __LINE__, "Error allocating pinned host memory. Terminating.");
    return cuerr;
}

void gpuFree(void *ptr)
{
    gpuFree(ptr);
}

#endif

