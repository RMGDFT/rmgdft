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

#include "make_conf.h"
#include "GpuAlloc.h"
#include "rmg_error.h"
#include "transition.h"
#include "ErrorFuncs.h"

#if GPU_ENABLED

#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cublas_v2.h>


// We also need a way to manage blocks of page locked memory for moving
// data into and out of the GPU since transfers to this type of memory
// are much faster. Certain data structures (e.g. the wavefunction arrays
// and non-local projectors) are allocated in this fashion separately
// but there are other methods/functions that need such memory on a
// temporary basis.


#define         GPU_ALIGNMENT   1024  
#define         MAX_HOSTGPU_BLOCKS  20


static void *host_gpubuffer;
static unsigned char *host_curptr;
static unsigned char *host_block_ptrs[MAX_HOSTGPU_BLOCKS];
static size_t host_block_sizes[MAX_HOSTGPU_BLOCKS];
static int host_allocated_blocks = 0;
static size_t host_cur_size = 0;
static size_t host_max_size;

void InitGpuMallocHost(size_t bufsize)
{

#if LINUX
    if (ct.require_huge_pages) return;
#endif

    cudaError_t custat;
    bufsize += GPU_ALIGNMENT * MAX_HOSTGPU_BLOCKS;
    custat = cudaMallocHost((void **)&host_gpubuffer , bufsize );
    RmgCudaError(__FILE__, __LINE__, custat, "Error: cudaHostMalloc failed in InitGpuHostMalloc\n");
    host_max_size = bufsize;
    host_curptr = (unsigned char*)host_gpubuffer;

}


void *DGpuMallocHost(size_t size, const char *fname, size_t line)
{

#if LINUX
    if(ct.require_huge_pages)
    {
        void *tptr = malloc(size);
        if(!tptr)
            rmg_error_handler (fname, line, "Error: Cannot allocate memory in GpuMallocHost.\n");

        cudaHostRegister( tptr, size, cudaHostRegisterPortable);
        return tptr;
    }
#endif
    size_t new_size, new_block;

    if(host_allocated_blocks == MAX_HOSTGPU_BLOCKS) {
        rmg_error_handler (fname, line, "Error: Too many blocks. Consider increasing MAX_HOSTGPU_BLOCKS.\n");
    }

    new_block = size / GPU_ALIGNMENT;
    new_block = (new_block + 1) * GPU_ALIGNMENT;
    new_size = host_cur_size + new_block; 
    if(new_size > host_max_size) {
        rmg_printf("GPU memory of %d bytes exceeds reserved size of %d.\n", new_size, host_max_size);
        rmg_error_handler (fname, line, "Error: Reservation too large. Consider increasing reserved GPU host memory.\n");
    }

    host_block_ptrs[host_allocated_blocks] = host_curptr; 
    host_block_sizes[host_allocated_blocks] = new_block; 
    host_curptr += new_block;
    host_cur_size += new_block;
    host_allocated_blocks++;
    return host_block_ptrs[host_allocated_blocks - 1];
}

void DGpuFreeHost(void *ptr, const char *fname, size_t line)
{
#if LINUX
    if(ct.require_huge_pages)
    {
        cudaHostUnregister( ptr );
        free(ptr);
        return;
    }
#endif

    if(host_allocated_blocks == 0) {
        rmg_printf("DEBUG: allocated_blocks = %d\n", host_allocated_blocks);
        rmg_error_handler (fname, line, "Error: Attempt to release non reserved block.\n");
    }

    if(ptr != host_block_ptrs[host_allocated_blocks-1]) {
        rmg_printf("DEBUG: ptr = %p    allocated_ptr = %p\n", ptr, host_block_ptrs[host_allocated_blocks-1]);
        rmg_error_handler (fname, line, "Error: Attempt to release non reserved block.\n");
    }

    host_allocated_blocks--;
    host_curptr -= host_block_sizes[host_allocated_blocks];
    host_cur_size -= host_block_sizes[host_allocated_blocks];
}

#endif
