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


// Since cudamalloc is very slow we grab a large buffer of GPU ram
// at initialization time and use GpuMalloc to return pointers into
// it. This is not a sophisticated allocater. One must know how much
// gpu memory the program needs at compile time and set the initial 
// block size large enough. Also blocks must be allocated and freed
// in FIFO fashion. Considering the way that RMG uses GPU memory
// at the present time with a few large blocks in subdiag this is 
// not a serious restriction.

#define         GPU_ALIGNMENT   1024  
#define         MAX_GPU_BLOCKS      12
static void *gpubuffer;
static unsigned char *curptr;
static unsigned char *block_ptrs[MAX_GPU_BLOCKS];
static size_t block_sizes[MAX_GPU_BLOCKS];
static int allocated_blocks = 0;
static size_t cur_size = 0;
static size_t max_size;


void InitGpuMalloc(size_t bufsize)
{

    cudaError_t custat;

    // Add alignment factors
    bufsize += GPU_ALIGNMENT * MAX_GPU_BLOCKS;
    custat = cudaMalloc((void **)&gpubuffer , bufsize  );
    RmgCudaError(__FILE__, __LINE__, custat, "Error: cudaMalloc failed in InitGpuMalloc\n");
    max_size = bufsize;
    curptr = (unsigned char*)gpubuffer;

}


void *DGpuMallocDevice(size_t size, const char *fname, size_t line)
{
    size_t new_size, new_block;

    if(allocated_blocks == MAX_GPU_BLOCKS) {
        rmg_error_handler (fname, line, "Error: Too many blocks. Consider increasing MAX_GPU_BLOCKS.\n");
    }

    new_block = size / GPU_ALIGNMENT;
    new_block = (new_block + 1) * GPU_ALIGNMENT;
    new_size = cur_size + new_block; 
    if(new_size > max_size) {
        rmg_printf("GPU memory of %d bytes exceeds reserved size of %d.\n", new_size, max_size);
        rmg_error_handler (fname, line, "Error: Reservation too large. Consider increasing reserved GPU memory.\n");
    }

    block_ptrs[allocated_blocks] = curptr; 
    block_sizes[allocated_blocks] = new_block; 
    curptr += new_block;
    cur_size += new_block;
    allocated_blocks++;
    return block_ptrs[allocated_blocks - 1];
}

void DGpuFreeDevice(void *ptr, const char *fname, size_t line)
{
    if(allocated_blocks == 0) {
        rmg_printf("DEBUG: allocated_blocks = %d\n", allocated_blocks);
        rmg_error_handler (fname, line, "Error: Attempt to release non reserved block.\n");
    }

    if(ptr != block_ptrs[allocated_blocks-1]) {
        rmg_printf("DEBUG: ptr = %p    allocated_ptr = %p\n", ptr, block_ptrs[allocated_blocks-1]);
        rmg_error_handler (fname, line, "Error: Attempt to release non reserved block.\n");
    }

    allocated_blocks--;
    curptr -= block_sizes[allocated_blocks];
    cur_size -= block_sizes[allocated_blocks];
}

#endif
