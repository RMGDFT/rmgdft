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
    custat = cudaMalloc((void **)&gpubuffer , bufsize );
    RmgCudaError(__FILE__, __LINE__, custat, "Error: cudaMalloc failed in InitGpuMalloc\n");
    max_size = bufsize;
    curptr = (unsigned char*)gpubuffer;

}


void *GpuMalloc(size_t size)
{
    size_t new_size, new_block;

    if(allocated_blocks == MAX_GPU_BLOCKS) {
        rmg_error_handler (__FILE__, __LINE__, "Error: Too many blocks. Consider increasing MAX_GPU_BLOCKS.\n");
    }

    new_block = size / GPU_ALIGNMENT;
    new_block = (new_block + 1) * GPU_ALIGNMENT;
    new_size = cur_size + new_block; 
    if(new_size > max_size) {
        rmg_printf("GPU memory of %d bytes exceeds reserved size.\n", size);
        rmg_error_handler (__FILE__, __LINE__, "Error: Reservation too large. Consider increasing reserved GPU memory.\n");
    }

    block_ptrs[allocated_blocks] = curptr; 
    block_sizes[allocated_blocks] = new_block; 
    curptr += new_block;
    cur_size += new_block;
    allocated_blocks++;
    return block_ptrs[allocated_blocks - 1];
}

void GpuFree(void *ptr)
{
    if(allocated_blocks == 0) {
        rmg_printf("DEBUG: allocated_blocks = %d\n", allocated_blocks);
        rmg_error_handler (__FILE__, __LINE__, "Error: Attempt to release non reserved block.\n");
    }

    if(ptr != block_ptrs[allocated_blocks-1]) {
        rmg_printf("DEBUG: ptr = %p    allocated_ptr = %p\n", ptr, block_ptrs[allocated_blocks-1]);
        rmg_error_handler (__FILE__, __LINE__, "Error: Attempt to release non reserved block.\n");
    }

    allocated_blocks--;
    curptr -= block_sizes[allocated_blocks];
    cur_size -= block_sizes[allocated_blocks];
}

#endif
