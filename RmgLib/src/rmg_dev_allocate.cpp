/*
 *
 * Copyright (c) 2026, Emil Briggs
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the <organization> nor the
 *       names of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.

 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * 
 */

#include "rmgtypedefs.h"
#include "rmg_dev_allocate.h"
#include "GpuAlloc.h"
#include <sys/mman.h>
#include <sys/sysinfo.h>
#include <sys/resource.h>
#include <complex>

#if CUDA_ENABLED
#include <cuda_runtime_api.h>
#include <cuda.h>
#endif

// This allocator is used for pinned host memory. The initial mmap call reserves
// a large contiguous address space but the pages are not initially mapped and pinned.
// When a malloc call is made we check to see if the new allocation would exceed the
// number of currently mapped and pinned pages. If not then we record the allocation
// and return a pointer to it. If it does exceed the number of currently mapped and
// pinned pages then we have to expand the space.
// 

namespace rmg
{
    size_t dev_allocate::allocated_pages;
    size_t dev_allocate::mapped_pages;
    std::byte *dev_allocate::baseptr;
    size_t gpu_pagesize;
    std::stack<std::pair<std::byte *, size_t>> dev_allocate::plist;
    int dev_allocate::device_id;
    size_t dev_allocate::deviceMem;

    dev_allocate::dev_allocate(int deviceId, size_t initial_size)
    {

#if CUDA_ENABLED
        device_id = deviceId;    // save for later
        rmg::error(cudaSetDevice(device_id));
        rmg::error(cuDeviceTotalMem( &deviceMem, device_id));

        CUmemAllocationProp props = {};
        props.type = CU_MEM_ALLOCATION_TYPE_PINNED;
        props.location.type = CU_MEM_LOCATION_TYPE_DEVICE;
        props.location.id = device_id;
        rmg::error(cuMemGetAllocationGranularity(&gpu_pagesize, &props, CU_MEM_ALLOC_GRANULARITY_MINIMUM));

        // Construction reserve a large contigous space
        if(!baseptr)
        {
            size_t rmem = initial_size / gpu_pagesize;
            rmem *= gpu_pagesize;
            size_t dmem = deviceMem / gpu_pagesize;
            dmem *= gpu_pagesize;
            CUmemGenericAllocationHandle_v1 handle1;
            rmg::error(cuMemAddressReserve((CUdeviceptr *)&baseptr, dmem, 0, 0, 0));
            rmg::error(cuMemCreate(&handle1, rmem, &props, 0));
            rmg::error(cuMemMap((CUdeviceptr)baseptr, rmem, 0, handle1, 0));

            CUmemAccessDesc accessDesc = {};
            accessDesc.location.type = CU_MEM_LOCATION_TYPE_DEVICE;
            accessDesc.location.id = device_id;
            accessDesc.flags = CU_MEM_ACCESS_FLAGS_PROT_READWRITE; 

            rmg::error(cuMemSetAccess((CUdeviceptr)baseptr, rmem, (CUmemAccessDesc *)&accessDesc, 1));
            rmg::error(cuMemRelease(handle1));
            mapped_pages = rmem / gpu_pagesize;
            plist.push(std::pair<std::byte *, size_t>(baseptr, 0));
        }
#elif HIP_ENABLED
        device_id = deviceId;    // save for later
        rmg::error(hipSetDevice(device_id));
        rmg::error(hipDeviceTotalMem( &deviceMem, device_id));

        hipMemAllocationProp props = {};
        props.type = hipMemAllocationTypePinned;
        props.location.type = hipMemLocationTypeDevice;
        props.location.id = device_id;
        rmg::error(hipMemGetAllocationGranularity(&gpu_pagesize, &props, hipMemAllocationGranularityRecommended));

        // Construction reserve a large contigous space
        if(!baseptr)
        {
            size_t rmem = initial_size / gpu_pagesize;
            rmem *= gpu_pagesize;
            hipMemGenericAllocationHandle_t handle1;
            rmg::error(hipMemAddressReserve((void **)&baseptr, deviceMem, 0, 0, 0));
            rmg::error(hipMemCreate(&handle1, rmem, &props, 0));
            rmg::error(hipMemMap(baseptr, rmem, 0, handle1, 0));
            hipMemAccessDesc accessDesc = {};
            accessDesc.location.type = hipMemLocationTypeDevice;
            accessDesc.location.id = device_id;
            accessDesc.flags = hipMemAccessFlagsProtReadWrite;
            rmg::error(hipMemSetAccess(baseptr, rmem, &accessDesc, 1));
            rmg::error(hipMemRelease(handle1));
            mapped_pages = rmem / gpu_pagesize;
            plist.push(std::pair<std::byte *, size_t>(baseptr, 0));
        }
#endif
    }

    dev_allocate::~dev_allocate(void)
    {
    }

    template<typename T>
    void dev_allocate::malloc(T **ptr, size_t size)
    {
        size_t nsize = sizeof(T) * size;
        size_t newpages = nsize / gpu_pagesize;
        if(nsize % gpu_pagesize) newpages++;
        size_t totalpages = allocated_pages + newpages;
        if(totalpages > mapped_pages)
        {
#if CUDA_ENABLED
            rmg::error(cudaSetDevice(device_id));
            CUmemAllocationProp props = {};
            props.type = CU_MEM_ALLOCATION_TYPE_PINNED;
            props.location.type = CU_MEM_LOCATION_TYPE_DEVICE;
            props.location.id = device_id;
            CUmemGenericAllocationHandle_v1 handle1;
            rmg::error(cuMemCreate(&handle1, newpages*gpu_pagesize, &props, 0));
            rmg::error(cuMemMap((CUdeviceptr)baseptr+mapped_pages*gpu_pagesize, newpages*gpu_pagesize, 0, handle1, 0));

            CUmemAccessDesc accessDesc = {};
            accessDesc.location.type = CU_MEM_LOCATION_TYPE_DEVICE;
            accessDesc.location.id = device_id;
            accessDesc.flags = CU_MEM_ACCESS_FLAGS_PROT_READWRITE; 
            mapped_pages += newpages;
            rmg::error(cuMemSetAccess((CUdeviceptr)baseptr, mapped_pages*gpu_pagesize, (CUmemAccessDesc *)&accessDesc, 1));
            rmg::error(cuMemRelease(handle1));
#elif HIP_ENABLED
            rmg::error(hipSetDevice(device_id));
            hipMemAllocationProp props = {};
            props.type = hipMemAllocationTypePinned;
            props.location.type = hipMemLocationTypeDevice;
            props.location.id = device_id;
            hipMemGenericAllocationHandle_t handle1;
            rmg::error(hipMemCreate(&handle1, newpages*gpu_pagesize, &props, 0));
            rmg::error(hipMemMap(baseptr+mapped_pages*gpu_pagesize, newpages*gpu_pagesize, 0, handle1, 0));

            hipMemAccessDesc accessDesc = {};
            accessDesc.location.type = hipMemLocationTypeDevice;
            accessDesc.location.id = device_id;
            accessDesc.flags = hipMemAccessFlagsProtReadWrite;
            mapped_pages += newpages;
            rmg::error(hipMemSetAccess(baseptr, mapped_pages*gpu_pagesize, &accessDesc, 1));
            rmg::error(hipMemRelease(handle1));
#endif
        }
//printf("AAAA  %lu  %lu  %lu  %lu\n",gpu_pagesize,totalpages, newpages, mapped_pages);
        
        std::pair<std::byte *, size_t> tpair = plist.top();
        T *nptr = (T *)tpair.first;
        *ptr = nptr; 
        std::byte *nptr1 = (std::byte *)nptr + newpages*gpu_pagesize;
        plist.push(std::pair<std::byte *, size_t>(nptr1, newpages));
        allocated_pages += newpages;
    }

    void dev_allocate::free(void *ptr)
    {
        if(plist.empty())
            rmg::error("Attempt to free non-existent allocation.");
        std::pair<std::byte *, size_t> tpair = plist.top();
        plist.pop();
        std::byte *nptr = tpair.first;
        std::byte *nptr1 = (std::byte *)ptr + tpair.second*gpu_pagesize;
        if(nptr != nptr1)
            rmg::error("Out of order allocation or free."); 
        allocated_pages -= tpair.second; 
    }

    template void dev_allocate::malloc(double **, size_t);
    template void dev_allocate::malloc(float **, size_t);
    template void dev_allocate::malloc(std::complex<double> **, size_t);
    template void dev_allocate::malloc(std::complex<float> **, size_t);
    template void dev_allocate::malloc(int **, size_t);
    template void dev_allocate::malloc(size_t **, size_t);
}

rmg::dev_allocate *rmg_device_pool;
