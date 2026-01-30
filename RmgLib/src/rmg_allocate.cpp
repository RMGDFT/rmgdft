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
#include "rmg_allocate.h"
#include "GpuAlloc.h"
#include <sys/mman.h>
#include <complex>

// This allocator is used for pinned host memory. The initial mmap call reserves
// a large contiguous address space but the pages are not initially mapped and pinned.
// When a malloc call is made we check to see if the new allocation would exceed the
// number of currently mapped and pinned pages. If not then we record the allocation
// and return a pointer to it. If it does exceed the number of currently mapped and
// pinned pages then we have to expand the space.
// 
namespace rmg
{
    size_t allocate::reserved_pages;
    size_t allocate::allocated_pages;
    size_t allocate::mapped_pages;
    std::byte *allocate::baseptr;
    size_t pagesize = 4096;
    std::stack<std::pair<std::byte *, size_t>> allocate::plist;

    allocate::allocate(void)
    {
        // Construction mmap a large contigous space
        if(!baseptr)
        {
            baseptr = (std::byte *)mmap(NULL, MAX_ALLOCATOR_SIZE, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
            if(baseptr == MAP_FAILED) rmg::error("Error reserving memory");
            reserved_pages = MAX_ALLOCATOR_SIZE / pagesize;
            plist.push(std::pair<std::byte *, size_t>(baseptr, 0));
        }
    }

    allocate::~allocate(void)
    {
    }

    template<typename T>
    void allocate::malloc(T **ptr, size_t size)
    {
        size_t nsize = sizeof(T) * size;
        size_t newpages = nsize / pagesize;
        if(nsize % pagesize) newpages++;
        size_t totalpages = allocated_pages + newpages;
        if(totalpages > mapped_pages)
        {
#if CUDA_ENABLED || HIP_ENABLED
            if(mapped_pages)
            {
                gpuHostUnregister(baseptr);
            }
            gpuHostRegister(baseptr, totalpages*pagesize, gpuHostRegisterPortable);
#endif
        }
        
        std::pair<std::byte *, size_t> tpair = plist.top();
        T *nptr = (T *)tpair.first;
        *ptr = nptr; 
        std::byte *nptr1 = (std::byte *)nptr + newpages*pagesize;
        plist.push(std::pair<std::byte *, size_t>(nptr1, newpages));
        allocated_pages += newpages;
    }

    void allocate::free(void *ptr)
    {
        if(plist.empty())
            rmg::error("Attempt to free non-existent allocation.");
        std::pair<std::byte *, size_t> tpair = plist.top();
        plist.pop();
        std::byte *nptr = tpair.first;
        std::byte *nptr1 = (std::byte *)ptr + tpair.second*pagesize;
        if(nptr != nptr1)
            rmg::error("Out of order allocation or free."); 
        allocated_pages -= tpair.second; 
    }
    template void allocate::malloc(double **, size_t);
    template void allocate::malloc(float **, size_t);
    template void allocate::malloc(std::complex<double> **, size_t);
    template void allocate::malloc(std::complex<float> **, size_t);
    template void allocate::malloc(int **, size_t);
    template void allocate::malloc(size_t **, size_t);
}

rmg::allocate hpool;
