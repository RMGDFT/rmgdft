/*
 *
 * Copyright 2022 The RMG Project Developers. See the COPYRIGHT file 
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


/*

    getGpuStream    Returns and if necessary initializes a per thread GPU stream handle for use by the caller.

*/

#include "BaseThread.h"
#include "Gpufuncs.h"

#if HIP_ENABLED
#include <hip/hip_runtime.h>
#include <hip/hip_ext.h>

hipStream_t getGpuStream(void)
{
    BaseThread *T = BaseThread::getBaseThread(0);
    int tid = T->get_thread_tid();
    if(tid < 0) tid = 0;
    if(!T->streams[tid])
    {
        hipStreamCreateWithFlags(&T->streams[tid], hipStreamNonBlocking);
    }
    return T->streams[tid];
}


#endif

#if HIP_ENABLED
#include <hip/hip_runtime.h>
#endif

#if CUDA_ENABLED
#include <cuda.h>
#include <cuda_runtime_api.h>

cudaStream_t getGpuStream(void)
{
    BaseThread *T = BaseThread::getBaseThread(0);
    int tid = T->get_thread_tid();
    if(tid < 0) tid = 0;
    if(!T->streams[tid])
    {
        cudaStreamCreateWithFlags(&T->streams[tid], cudaStreamNonBlocking);
    }
    return T->streams[tid];
}
#endif

