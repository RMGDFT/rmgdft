/*
 *
 * Copyright 2018 The RMG Project Developers. See the COPYRIGHT file 
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


#if HIP_ENABLED

#include <hip/hip_runtime.h>
#include <hip/hip_runtime.h>
#include <hip/hip_runtime_api.h>
#include <hip/hip_ext.h>
#include <stdio.h>
#include <unistd.h>
#include "Gpufuncs.h"
#include <iostream>
#include <vector>

//extern "C" void *BeginRmgTimer(const char *what);
//extern "C" void EndRmgTimer(void *ptr);

void GetPrimeFactors(std::vector<int>& factors, int val, int stop);

#define IMAGES 4


// Version with shared memory slice
template <typename T>
__global__ void app8_del2_kernel(const T * __restrict__ a, 
                                                T *b, 
                                                const int dimx,
                                                const int dimy,
                                                const int dimz,
                                                const fdparms_o8<T> c)

{

    extern __shared__ __align__(sizeof(T)) unsigned char sbuf[];
    T *slice = reinterpret_cast<T *>(sbuf);

    // iy and iz map to the x and y coordinates of the thread
    // within a block
    const int iy = blockIdx.x * (blockDim.x-2*IMAGES) + threadIdx.x;
    const int iz = blockIdx.y * (blockDim.y-2*IMAGES) + threadIdx.y;

    // thread y-index into shared memory tile
    const int ty = threadIdx.x;
    const int tz = threadIdx.y;

    int in_idx = iy*(dimz + 2*IMAGES) + iz;
    int in_stride = (dimy + 2*IMAGES) * (dimz + 2*IMAGES);

    int out_idx = (iy-IMAGES)*dimz + (iz-IMAGES);

    // All threads load array to shared memory. This bool is
    // used to identify which threads compute and write results
    bool execute = ((threadIdx.x >= IMAGES) &&
           (threadIdx.y >= IMAGES) &&
           (threadIdx.x < (blockDim.x-IMAGES)) &&
           (threadIdx.y < (blockDim.y-IMAGES)));


    T acc_b4 = 0.0;
    T acc_b3 = a[in_idx];in_idx += in_stride;
    T acc_b2 = a[in_idx];in_idx += in_stride;
    T acc_b1 = a[in_idx];in_idx += in_stride;
    T acc    = a[in_idx];in_idx += in_stride;
    T acc_a1 = a[in_idx];in_idx += in_stride;
    T acc_a2 = a[in_idx];in_idx += in_stride;
    T acc_a3 = a[in_idx];in_idx += in_stride;
    T acc_a4 = a[in_idx];in_idx += in_stride;

    __syncthreads();

    for (int ix = 0; ix < dimx; ix++)
    {
        // Advance the central point data
        acc_b4 = acc_b3;
        acc_b3 = acc_b2;
        acc_b2 = acc_b1;
        acc_b1 = acc;
        acc = acc_a1;
        acc_a1 = acc_a2;
        acc_a2 = acc_a3;
        acc_a3 = acc_a4;
        acc_a4 = a[in_idx];

        __syncthreads();
        slice[threadIdx.x*blockDim.y + threadIdx.y] = acc;
        __syncthreads();

        in_idx += in_stride;
        if(execute)
        {
                T tsum =      c.gmt4z*slice[ty*blockDim.y+tz-4] + c.gpt4z*slice[ty*blockDim.y+tz+4] +
                              c.gmt3z*slice[ty*blockDim.y+tz-3] + c.gpt3z*slice[ty*blockDim.y+tz+3] +
                              c.gmt2z*slice[ty*blockDim.y+tz-2] + c.gpt2z*slice[ty*blockDim.y+tz+2] +
                              c.gmt1z*slice[ty*blockDim.y+tz-1] + c.gpt1z*slice[ty*blockDim.y+tz+1] +
                              c.gmt4y*slice[(ty-4)*blockDim.y+tz] + c.gpt4y*slice[(ty+4)*blockDim.y+tz] +
                              c.gmt3y*slice[(ty-3)*blockDim.y+tz] + c.gpt3y*slice[(ty+3)*blockDim.y+tz] +
                              c.gmt2y*slice[(ty-2)*blockDim.y+tz] + c.gpt2y*slice[(ty+2)*blockDim.y+tz] +
                              c.gmt1y*slice[(ty-1)*blockDim.y+tz] + c.gpt1y*slice[(ty+1)*blockDim.y+tz];
                // Write back the results
                b[out_idx] = tsum +
                              c.gmt4x * acc_b4 +
                              c.gmt3x * acc_b3 +
                              c.gmt2x * acc_b2 +
                              c.gmt1x * acc_b1 +
                              c.a0 * acc +
                              c.gpt1x * acc_a1 +
                              c.gpt2x * acc_a2 +
                              c.gpt3x * acc_a3 +
                              c.gpt4x * acc_a4;
        }
        out_idx += dimy*dimz;
    }                   /* end for */

}
#endif

static std::vector<double *> abufs;
static std::vector<double *> bbufs;

void init_hip_fd(int max_threads, size_t bufsize)
{
    abufs.resize(max_threads);
    bbufs.resize(max_threads);

    for(int i=0;i < max_threads;i++)
    {
        hipMalloc((void **)&abufs[i], bufsize);
        hipMalloc((void **)&bbufs[i], bufsize);
    }
}

template void app8_del2_gpu(float * , float *, int, int, int, const fdparms_o8<float> &);
template void app8_del2_gpu(double * , double *, int, int, int, const fdparms_o8<double> &);
template void app8_del2_gpu(std::complex<float> * , std::complex<float> *, int, int, int, const fdparms_o8<std::complex<float>> &);
template void app8_del2_gpu(std::complex<double> * , std::complex<double> *, int, int, int, const fdparms_o8<std::complex<double>> &);



template <typename T>
void app8_del2_gpu(T *a, 
                   T *b, 
                   const int dimx,
                   const int dimy,
                   const int dimz,
                   const fdparms_o8<T> &c)
{
    dim3 Grid, Block;
    hipStream_t stream = getGpuStream();
    int tid = getThreadId();
    std::vector<int> yf, zf;
    GetPrimeFactors(yf, dimy, 19);
    GetPrimeFactors(zf, dimy, 19);

    // Block and grid layout can definitely use some optimization
    int xchunk = yf.back();
    if(yf[0] == 2 && yf[1] == 2) xchunk = 4;
    if(yf[0] == 2 && yf[1] == 3) xchunk = 6;
    int ychunk = zf.back();
    if(zf[0] == 2 && zf[1] == 2) ychunk = 4;
    if(zf[0] == 2 && zf[1] == 3) ychunk = 6;
    Grid.x = dimy / xchunk;
    Block.x = xchunk + 2*IMAGES;
    Grid.y = dimz / ychunk;
    Block.y = ychunk+2*IMAGES;
    int smem_siz = Block.x*Block.y*sizeof(T);


    hipStreamSynchronize(stream);
    hipLaunchKernelGGL(HIP_KERNEL_NAME(app8_del2_kernel<T>), Grid, Block, smem_siz, stream,
               a, b, dimx, dimy, dimz, c);
    hipStreamSynchronize(stream);
    return;

}
