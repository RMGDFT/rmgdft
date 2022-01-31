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

void GetPrimeFactors(std::vector<int>& factors, int val, int stop);

#define IMAGES 4


#if 0
// Version without shared memory slice
template <typename T>
__global__ void app8_del2_kernel(T * __restrict__ a, 
                                                T *b, 
                                                const int dimx,
                                                const int dimy,
                                                const int dimz,
                                                const fdparms_o8<T> c)

{


    // iy and iz map to the x and y coordinates of the thread
    // within a block
    int iy = blockIdx.x * blockDim.x + threadIdx.x + IMAGES;
    int iz = blockIdx.y * blockDim.y + threadIdx.y + IMAGES;

    int in_idx = iy*(dimz + 2*IMAGES) + iz;
    int in_stride = (dimy + 2*IMAGES) * (dimz + 2*IMAGES);

    int out_idx = (iy-IMAGES)*dimz + iz-IMAGES;
    int ys = dimz + 2*IMAGES;
    T *p=a;
    T *pp=&a[IMAGES*in_stride];

    T acc_b4 = 0.0;
    T acc_b3 = p[in_idx];in_idx += in_stride;
    T acc_b2 = p[in_idx];in_idx += in_stride;
    T acc_b1 = p[in_idx];in_idx += in_stride;
    T acc    = p[in_idx];in_idx += in_stride;
    T acc_a1 = p[in_idx];in_idx += in_stride;
    T acc_a2 = p[in_idx];in_idx += in_stride;
    T acc_a3 = p[in_idx];in_idx += in_stride;
    T acc_a4 = p[in_idx];in_idx += in_stride;
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
        acc_a4 = p[in_idx];
        in_idx += in_stride;

                T tsum =      c.gmt4z*pp[iy*ys+iz-4] + c.gpt4z*pp[iy*ys+iz+4] +
                              c.gmt3z*pp[iy*ys+iz-3] + c.gpt3z*pp[iy*ys+iz+3] +
                              c.gmt2z*pp[iy*ys+iz-2] + c.gpt2z*pp[iy*ys+iz+2] +
                              c.gmt1z*pp[iy*ys+iz-1] + c.gpt1z*pp[iy*ys+iz+1] +
                              c.gmt4y*pp[(iy-4)*ys+iz] + c.gpt4y*pp[(iy+4)*ys+iz] +
                              c.gmt3y*pp[(iy-3)*ys+iz] + c.gpt3y*pp[(iy+3)*ys+iz] +
                              c.gmt2y*pp[(iy-2)*ys+iz] + c.gpt2y*pp[(iy+2)*ys+iz] +
                              c.gmt1y*pp[(iy-1)*ys+iz] + c.gpt1y*pp[(iy+1)*ys+iz];
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
        out_idx += dimy*dimz;
        pp += in_stride;
    }                   /* end for */
}

#else
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

    T acc_b4 = 0.0;
    T acc_b3 = a[in_idx];in_idx += in_stride;
    T acc_b2 = a[in_idx];in_idx += in_stride;
    T acc_b1 = a[in_idx];in_idx += in_stride;
    T acc    = a[in_idx];in_idx += in_stride;
    T acc_a1 = a[in_idx];in_idx += in_stride;
    T acc_a2 = a[in_idx];in_idx += in_stride;
    T acc_a3 = a[in_idx];in_idx += in_stride;
    T acc_a4 = a[in_idx];in_idx += in_stride;

    // All threads load array to shared memory. This bool is
    // used to identify which threads compute and write results
    bool execute = ((threadIdx.x >= IMAGES) &&
           (threadIdx.y >= IMAGES) &&
           (threadIdx.x < (blockDim.x-IMAGES)) &&
           (threadIdx.y < (blockDim.y-IMAGES)));

    __syncthreads();
    T *p = &slice[threadIdx.x*blockDim.y + threadIdx.y];
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

// This section is for overlapping transfer/computation
#define MAX_CHUNKS 4
void *abufs[64];
void *bbufs[64];
static hipStream_t streams[64][MAX_CHUNKS];

void init_hip_fd(int max_threads)
{
    for(int i=0;i < max_threads;i++)
    {
        hipMalloc((void **)&abufs[i], 8000000);
        hipMalloc((void **)&bbufs[i], 8000000);
        for(int j=0;j < MAX_CHUNKS;j++)
        {
            hipStreamCreateWithPriority(&streams[i][j], hipStreamNonBlocking, i);
        } 
    }
}

template void app8_del2_gpu(float * , float *, int, int, int, const fdparms_o8<float> &, int);
template void app8_del2_gpu(double * , double *, int, int, int, const fdparms_o8<double> &, int);
template void app8_del2_gpu(std::complex<float> * , std::complex<float> *, int, int, int, const fdparms_o8<std::complex<float>> &, int);
template void app8_del2_gpu(std::complex<double> * , std::complex<double> *, int, int, int, const fdparms_o8<std::complex<double>> &, int);


void FD_barrier(void);

template <typename T>
void app8_del2_gpu(T *a, 
                   T *b, 
                   const int dimx,
                   const int dimy,
                   const int dimz,
                   const fdparms_o8<T> &c,
                   int tid)
{
    dim3 Grid, Block;
    std::vector<int> xf, yf;
    GetPrimeFactors(xf, dimx, 19);
    GetPrimeFactors(yf, dimy, 19);
#if 1
    Grid.x = dimy/8;
    Block.x = 16;
    Grid.y = dimz / 8;
    Block.y = 16;
    int smem_siz = Block.x*Block.y*sizeof(T);
#else
    Grid.x = dimy;
    Block.x = 1;
    Grid.y = 1;
    Block.y = dimz;
    int smem_siz = 256;
#endif

#if 0
    // Managed memory version
    hipLaunchKernelGGL(HIP_KERNEL_NAME(app8_del2_kernel<T>), Grid, Block, smem_siz, streams[tid][0],
               a, b, dimx, dimy, dimz, c);
    hipStreamSynchronize(streams[tid][0]);
    return;
#else

    int chunks  = 8;
    if(dimx % chunks) chunks = 1;
chunks=4;
    int chunk_in = (dimx/chunks)*(dimy+2*IMAGES)*(dimz+2*IMAGES);
    int chunk_out = (dimx/chunks)*dimy*dimz;
    int copy_in, offset_in, offset_out;

    T *d = (T *)abufs[tid];
    T *p = (T *)bbufs[tid];

    for(int i=0;i < chunks;i++) hipStreamSynchronize(streams[tid][i]);
     
    copy_in = chunk_in + 2*IMAGES*(dimy+2*IMAGES)*(dimz+2*IMAGES);
    hipMemcpyAsync(d, a, copy_in*sizeof(T), hipMemcpyHostToDevice, streams[tid][0]);
    offset_in = copy_in;
    for(int i=1;i < chunks;i++)
    {

        hipMemcpyAsync(&d[offset_in], &a[offset_in], chunk_in*sizeof(T), hipMemcpyHostToDevice, streams[tid][i]);
        offset_in += chunk_in;
    }

    hipLaunchKernelGGL(HIP_KERNEL_NAME(app8_del2_kernel<T>), Grid, Block, smem_siz, streams[tid][0],
               d, p, dimx/chunks, dimy, dimz, c);

    offset_in = chunk_in;
    offset_out = chunk_out;
    for(int i=1;i < chunks;i++)
    {
        hipStreamSynchronize(streams[tid][i-1]);
        hipMemcpyAsync(&b[(i-1)*chunk_out], &p[(i-1)*chunk_out], chunk_out*sizeof(T), hipMemcpyDeviceToHost, streams[tid][i-1]);
        hipLaunchKernelGGL(HIP_KERNEL_NAME(app8_del2_kernel<T>), Grid, Block, smem_siz, streams[tid][i],
               &d[offset_in], &p[offset_out], dimx/chunks, dimy, dimz, c);
        offset_in += chunk_in;
        offset_out += chunk_out;
    }
    hipMemcpyAsync(&b[(chunks-1)*chunk_out], &p[(chunks-1)*chunk_out], chunk_out*sizeof(T), hipMemcpyDeviceToHost, streams[tid][chunks-1]);

    for(int i=0;i < chunks;i++) hipStreamSynchronize(streams[tid][i]);

#endif

}
#endif
