#include "hip/hip_runtime.h"
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

#define IMAGES 4


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
    int out_stride = dimy*dimz;

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

        if((threadIdx.x >= IMAGES) &&
           (threadIdx.y >= IMAGES) &&
           (threadIdx.x < (blockDim.x-IMAGES)) &&
           (threadIdx.y < (blockDim.y-IMAGES)))
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
        //__syncthreads();
        out_idx += out_stride;
    }                   /* end for */

}


template double app8_del2_gpu(const float * , float *, int, int, int, double, double, double, hipStream_t);
template double app8_del2_gpu(const double * , double *, int, int, int, double, double, double, hipStream_t);
template double app8_del2_gpu(const std::complex<float> * , std::complex<float> *, int, int, int, double, double, double, hipStream_t);
template double app8_del2_gpu(const std::complex<double> * , std::complex<double> *, int, int, int, double, double, double, hipStream_t);

template <typename T>
double app8_del2_gpu(const T * __restrict__ a, 
                   T *b, 
                   const int dimx,
                   const int dimy,
                   const int dimz,
                   double h2x,
                   double h2y,
                   double h2z,
                   hipStream_t cstream)
{
    dim3 Grid, Block;
    T ONE_T(1.0);
    T fac(205.0/72.0);
    T h2x_t(h2x);
    T h2y_t(h2y);
    T h2z_t(h2z);
    h2x_t = ONE_T/h2x_t;
    h2y_t = ONE_T/h2y_t;
    h2z_t = ONE_T/h2z_t;
    double retval = -std::real(fac * (h2x_t + h2y_t + h2z_t));

    fdparms_o8<T> c;
    T a0 = -T(205.0/72.0) * (h2x_t + h2y_t + h2z_t);
    T c1x = T(8.0) * h2x_t / T(5.0);
    T c2x = -T(1.0) * h2x_t / T(5.0);
    T c3x = T(8.0) * h2x_t / T(315.0);
    T c4x = -T(1.0) * h2x_t / T(560.0);

    T c1y = T(8.0) * h2y_t / T(5.0);
    T c2y = -T(1.0) * h2y_t / T(5.0);
    T c3y = T(8.0) * h2y_t / T(315.0);
    T c4y = -T(1.0) * h2y_t / T(560.0);

    T c1z = T(8.0) * h2z_t / T(5.0);
    T c2z = -T(1.0) * h2z_t / T(5.0);
    T c3z = T(8.0) * h2z_t / T(315.0);
    T c4z = -T(1.0) * h2z_t / T(560.0);

    c.a0 = a0;
    c.gpt1x = c.gmt1x = c1x;
    c.gpt2x = c.gmt2x = c2x;
    c.gpt3x = c.gmt3x = c3x;
    c.gpt4x = c.gmt4x = c4x;

    c.gpt1y = c.gmt1y = c1y;
    c.gpt2y = c.gmt2y = c2y;
    c.gpt3y = c.gmt3y = c3y;
    c.gpt4y = c.gmt4y = c4y;

    c.gpt1z = c.gmt1z = c1z;
    c.gpt2z = c.gmt2z = c2z;
    c.gpt3z = c.gmt3z = c3z;
    c.gpt4z = c.gmt4z = c4z;

    Grid.x = dimy/8;
    Block.x = 16;
    Grid.y = dimz / 8;
    Block.y = 16;
    int smem_siz = (Block.x + 2*IMAGES)*(Block.y + 2*IMAGES)*sizeof(T);
    hipLaunchKernelGGL(HIP_KERNEL_NAME(app8_del2_kernel<T>), Grid, Block, smem_siz, cstream, a, b, dimx, dimy, dimz,c);

    return retval;
}
#endif
