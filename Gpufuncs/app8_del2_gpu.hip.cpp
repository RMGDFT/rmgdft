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

template <typename T, int BLOCKX, int BLOCKY, int n>
__global__ void app8_del2_kernel(const T * __restrict__ a, 
                                                T *b, 
                                                const int dimx,
                                                const int dimy,
                                                const int dimz,
                                                const T h2x,
                                                const T h2y,
                                                const T h2z)

{

    __shared__ double aslice[BLOCKX + 2*IMAGES][BLOCKY + 2*IMAGES];
//    extern __shared__ __align__(sizeof(T)) unsigned char sbuf[];
    T *slice = reinterpret_cast<T *>(&aslice[0][0]);
    int islice = BLOCKY + 2*IMAGES;
    int limit = dimx * dimy * dimz;

    T a1(205.0/72.0);
    T ONE_T(1.0);
    T EIGHT_T(8.0);
    T FIVE_T(5.0);
    T T15_T(315.0);
    T T560_T(560.0);

    T c0 = -(a1) * (h2x + h2y + h2z);
    T c1x = EIGHT_T * h2x / FIVE_T;
    T c2x = -ONE_T * h2x / FIVE_T;;
    T c3x = EIGHT_T * h2x / T15_T;
    T c4x = -ONE_T * h2x / T560_T;

    T c1y = EIGHT_T * h2y / FIVE_T;
    T c2y = -ONE_T * h2y / FIVE_T;
    T c3y = EIGHT_T * h2y / T15_T;
    T c4y = -ONE_T * h2y / T560_T;

    T c1z = EIGHT_T * h2z / FIVE_T;
    T c2z = -ONE_T * h2z / FIVE_T;
    T c3z = EIGHT_T * h2z / T15_T;
    T c4z = -ONE_T * h2z / T560_T;

    // iy and iz map to the x and y coordinates of the thread
    // within a block
    int iy = blockIdx.x * blockDim.x + threadIdx.x;
    int iz = blockIdx.y * blockDim.y + threadIdx.y;

    // thread y-index into shared memory tile
    int ty = threadIdx.x + IMAGES;
    int tz = threadIdx.y + IMAGES;

    int in_idx = (iy + IMAGES)*(dimz + 2*IMAGES) + iz + IMAGES;
    int in_stride = (dimy + 2*IMAGES) * (dimz + 2*IMAGES);
    int out_idx = iy*dimz + iz;
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

        // Update the data slice in shared memory
        if(threadIdx.x < IMAGES)
        {
            slice[threadIdx.x*islice + tz] = a[in_idx - IMAGES*in_stride - IMAGES*(dimz+2*IMAGES)];
            slice[(threadIdx.x+blockDim.x+IMAGES)*islice + tz] = a[in_idx - IMAGES*in_stride + blockDim.x*(dimz+2*IMAGES)];
        }
        if(threadIdx.y < IMAGES)
        {
            slice[ty*islice + threadIdx.y] = a[in_idx - IMAGES*in_stride - IMAGES];
            slice[ty*islice + threadIdx.y+blockDim.y+IMAGES] = a[in_idx - IMAGES*in_stride + blockDim.y];
        }

        slice[ty*islice + tz] = acc;
        __syncthreads();

        in_idx += in_stride;

        T tsum =      c4z*slice[ty*islice+tz-4] + c4z*slice[ty*islice+tz+4] +
                      c3z*slice[ty*islice+tz-3] + c3z*slice[ty*islice+tz+3] +
                      c2z*slice[ty*islice+tz-2] + c2z*slice[ty*islice+tz+2] +
                      c1z*slice[ty*islice+tz-1] + c1z*slice[ty*islice+tz+1] +
                      c4y*slice[(ty-4)*islice+tz] + c4y*slice[(ty+4)*islice+tz] +
                      c3y*slice[(ty-3)*islice+tz] + c3y*slice[(ty+3)*islice+tz] +
                      c2y*slice[(ty-2)*islice+tz] + c2y*slice[(ty+2)*islice+tz] +
                      c1y*slice[(ty-1)*islice+tz] + c1y*slice[(ty+1)*islice+tz];

        // Write back the results

        b[out_idx] = tsum +
                      c4x * acc_b4 +
                      c3x * acc_b3 +
                      c2x * acc_b2 +
                      c1x * acc_b1 +
                      c0 * acc +
                      c1x * acc_a1 +
                      c2x * acc_a2 +
                      c3x * acc_a3 +
                      c4x * acc_a4;
        out_idx += out_stride;
    }                   /* end for */

}


template double app8_del2_gpu(const float * , float *, int, int, int, float, float, float, hipStream_t);
template double app8_del2_gpu(const double * , double *, int, int, int, double, double, double, hipStream_t);
template double app8_del2_gpu(const std::complex<float> * , std::complex<float> *, int, int, int, std::complex<float>, std::complex<float>, std::complex<float>, hipStream_t);
template double app8_del2_gpu(const std::complex<double> * , std::complex<double> *, int, int, int, std::complex<double>, std::complex<double>, std::complex<double>, hipStream_t);

template <typename T>
double app8_del2_gpu(const T * __restrict__ a, 
                   T *b, 
                   const int dimx,
                   const int dimy,
                   const int dimz,
                   T h2x,
                   T h2y,
                   T h2z,
                   hipStream_t cstream)
{
    dim3 Grid, Block;
    T ONE_T(1.0);
    T fac(205.0/72.0);
h2x = ONE_T/h2x;
h2y = ONE_T/h2y;
h2z = ONE_T/h2z;
    double retval = -std::real(fac * (h2x + h2y + h2z));

    if(!(dimy % 16) && !(dimz % 32))
    {
        Grid.x = dimy / 16;
        Block.x = 16;
        Grid.y = dimz / 32;
        Block.y = 32;
        hipLaunchKernelGGL(HIP_KERNEL_NAME(app8_del2_kernel<T, 16, 32, 40*40*8 >), Grid, Block, 0, cstream, a, b, dimx, dimy, dimz, h2x, h2y, h2z);
    }
    if(!(dimy % 4) && !(dimz % 64))
    {
        Grid.x = dimy / 4;
        Block.x = 4;
        Grid.y = dimz / 64;
        Block.y = 64;
        hipLaunchKernelGGL(HIP_KERNEL_NAME(app8_del2_kernel<T, 4, 64, 12*72*8 >), Grid, Block, 0, cstream, a, b, dimx, dimy, dimz, h2x, h2y, h2z);
        return retval;
    }
#if 0
    if(!(dimy % 8) && !(dimz % 64))
    {
        Grid.x = dimy / 8;
        Block.x = 8;
        Grid.y = dimz / 64;
        Block.y = 64;
        hipLaunchKernelGGL(HIP_KERNEL_NAME(app8_del2_kernel<T, 8, 64>), Grid, Block, 0, cstream, a, b, dimx, dimy, dimz, h2x, h2y, h2z);
        return retval;
    }
    if(!(dimy % 8) && !(dimz % 32))
    {
        Grid.x = dimy / 8;
        Block.x = 8;
        Grid.y = dimz / 32;
        Block.y = 32;
        hipLaunchKernelGGL(HIP_KERNEL_NAME(app8_del2_kernel<T, 8, 32>), Grid, Block, 0, cstream, a, b, dimx, dimy, dimz, h2x, h2y, h2z);
        return retval;
    }
    if(!(dimy % 16) && !(dimz % 16))
    {
        Grid.x = dimy / 16;
        Block.x = 16;
        Grid.y = dimz / 16;
        Block.y = 16;
        hipLaunchKernelGGL(HIP_KERNEL_NAME(app8_del2_kernel<T, 16, 16>), Grid, Block, 0, cstream, a, b, dimx, dimy, dimz, h2x, h2y, h2z);
        return retval;
    }
    if(!(dimy % 16) && !(dimz % 24))
    {
        Grid.x = dimy / 16;
        Block.x = 16;
        Grid.y = dimz / 24;
        Block.y = 24;
        hipLaunchKernelGGL(HIP_KERNEL_NAME(app8_del2_kernel<T, 16, 24>), Grid, Block, 0, cstream, a, b, dimx, dimy, dimz, h2x, h2y, h2z);
        return retval;
    }
    if(!(dimy % 24) && !(dimz % 16))
    {
        Grid.x = dimy / 24;
        Block.x = 24;
        Grid.y = dimz / 16;
        Block.y = 16;
        hipLaunchKernelGGL(HIP_KERNEL_NAME(app8_del2_kernel<T, 24, 16>), Grid, Block, 0, cstream, a, b, dimx, dimy, dimz, h2x, h2y, h2z);
        return retval;
    }
    if(!(dimy % 24) && !(dimz % 24))
    {
        Grid.x = dimy / 24;
        Block.x = 24;
        Grid.y = dimz / 24;
        Block.y = 24;
        hipLaunchKernelGGL(HIP_KERNEL_NAME(app8_del2_kernel<T, 24, 24>), Grid, Block, 0, cstream, a, b, dimx, dimy, dimz, h2x, h2y, h2z);
        return retval;
    }
    if(!(dimy % 10) && !(dimz % 10))
    {
        Grid.x = dimy / 10;
        Block.x = 10;
        Grid.y = dimz / 10;
        Block.y = 10;
        hipLaunchKernelGGL(HIP_KERNEL_NAME(app8_del2_kernel<T, 10, 10>), Grid, Block, 0, cstream, a, b, dimx, dimy, dimz, h2x, h2y, h2z);
        return retval;
    }
    if(!(dimy % 8) && !(dimz % 8))
    {
        Grid.x = dimy / 8;
        Block.x = 8;
        Grid.y = dimz / 8;
        Block.y = 8;
        hipLaunchKernelGGL(HIP_KERNEL_NAME(app8_del2_kernel<T, 8, 8>), Grid, Block, 0, cstream, a, b, dimx, dimy, dimz, h2x, h2y, h2z);
        return retval;
    }
    if(!(dimy % 12) && !(dimz % 12))
    {
        Grid.x = dimy / 12;
        Block.x = 12;
        Grid.y = dimz / 12;
        Block.y = 12;
        hipLaunchKernelGGL(HIP_KERNEL_NAME(app8_del2_kernel<T, 12, 12>), Grid, Block, 0, cstream, a, b, dimx, dimy, dimz, h2x, h2y, h2z);
        return retval;
    }
#endif
    if(!(dimy % 15) && !(dimz % 15))
    {
        Grid.x = dimy / 15;
        Block.x = 15;
        Grid.y = dimz / 15;
        Block.y = 15;
        hipLaunchKernelGGL(HIP_KERNEL_NAME(app8_del2_kernel<T, 15, 35, 24*24*8 >), Grid, Block, 0, cstream, a, b, dimx, dimy, dimz, h2x, h2y, h2z);
        return retval;
    }
//printf("Missed  %d  %d\n", dimy, dimz);
    return retval;
}
#endif
