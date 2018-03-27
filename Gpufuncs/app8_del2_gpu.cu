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


#if GPU_ENABLED

#include <cuda.h>
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>
#include <cuda_device_runtime_api.h>
#include <crt/host_runtime.h>
#include <stdio.h>

#define IMAGES 4

template <int BLOCKX, int BLOCKY>
__global__ void app8_del2_kernel(const double * __restrict__ a, 
                                                double *b, 
                                                const int dimx,
                                                const int dimy,
                                                const int dimz,
                                                const double h2x,
                                                const double h2y,
                                                const double h2z)

{

    __shared__ double slice[BLOCKX + 2*IMAGES][BLOCKY + 2*IMAGES];
    double c0 = -(205.0 / 72.0) * (h2x + h2y + h2z);
    double c1x = 8.0 * h2x / 5.0;
    double c2x = -1.0 * h2x / 5.0;
    double c3x = 8.0 * h2x / 315.0;
    double c4x = -1.0 * h2x / 560.0;
    double c1y = 8.0 * h2y / 5.0;
    double c2y = -1.0 * h2y / 5.0;
    double c3y = 8.0 * h2y / 315.0;
    double c4y = -1.0 * h2y / 560.0;

    double c1z = 8.0 * h2z / 5.0;
    double c2z = -1.0 * h2z / 5.0;
    double c3z = 8.0 * h2z / 315.0;
    double c4z = -1.0 * h2z / 560.0;

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

    double acc_b4 = 0.0;
    double acc_b3 = a[in_idx];in_idx += in_stride;
    double acc_b2 = a[in_idx];in_idx += in_stride;
    double acc_b1 = a[in_idx];in_idx += in_stride;
    double acc    = a[in_idx];in_idx += in_stride;
    double acc_a1 = a[in_idx];in_idx += in_stride;
    double acc_a2 = a[in_idx];in_idx += in_stride;
    double acc_a3 = a[in_idx];in_idx += in_stride;
    double acc_a4 = a[in_idx];in_idx += in_stride;

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
            slice[threadIdx.x][tz] = a[in_idx - IMAGES*in_stride - IMAGES*(dimz+2*IMAGES)];
            slice[threadIdx.x+blockDim.x+IMAGES][tz] = a[in_idx - IMAGES*in_stride + blockDim.x*(dimz+2*IMAGES)];
        }
        if(threadIdx.y < IMAGES)
        {
            slice[ty][threadIdx.y] = a[in_idx - IMAGES*in_stride - IMAGES];
            slice[ty][threadIdx.y+blockDim.y+IMAGES] = a[in_idx - IMAGES*in_stride + blockDim.y];
        }

        slice[ty][tz] = acc;
        __syncthreads();

        in_idx += in_stride;

        double tsum = c4z*slice[ty][tz-4] + c4z*slice[ty][tz+4] +
                      c3z*slice[ty][tz-3] + c3z*slice[ty][tz+3] +
                      c2z*slice[ty][tz-2] + c2z*slice[ty][tz+2] +
                      c1z*slice[ty][tz-1] + c1z*slice[ty][tz+1] +
                      c4y*slice[ty-4][tz] + c4y*slice[ty+4][tz] +
                      c3y*slice[ty-3][tz] + c3y*slice[ty+3][tz] +
                      c2y*slice[ty-2][tz] + c2y*slice[ty+2][tz] +
                      c1y*slice[ty-1][tz] + c1y*slice[ty+1][tz];

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


double app8_del2_gpu(const double * __restrict__ a, 
                   double *b, 
                   const int dimx,
                   const int dimy,
                   const int dimz,
                   double h2x,
                   double h2y,
                   double h2z,
                   cudaStream_t cstream)
{
    dim3 Grid, Block;
    double retval = -(205.0 / 72.0) * (h2x + h2y + h2z);

    if(!(dimy % 16) && !(dimz % 32))
    {
        Grid.x = dimy / 16;
        Block.x = 16;
        Grid.y = dimz / 32;
        Block.y = 32;
        app8_del2_kernel<16, 32><<<Grid, Block, 0, cstream>>>(a, b, dimx, dimy, dimz, h2x, h2y, h2z);
        return retval;
    }
    if(!(dimy % 4) && !(dimz % 64))
    {
        Grid.x = dimy / 4;
        Block.x = 4;
        Grid.y = dimz / 64;
        Block.y = 64;
        app8_del2_kernel<4, 64><<<Grid, Block, 0, cstream>>>(a, b, dimx, dimy, dimz, h2x, h2y, h2z);
        return retval;
    }
    if(!(dimy % 8) && !(dimz % 64))
    {
        Grid.x = dimy / 8;
        Block.x = 8;
        Grid.y = dimz / 64;
        Block.y = 64;
        app8_del2_kernel<8, 64><<<Grid, Block, 0, cstream>>>(a, b, dimx, dimy, dimz, h2x, h2y, h2z);
        return retval;
    }
    if(!(dimy % 8) && !(dimz % 32))
    {
        Grid.x = dimy / 8;
        Block.x = 8;
        Grid.y = dimz / 32;
        Block.y = 32;
        app8_del2_kernel<8, 32><<<Grid, Block, 0, cstream>>>(a, b, dimx, dimy, dimz, h2x, h2y, h2z);
        return retval;
    }
    if(!(dimy % 16) && !(dimz % 16))
    {
        Grid.x = dimy / 16;
        Block.x = 16;
        Grid.y = dimz / 16;
        Block.y = 16;
        app8_del2_kernel<16, 16><<<Grid, Block, 0, cstream>>>(a, b, dimx, dimy, dimz, h2x, h2y, h2z);
        return retval;
    }
    if(!(dimy % 16) && !(dimz % 24))
    {
        Grid.x = dimy / 16;
        Block.x = 16;
        Grid.y = dimz / 24;
        Block.y = 24;
        app8_del2_kernel<16, 24><<<Grid, Block, 0, cstream>>>(a, b, dimx, dimy, dimz, h2x, h2y, h2z);
        return retval;
    }
    if(!(dimy % 24) && !(dimz % 16))
    {
        Grid.x = dimy / 24;
        Block.x = 24;
        Grid.y = dimz / 16;
        Block.y = 16;
        app8_del2_kernel<24, 16><<<Grid, Block, 0, cstream>>>(a, b, dimx, dimy, dimz, h2x, h2y, h2z);
        return retval;
    }
    if(!(dimy % 24) && !(dimz % 24))
    {
        Grid.x = dimy / 24;
        Block.x = 24;
        Grid.y = dimz / 24;
        Block.y = 24;
        app8_del2_kernel<24, 24><<<Grid, Block, 0, cstream>>>(a, b, dimx, dimy, dimz, h2x, h2y, h2z);
        return retval;
    }
    if(!(dimy % 10) && !(dimz % 10))
    {
        Grid.x = dimy / 10;
        Block.x = 10;
        Grid.y = dimz / 10;
        Block.y = 10;
        app8_del2_kernel<10, 10><<<Grid, Block, 0, cstream>>>(a, b, dimx, dimy, dimz, h2x, h2y, h2z);
        return retval;
    }
    if(!(dimy % 8) && !(dimz % 8))
    {
        Grid.x = dimy / 8;
        Block.x = 8;
        Grid.y = dimz / 8;
        Block.y = 8;
        app8_del2_kernel<8, 8><<<Grid, Block, 0, cstream>>>(a, b, dimx, dimy, dimz, h2x, h2y, h2z);
        return retval;
    }
    if(!(dimy % 12) && !(dimz % 12))
    {
        Grid.x = dimy / 12;
        Block.x = 12;
        Grid.y = dimz / 12;
        Block.y = 12;
        app8_del2_kernel<12, 12><<<Grid, Block, 0, cstream>>>(a, b, dimx, dimy, dimz, h2x, h2y, h2z);
        return retval;
    }

    return retval;
}
#endif
