/*
 *
 * Copyright 2026 The RMG Project Developers. See the COPYRIGHT file 
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
#include "Gpufuncs.h"
#include <complex>

__global__ void double_to_float(double* d_in, float* d_out, size_t n) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < n) {
        d_out[idx] = (float)d_in[idx];
    }
}

__global__ void float_to_double(float* d_in, double* d_out, size_t n) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < n) {
        d_out[idx] = (double)d_in[idx];
    }
}

__global__ void pack_diagonal(double *matrix, double *diagonal, size_t N) {
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < N) {
        diagonal[i] = matrix[i * (N + 1)];
    }
}

__global__ void scatter_diagonal(double *matrix, double *diagonal, size_t N) {
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < N) {
        matrix[i * (N + 1)] = diagonal[i];
    }
}

// HIP kernel using 1D mapping for 100% thread utilization
__global__ void packUpperTriangular(double *matrix, float *packed, size_t N) {
    size_t total_dense_elements = N * N;

    for (size_t idx = (size_t)blockIdx.x * blockDim.x + threadIdx.x; 
         idx < total_dense_elements; 
         idx += (size_t)gridDim.x * blockDim.x) {
        
        size_t i = idx % N; // Row index (swapped for column-major)
        size_t j = idx / N; // Column index

        if (i <= j) { // Upper triangular condition
            // Exact BLAS column-packed index calculation
            size_t packed_idx = (j * (j + 1)) / 2 + i;
            packed[packed_idx] = (double)matrix[idx];
        }
    }
}


__global__ void unpackUpperTriangular(double *matrix, float *packed, size_t N) {
    size_t total_dense_elements = N * N;

    for (size_t idx = (size_t)blockIdx.x * blockDim.x + threadIdx.x; 
         idx < total_dense_elements; 
         idx += (size_t)gridDim.x * blockDim.x) {
        
        size_t i = idx % N;
        size_t j = idx / N;

        if (i <= j) {
            size_t packed_idx = (j * (j + 1)) / 2 + i;
            matrix[idx] = (double)packed[packed_idx];
        } else {
            //A[idx] = 0.0f;
        }
    }
}

void gpuPackSqToUpTr(double *matrix, float *packed, size_t N)
{
    int blockSize = 256;
    
    // Bounded grid size safely handles any value of N up to hardware limits
    int gridSize = 4096; 

    hipLaunchKernelGGL(packUpperTriangular, dim3(gridSize), dim3(blockSize), 0, 0, matrix, packed, N);
    //hipDeviceSynchronize();
}

void gpuUnPackTrToSq(double *matrix, float *packed, size_t N)
{
    int blockSize = 256;
    
    // Bounded grid size safely handles any value of N up to hardware limits
    int gridSize = 4096;
            
    hipLaunchKernelGGL(unpackUpperTriangular, dim3(gridSize), dim3(blockSize), 0, 0, matrix, packed, N);
    //hipDeviceSynchronize();
}

void gpu_pack_diagonal(double *matrix, double *diagonal, size_t n)
{
    int block_size = 256;
    int grid_size = (n + block_size - 1) / block_size;
    pack_diagonal<<<grid_size, block_size>>>(matrix, diagonal, n);
}

void gpu_scatter_diagonal(double *matrix, double *diagonal, size_t n)
{
    int block_size = 256;
    int grid_size = (n + block_size - 1) / block_size;
    scatter_diagonal<<<grid_size, block_size>>>(matrix, diagonal, n);
}

void gpu_copy_and_convert(double* d_in, float* d_out, size_t n) {
    int block_size = 256;
    int grid_size = (n + block_size - 1) / block_size;

    double_to_float<<<grid_size, block_size>>>(d_in, d_out, n);
}

void gpu_copy_and_convert(float* d_in, double* d_out, size_t n) {
    int block_size = 256;
    int grid_size = (n + block_size - 1) / block_size;

    float_to_double<<<grid_size, block_size>>>(d_in, d_out, n);
}

void gpu_copy_and_convert(std::complex<double> *d_in, std::complex<float> *d_out, size_t n) {
    int block_size = 256;
    int grid_size = (n + block_size - 1) / block_size;
    grid_size *= 2;

    double_to_float<<<grid_size, block_size>>>((double *)d_in, (float *)d_out, 2*n);
}

void gpu_copy_and_convert(std::complex<float> *d_in, std::complex<double> *d_out, size_t n) {
    int block_size = 256;
    int grid_size = (n + block_size - 1) / block_size;
    grid_size *= 2;

    float_to_double<<<grid_size, block_size>>>((float *)d_in, (double *)d_out, 2*n);
}

#endif

