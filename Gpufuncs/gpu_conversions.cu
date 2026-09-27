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

#if CUDA_ENABLED
#include "Gpufuncs.h"
#include <complex>

__global__ void double_to_float(const double* d_in, float* d_out, int n) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < n) {
        d_out[idx] = (float)d_in[idx];
    }
}

__global__ void float_to_double(const float* d_in, double* d_out, int n) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < n) {
        d_out[idx] = (double)d_in[idx];
    }
}


void gpu_copy_and_convert(const double* d_in, float* d_out, int n) {
    int block_size = 256;
    int grid_size = (n + block_size - 1) / block_size;

    double_to_float<<<grid_size, block_size>>>(d_in, d_out, n);
}

void gpu_copy_and_convert(const float* d_in, double* d_out, int n) {
    int block_size = 256;
    int grid_size = (n + block_size - 1) / block_size;

    float_to_double<<<grid_size, block_size>>>(d_in, d_out, n);
}

void gpu_copy_and_convert(const std::complex<double> *d_in, std::complex<float> *d_out, int n) {
    int block_size = 256;
    int grid_size = (n + block_size - 1) / block_size;
    grid_size *= 2;

    double_to_float<<<grid_size, block_size>>>((double *)d_in, (float *)d_out, 2*n);
}

void gpu_copy_and_convert(const std::complex<float> *d_in, std::complex<double> *d_out, int n) {
    int block_size = 256;
    int grid_size = (n + block_size - 1) / block_size;
    grid_size *= 2;

    float_to_double<<<grid_size, block_size>>>((float *)d_in, (double *)d_out, 2*n);
}

#endif

