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


#include <complex>
#include <hip/hip_ext.h>
#include "Gpufuncs.h"

__global__ void MulVec(double *dx, hipDoubleComplex *dy, int n)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    for (int i = idx; i < n; i += gridDim.x * blockDim.x) dy[i] = make_hipDoubleComplex(hipCreal(dy[i]) * dx[i], hipCimag(dy[i]) * dx[i]);
}

__global__ void MulVec1(double *dx, hipFloatComplex *dy, int n)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    for (int i = idx; i < n; i += gridDim.x * blockDim.x) dy[i] = make_hipFloatComplex(hipCrealf(dy[i]) * dx[i], hipCimagf(dy[i]) * dx[i]);
}

void GpuEleMul(double *dx, std::complex<double> *dy, int n, hipStream_t stream)
{
    hipStreamSynchronize(stream);
    int blockSize = 128;
    int numBlocks = (n + blockSize - 1) / n;
    hipLaunchKernelGGL(MulVec, dim3(numBlocks), dim3(blockSize), 0, stream, dx, (hipDoubleComplex *)dy, n);
    hipStreamSynchronize(stream);
}

void GpuEleMul(double *dx, std::complex<float> *dy, int n, hipStream_t stream)
{
    hipStreamSynchronize(stream);
    int blockSize = 128;
    int numBlocks = (n + blockSize - 1) / n;
    hipLaunchKernelGGL(MulVec1, dim3(numBlocks), dim3(blockSize), 0, stream, dx, (hipFloatComplex *)dy, n);
    hipStreamSynchronize(stream);
}
  

#endif
