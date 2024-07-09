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
#include <hip/hip_ext.h>
#include <complex>
#include "Gpufuncs.h"

template <typename RmgType>
__global__ void FillMatrix(RmgType *matrix, int n, RmgType val)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    for (int i = idx; i < n; i += gridDim.x * blockDim.x) matrix[i]=val;
}


// Fills a device pointer with a floating point value
void GpuFill(double *dptr, int n, double fillval)
{
    hipStream_t stream = getGpuStream();

    //hipDeviceSynchronize();
    int blockSize = 256;
    int numBlocks = n / blockSize;
    if(n % blockSize) numBlocks++;
    //hipLaunchKernelGGL(FillMatrix, dim3(numBlocks), dim3(blockSize), 0, 0, dptr, n, fillval);
    FillMatrix<<<numBlocks, blockSize, 0, stream>>>(dptr, n, fillval);

    //hipDeviceSynchronize();
}

#endif
