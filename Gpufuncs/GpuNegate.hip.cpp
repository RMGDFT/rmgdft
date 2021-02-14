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

#include "Gpufuncs.h"

__global__ void NegateDiag(double *dx, int incx, double *dy, int incy, int n)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    for (int i = idx; i < n; i += gridDim.x * blockDim.x) 
    {
        if(dx[i*incx + i] < 0.0) dy[i*incy] = -dy[i*incy];
    }
}


void GpuNegate(double *dx, int incx, double *dy, int incy, int n)
{
    int blockSize = 256;
    int numBlocks = (n + blockSize - 1) / n;
    hipLaunchKernelGGL(NegateDiag, dim3(numBlocks), dim3(blockSize), 0, 0, dx, incx, dy, incy, n);
}

#endif
