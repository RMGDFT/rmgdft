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

#define GRAM_GPU_BLOCK 64

__global__ void gramsch_update_psi_kernel(double *G,
                                     const double * __restrict__ V,
                                     const double * __restrict__ C,
                                     int N,
                                     int eig_start,
                                     int eig_stop)
{

//    __shared__ double sarr[GRAM_GPU_BLOCK];
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
#if 0
    // Get inverse of diagonal elements
    for(int ix = 0;ix < n;ix++) tarr[ix] = 1.0 / C[n * ix + ix];

//----------------------------------------------------------------
//----------------------------------------------------------------

    for(int idx = eig_start;idx < eig_stop;idx++)
    {
        for (int st = 0; st < n; st++) 
        {
            double dtmp = V[st*n + idx] * tarr[st];
            for (int st1 = st+1; st1 < n; st1++)
            {
                V[st1*n + idx] -= C[st1 + n*st] * dtmp;
            }
        }
    }

    for(int idx = eig_start;idx < eig_stop;idx++)
    {
        for (int st = 0; st < n; st++) G[idx*n + st] = V[st*n + idx];
    }
#endif
}


void gramsch_update_psi(double *V,
                        double *G,
                        double *C,
                        int N,
                        int eig_start,
                        int eig_stop)
{

  dim3 Grid, Block;
  Grid.x = eig_stop - eig_start + 1;
  Grid.y = 1;
  Block.x = N;
  Block.y = 1;

  gramsch_update_psi_kernel<<<Grid, Block, 0>>>(
                                           V,
                                           G,
                                           C,
                                           N,
                                           eig_start,
                                           eig_stop);
}


#endif
