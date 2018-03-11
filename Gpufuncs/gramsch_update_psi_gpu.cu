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
#include <cublas_v2.h>
#include <cublasXt.h>
#include "ErrorFuncs.h"

__global__ void gramsch_update_psi_kernel(
                                     const double * __restrict__ C,
                                     int N,
                                     double *darr,
                                     int eig)
{

    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    for(int st=0;st < N;st++)
    {

        if(tid == 0) darr[st] = darr[st] / C[st*N + st];
        __threadfence();
        for(int st1 = blockIdx.x * blockDim.x + threadIdx.x + st + 1;st1 < N;st1 += blockDim.x * gridDim.x) 
        {
            darr[st1] -= C[st1 + N*st] * darr[st];
        }
        __threadfence();

    }

}


void gramsch_update_psi(double *V,
                        double *G,
                        double *C,
                        int N,
                        int eig_start,
                        int eig_stop,
                        cublasHandle_t cublasH)
{

  int threads = 32;
  int blocks = N / threads + 1;
  double *darr;
  RmgCudaError(__FILE__, __LINE__, cudaMallocManaged ( &darr, N*sizeof(double), cudaMemAttachGlobal ), "Error: cudaMallocManaged failed.\n");

  for(int eig = eig_start;eig < eig_stop;eig++)
  {
      cublasDcopy(cublasH, N, &V[eig], N, darr, 1);
      gramsch_update_psi_kernel<<<blocks, threads>>>(C, N, darr, eig);
      cublasDcopy(cublasH, N, darr, 1, &G[eig*N], 1);
  }

  cudaFree(darr);

}


#endif
