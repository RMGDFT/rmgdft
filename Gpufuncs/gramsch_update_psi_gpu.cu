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
#include "blas.h"
#include "ErrorFuncs.h"

// In some cases this may be faster than the version using cublasDger inside a loop
// but it only uses a single SM per eigenstate so 
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
        __syncthreads();
        for(int st1 = blockIdx.x * blockDim.x + threadIdx.x + st + 1;st1 < N;st1 += blockDim.x * gridDim.x) 
        {
            darr[st1] -= C[st1 + N*st] * darr[st];
        }
        __syncthreads();

    }

}


void gramsch_update_psi(double *V,
                        double *C,
                        int N,
                        int eig_start,
                        int eig_stop,
                        cublasHandle_t cublasH)
{
    int eig_step = eig_stop - eig_start;
    int ione = 1;
    double alpha = -1.0;

    // We get the inverse of the diagonal elements here rather than inside the loop to avoid page faults
    double *darr;
    RmgCudaError(__FILE__, __LINE__, cudaMallocManaged ( &darr, N*sizeof(double), cudaMemAttachGlobal ), "Error: cudaMallocManaged failed.\n");
    //for(int i = 0;i < N;i++) darr[i] = 1.0 / C[i*N + i];
    cublasDcopy(cublasH, N, C, N + 1, darr, 1);
    cudaDeviceSynchronize();
    for(int i = 0;i < N;i++) darr[i] = 1.0 / darr[i];
    cudaDeviceSynchronize();
    /* apply inverse of cholesky factor to states */
    for (int st = 0; st < N; st++)
    {

        /* normalize V[st] */
        cublasDscal(cublasH, eig_step, &darr[st], &V[st * N + eig_start], ione);

        /* subtract the projection along c[st] from the remaining vectors */
        int idx = N - st - 1;
        if(idx)
        {
            cublasDger(cublasH, eig_step, idx, &alpha, &V[st * N + eig_start], ione,
               &C[(st+1) + N*st], ione, &V[(st+1) * N + eig_start], N);
        }

    } /* end of for */

    cudaFree(darr);

}


#endif
