/*
 *
 * Copyright 2014 The RMG Project Developers. See the COPYRIGHT file 
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

#include <complex>
#include <omp.h>
#include "FiniteDiff.h"
#include "const.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "rmg_error.h"
#include "RmgTimer.h"
#include "GlobalSums.h"
#include "Kpoint.h"
#include "Subdiag.h"
#include "RmgGemm.h"
#include "GpuAlloc.h"
#include "ErrorFuncs.h"
#include "blas.h"

#include "prototypes.h"
#include "common_prototypes.h"
#include "common_prototypes1.h"
#include "transition.h"

// Performs the iteration step for the folded spectrum method
//
//  X_i(step+1) = X_i(step) + alpha  * (A - lambda_i*I) * X_i(step)
//  
// By doing this over a block of X_i we can optimize memory bandwidth usage
//
// Inputs: A is a square matrix of size n
//         eigs is a vector of length k
//         X is a matrix of size kxn
//         alpha is the multiplicative factor
//         iterations is the number of iterations

void FoldedSpectrumIterator(double *A, int n, double *eigs, int k, double *X, double alpha, int iterations, int driver)
{

    double ONE_t = 1.0;
    double ZERO_t = 0.0;
    double *Agpu = NULL;
    double *Xgpu = NULL;
    double *Ygpu = NULL;
    char *trans_n = "n";
    bool usecuxt;


#if GPU_ENABLED
    int ione = 1;
    int sizr = n * k;
    double *Tgpu = NULL;
    double *eigs_gpu = NULL;
    usecuxt = true;
    if(n <= ct.cublasxt_block_size) usecuxt = false;
    double *Y = (double *)GpuMallocHost(n * k *  sizeof(double));
    for(int i = 0;i < n * k;i++) Y[i] = 0.0;
#else
    usecuxt = false; 
    double *Y = new double[n * k]();
#endif

#if GPU_ENABLED

    cublasStatus_t custat;
    if(!usecuxt) {

        Agpu = (double *)GpuMallocDevice(n * n * sizeof(double));
        Xgpu = (double *)GpuMallocDevice(n * k * sizeof(double));
        Ygpu = (double *)GpuMallocDevice(n * k * sizeof(double));
        Tgpu = (double *)GpuMallocDevice(n * k * sizeof(double));
        eigs_gpu = (double *)GpuMallocDevice(n * sizeof(double));
        custat = cublasSetVector(n * n , sizeof(double), A, ione, Agpu, ione );
        RmgCudaError(__FILE__, __LINE__, custat, "Problem transferreing A from system memory to gpu.");

        // Must be a better way of doing this, setting matrix to zero on CPU and transferring seems wasteful
        custat = cublasSetVector(n * k , sizeof(double), Y, ione, Ygpu, ione );     // Must be a better way of doing this
        RmgCudaError(__FILE__, __LINE__, custat, "Problem transferreing Y from system memory to gpu.");

        // Must be a better way of doing this, setting matrix to zero on CPU and transferring seems wasteful
        custat = cublasSetVector(n * k , sizeof(double), X, ione, Xgpu, ione );     // Must be a better way of doing this
        RmgCudaError(__FILE__, __LINE__, custat, "Problem transferreing X from system memory to gpu.");

        custat = cublasSetVector(n, sizeof(double), eigs, ione, eigs_gpu, ione );     // Must be a better way of doing this
        RmgCudaError(__FILE__, __LINE__, custat, "Problem transferreing Y from system memory to gpu.");

    }

#endif

    // outer loop over steps
    for(int step = 0;step < iterations;step++) {

        // Generate A * X for entire block
        RmgGemm(trans_n, trans_n, n, k, n, ONE_t, A, n, X, n, ZERO_t, Y, n, Agpu, Xgpu, Ygpu, false, false, false, false);

        // Subtract off lamda * I component. Gemm call is mainly for simplicity with GPU.
#if GPU_ENABLED
        if(!usecuxt) {

            double neg_rone = -1.0;
            custat = cublasDdgmm(ct.cublas_handle, CUBLAS_SIDE_RIGHT, n, k, Xgpu, n, eigs_gpu, ione, Tgpu, n);
            RmgCudaError(__FILE__, __LINE__, custat, "Problem executing cublasDdgmm.");
            custat = cublasDaxpy(ct.cublas_handle, sizr, &neg_rone, Tgpu, ione, Ygpu, ione);
            RmgCudaError(__FILE__, __LINE__, custat, "Problem executing cublasDaxpy.");
            custat = cublasDaxpy(ct.cublas_handle, sizr, &alpha, Ygpu, ione, Xgpu, ione);
            RmgCudaError(__FILE__, __LINE__, custat, "Problem executing cublasDaxpy.");

        }
#endif

        if(!GPU_ENABLED || (usecuxt == true)) {

            int kcol, ix;
#pragma omp parallel private(kcol, ix)
{
#pragma omp for schedule(static, 1) nowait
            for(kcol = 0;kcol < k;kcol++) {
               for(ix = 0;ix < n;ix++) {
                   Y[kcol*n + ix] -= eigs[kcol] * X[kcol*n + ix];
               }
            }
}
            //daxpy(&sizr, &alpha, Y, &ione, X, &ione);
            for(int ix = 0;ix < n*k;ix++) X[ix] += alpha * Y[ix];

        }



    }    

#if GPU_ENABLED
    if(!usecuxt) {

        custat = cublasGetVector(n * k, sizeof( double ), Xgpu, 1, X, 1 );
        RmgCudaError(__FILE__, __LINE__, custat, "Problem transferring X matrix from GPU to system memory.");

        GpuFreeDevice(eigs_gpu);
        GpuFreeDevice(Tgpu);
        GpuFreeDevice(Ygpu);
        GpuFreeDevice(Xgpu);
        GpuFreeDevice(Agpu);
    }

    GpuFreeHost(Y);
#else
    delete [] Y;
#endif
}


