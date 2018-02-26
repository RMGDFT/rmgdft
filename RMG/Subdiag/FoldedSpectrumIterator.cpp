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
    char *trans_n = "n";


#if GPU_ENABLED
    cublasStatus_t custat;
    int ione = 1;
    int sizr = n * k;
    double *Y = (double *)GpuMallocManaged(n * k *  sizeof(double));
    for(int i = 0;i < n * k;i++) Y[i] = 0.0;
    double *T = (double *)GpuMallocManaged(n * k * sizeof(double));
#else
    double *Y = new double[n * k]();
#endif

    // outer loop over steps
    double neg_rone = -1.0;
    for(int step = 0;step < iterations;step++) {

        // Generate A * X for entire block
        RmgGemm(trans_n, trans_n, n, k, n, ONE_t, A, n, X, n, ZERO_t, Y, n);

        // Subtract off lamda * I component. Gemm call is mainly for simplicity with GPU.
#if GPU_ENABLED
        custat = cublasDdgmm(ct.cublas_handle, CUBLAS_SIDE_RIGHT, n, k, X, n, eigs, ione, T, n);
        RmgCudaError(__FILE__, __LINE__, custat, "Problem executing cublasDdgmm.");
        custat = cublasDaxpy(ct.cublas_handle, sizr, &neg_rone, T, ione, Y, ione);
        RmgCudaError(__FILE__, __LINE__, custat, "Problem executing cublasDaxpy.");
        custat = cublasDaxpy(ct.cublas_handle, sizr, &alpha, Y, ione, X, ione);
        RmgCudaError(__FILE__, __LINE__, custat, "Problem executing cublasDaxpy.");
#else
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
#endif
    }    

#if GPU_ENABLED
    GpuFreeManaged(T);
    GpuFreeManaged(Y);
#else
    delete [] Y;
#endif
}


