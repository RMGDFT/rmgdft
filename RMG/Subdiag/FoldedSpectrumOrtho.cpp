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
#include "const.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "rmg_error.h"
#include "RmgTimer.h"
#include "Subdiag.h"
#include "RmgGemm.h"
#include "GpuAlloc.h"
#include "ErrorFuncs.h"
#include "blas.h"


#include "transition.h"

#if GPU_ENABLED
    #include <cuda.h>
    #include <cuda_runtime_api.h>
    #include <cublas_v2.h>
    #if MAGMA_LIBS
        #include <magma.h>
    #endif
#endif

// Gram-Schmidt ortho for eigenvectors.
// n = dimensions of the square matrix
// index of the starting eigenvector this node is responsible for
// stopping eigenvector this node is responsible for
// V = matrix of raw eigenvectors on input, orthonormal eigenvectors on output
// B if not NULL then the orthogonalization condition is  <i|B|j> = I
// else <i|j> = I
//

template void FoldedSpectrumOrtho<double> (int, int, int, int *, int *, double *, double *, int, MPI_Comm &);

template <typename KpointType>
void FoldedSpectrumOrtho(int n, int eig_start, int eig_stop, int *fs_eigcounts, int *fs_eigstart, KpointType *V, KpointType *B, int driver, MPI_Comm &fs_comm)
{
    KpointType ZERO_t(0.0);
    KpointType ONE_t(1.0);

    KpointType *NULLptr = NULL;
    KpointType alpha(1.0);
    KpointType beta(0.0);
#if GPU_ENABLED
    cublasStatus_t custat;
    KpointType *C = (KpointType *)GpuMallocHost(n * n * sizeof(KpointType));
    KpointType *G = (KpointType *)GpuMallocHost(n * n * sizeof(KpointType));
#else
    KpointType *C = new KpointType[n * n];
    KpointType *G = new KpointType[n * n];
#endif
    double *tarr = new double[n];
    int info = 0;
    int eig_step = eig_stop - eig_start;

    char *trans_t="t", *trans_n="n", *cuplo = "l";

    // For mpi routines. Transfer twice as much data for complex orbitals
    int factor = 2;
    if(ct.is_gamma) factor = 1;

    // Overlaps
    RmgTimer *RT1 = new RmgTimer("4-Diagonalization: fs: Gram-overlaps");
    if(!B) {
#if GPU_ENABLED
        cublasXtDsyrk(ct.cublasXt_handle, CUBLAS_FILL_MODE_LOWER, CUBLAS_OP_T, n, n, &alpha, V, n, &beta, C, n);
#else
        dsyrk (cuplo, trans_t, &n, &n, &alpha, V, &n, &beta, C, &n);
#endif
    }
    else {
        // Multiply G by V and leave result in C
       if(n < 128) {

            RmgSymm("l", cuplo, n, n, ONE_t, B, n, V, n, ZERO_t, G, n, NULLptr, NULLptr, NULLptr, true, true, false, false);
            RmgGemm(trans_t, trans_n, n, n, n, ONE_t, V, n, G, n, ZERO_t, C, n, NULLptr, NULLptr, NULLptr, false, false, false, false);

        }
        else {

            // split over PE's if n is large enough
            RmgGemm(trans_t, trans_n, n, eig_step, n, ONE_t, B, n, &V[eig_start*n], n, ZERO_t, &G[eig_start*n], n, 
                    NULLptr, NULLptr, NULLptr, false, false, false, false);

            MPI_Allgatherv(MPI_IN_PLACE, eig_step * n * factor, MPI_DOUBLE, G, fs_eigcounts, fs_eigstart, MPI_DOUBLE, fs_comm);

            RmgGemm(trans_t, trans_n, n, eig_step, n, ONE_t, G, n, &V[eig_start*n], n, ZERO_t, &C[eig_start*n], n, 
                                                NULLptr, NULLptr, NULLptr, false, false, false, false);

            MPI_Allgatherv(MPI_IN_PLACE, eig_step * n * factor, MPI_DOUBLE, C, fs_eigcounts, fs_eigstart, MPI_DOUBLE, fs_comm);

        }

    }
    delete(RT1);

    // Cholesky factorization
    RT1 = new RmgTimer("4-Diagonalization: fs: Gram-cholesky");
#if GPU_ENABLED && MAGMA_LIBS
    magma_dpotrf(MagmaLower, n, C, n, &info);
#else
    dpotrf(cuplo, &n, C, &n, &info);
#endif
    delete(RT1);


    RT1 = new RmgTimer("4-Diagonalization: fs: Gram-update");
    // Get inverse of diagonal elements
    for(int ix = 0;ix < n;ix++) tarr[ix] = 1.0 / C[n * ix + ix];

//----------------------------------------------------------------
    for(int idx = 0;idx < n*n;idx++)G[idx] = ZERO_t;

#pragma omp parallel
{
    double *darr = new double[n];
#pragma omp barrier

#pragma omp for schedule(static, 1) nowait
    for(int idx = eig_start;idx < eig_stop;idx++) {

        for (int st = 0; st < n; st++) darr[st] = V[st*n + idx];

        for (int st = 0; st < n; st++) {

            darr[st] *= tarr[st];

            for (int st1 = st+1; st1 < n; st1++) {
                darr[st1] -= C[st1 + n*st] * darr[st];
            }

        }

        for (int st = 0; st < n; st++) G[st*n + idx] = darr[st];

    }

    delete [] darr;

} // end omp section

    delete(RT1);

    // The matrix transpose here lets us use an Allgatherv instead of an Allreduce which
    // greatly reduces the network bandwith required at the cost of doing local transposes.
    RT1 = new RmgTimer("4-Diagonalization: fs: Gram-allreduce");
//    MPI_Allreduce(MPI_IN_PLACE, G, n*n * factor, MPI_DOUBLE, MPI_SUM, fs_comm);

    for(int st1 = eig_start;st1 < eig_stop;st1++) {
        for(int st2 = 0;st2 < n;st2++) {
            V[st1*n + st2] = G[st1 + st2*n];
        }
    }

    MPI_Allgatherv(MPI_IN_PLACE, eig_step * n * factor, MPI_DOUBLE, V, fs_eigcounts, fs_eigstart, MPI_DOUBLE, fs_comm);


    delete(RT1);

    delete [] tarr;
#if GPU_ENABLED
    GpuFreeHost(G);
    GpuFreeHost(C);
#else
    delete [] G;
    delete [] C;
#endif
}


