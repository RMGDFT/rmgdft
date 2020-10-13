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
#include "Gpufuncs.h"
#include "blas.h"

#include "common_prototypes.h"
#include "common_prototypes1.h"
#include "transition.h"

#if CUDA_ENABLED
    #include <cuda.h>
    #include <cuda_runtime_api.h>
    #include <cublas_v2.h>
#endif

// Used to convert the generalized eigenvalue problem A*x = lambda*B*x into standard form
// B^(-1)*A*x = lambda*x
//
// Inputs are A and B matrices. Output is Z = B^(-1)*A. Parallelized over nodes, istart
// and istop delineate the rows of the matrix Z that this nodes is responsible for.
// fs_eigcounts and fs_eigstart are used to collect the data in the MPI_Allgatherv call
// at the end. A and B are both overwritten.
//
template void FoldedSpectrumGSE<double> (double *, double *, double *, int, int, int, int *, int *, int, int, MPI_Comm &);
template <typename DataType>
void FoldedSpectrumGSE(DataType *A, DataType *B, DataType *Z, int n, int istart, int istop, int *fs_eigcounts, int *fs_eigstart, int iterations, int driver, MPI_Comm &fs_comm)
{
    RmgTimer RT0("4-Diagonalization: fs: GSE");
    DataType ZERO_t(0.0);
    DataType ONE_t(1.0);

    int istep = istop - istart;

    // For mpi routines. Transfer twice as much data for complex orbitals
    int factor = 2;
    if(ct.is_gamma) factor = 1;

    char *trans_n = "n";



#if CUDA_ENABLED

    cublasStatus_t custat;
    RmgTimer *RT1 = new RmgTimer("4-Diagonalization: fs: GSE-setup");
    DataType *D = (DataType *)GpuMallocManaged(n * sizeof(DataType));
    DataType *T1 = (DataType *)GpuMallocManaged(n * n * sizeof(DataType));
    double *unitvector = (double *)GpuMallocManaged(n*sizeof(double));
    double *tvector = (double *)GpuMallocManaged(n*sizeof(double)); 
    cublasDcopy(ct.cublas_handle, n, B, n+1, tvector, 1);
    DeviceSynchronize();
    for(int ix = 0;ix < n;ix++) D[ix] = 1.0 / tvector[ix]; 
    DeviceSynchronize();


    // Set up D^(-1) and transfer it to the GPU
    //for(int ix = 0;ix < n;ix++) D[ix] = 1.0 / B[ix*n + ix];

    // Initial starting guess goes in Z and is just the identity
    GpuFill((double *)&Z[istart*n], istep*n, 0.0);
    GpuFill((double *)unitvector, n, 1.0);
    cublasDcopy(ct.cublas_handle, istep, unitvector, 1, (double *)&Z[istart*n], istep+1);

    //for(int st1 = istart;st1 < istop;st1++){
    //    for(int st2 = 0;st2 < n;st2++) Z[st1*n + st2] = ZERO_t;
    //}
    //for(int ix = istart;ix < istop;ix++) Z[ix*n + ix] = ONE_t;
    

    // T1 starts out as the identity but will hold (I - D-1 * B)

    // (I - D-1 * B)
    double negrone = -1.0;
    double rone = 1.0;
    custat = cublasDdgmm(ct.cublas_handle, CUBLAS_SIDE_LEFT, n, n, B, n, D, 1, T1, n);
    custat = cublasDscal(ct.cublas_handle, n*n, &negrone, T1, 1);
    custat = cublasDaxpy(ct.cublas_handle, n, &rone, unitvector, 1, T1, n+1);

//#pragma omp for schedule(static, 1) nowait
//    for(int st1 = 0;st1 < n;st1++){
//        for(int st2 = 0;st2 < n;st2++){
//            T1[st1*n + st2] = -D[st2] * B[st1*n + st2];
//        }
//    }
//    for(int ix = 0;ix < n;ix++) T1[ix*n + ix] += 1.0;

    GpuFreeManaged(unitvector);
    delete(RT1);

    RT1 = new RmgTimer("4-Diagonalization: fs: GSE-Second term");
    // Compute D^(-1) * B * I and store in B
    custat = cublasDdgmm(ct.cublas_handle, CUBLAS_SIDE_LEFT, n, istep, &A[istart*n], n, D, 1, &B[istart*n], n);
//#pragma omp for schedule(static, 1) nowait
//    for(int st1 = istart;st1 < istop;st1++){
//        for(int st2 = 0;st2 < n;st2++){
//              B[st1*n + st2] = D[st2] * A[st1*n + st2];
//        }
//    }


    delete(RT1);
    RT1 = new RmgTimer("4-Diagonalization: fs: GSE-First term");


    // outer loop over steps
    int device = -1;
    cudaGetDevice(&device);
    DeviceSynchronize();
    for(int step = 0;step < iterations;step++) {

            RmgGemm(trans_n, trans_n, n, istep, n, ONE_t, T1, n, &Z[istart*n], n, ZERO_t, &A[istart*n], n);
            // Finally generate Z(step+1) = (I - D-1 * B) * Z(step) + D^(-1) * B * X 
            //for(int ix=0;ix < n*n;ix++) Z[ix] = A[ix] + B[ix];
            custat = cublasDgeam(ct.cublas_handle, CUBLAS_OP_N, CUBLAS_OP_N, n, istep, &ONE_t, 
                                 &A[istart*n], n, &ONE_t, &B[istart*n], n, &Z[istart*n], n);
            RmgGpuError(__FILE__, __LINE__, custat, "Problem executing cublasDgeam.");
//#pragma omp for schedule(static, 1) nowait
//            for(int st1 = istart;st1 < istop;st1++){
//                for(int st2 = 0;st2 < n;st2++){
//                    Z[st1*n + st2] =  A[st1*n + st2] +  B[st1*n + st2];
//                }
//            }


    }
    DeviceSynchronize();
    gpuMemPrefetchAsync ( Z, n*n*sizeof(double), cudaCpuDeviceId, NULL);
    GpuFreeManaged(T1);
    GpuFreeManaged(D);

    delete(RT1);
#else

    DataType *D = new DataType[n];
    RmgTimer *RT1 = new RmgTimer("4-Diagonalization: fs: GSE-setup");

    DataType *T1 = new DataType[n*n]();

    // Set up D^(-1)
    for(int ix = 0;ix < n;ix++) D[ix] = 1.0 / B[ix*n + ix];

    // Initial starting guess is just the identity
    //for(int ix = 0;ix < n*n;ix++) Z[ix] = ZERO_t;
    for(int st1 = istart;st1 < istop;st1++){
        for(int st2 = 0;st2 < n;st2++) Z[st1*n + st2] = ZERO_t;
    }
    for(int ix = istart;ix < istop;ix++) Z[ix*n + ix] = ONE_t;

    // T1 starts out as the identity but will hold (I - D-1 * B)
    for(int ix = 0;ix < n;ix++) T1[ix*n + ix] = 1.0;

#pragma omp for schedule(static, 1) nowait
    // (I - D-1 * B)
    for(int st1 = 0;st1 < n;st1++){
        for(int st2 = 0;st2 < n;st2++){
            T1[st1*n + st2] -= D[st2] * B[st1*n + st2];
        }
    }

    delete(RT1);


    RT1 = new RmgTimer("4-Diagonalization: fs: GSE-Second term");
#pragma omp for schedule(static, 1) nowait
    // Compute D^(-1) * B * X and store in B
    for(int st1 = istart;st1 < istop;st1++){
        for(int st2 = 0;st2 < n;st2++){
            B[st1*n + st2] = D[st2] * A[st1*n + st2];
        }
    }
    delete(RT1);

    RT1 = new RmgTimer("4-Diagonalization: fs: GSE-First term");

    // outer loop over steps
    for(int step = 0;step < iterations;step++) {

        // Compute (I - D-1 * B) * Z(step) and store in A
        RmgGemm(trans_n, trans_n, n, istep, n, ONE_t, T1, n, &Z[istart*n], n, ZERO_t, &A[istart*n], n);

        // Finally generate Z(step+1) = (I - D-1 * B) * Z(step) + D^(-1) * B * X 
        //for(int ix=0;ix < n*n;ix++) Z[ix] = A[ix] + B[ix];
#pragma omp for schedule(static, 1) nowait
        for(int st1 = istart;st1 < istop;st1++){
            for(int st2 = 0;st2 < n;st2++){
                Z[st1*n + st2] =  A[st1*n + st2] +  B[st1*n + st2];
            }
        }

    }    

    delete(RT1);
    delete [] T1;
    delete [] D;
#endif

    // Make sure everybody has a copy
    RT1 = new RmgTimer("4-Diagonalization: fs: GSE-Allgatherv");
    MPI_Allgatherv(MPI_IN_PLACE, istep * n * factor, MPI_DOUBLE, Z, fs_eigcounts, fs_eigstart, MPI_DOUBLE, fs_comm);
    delete(RT1);

}
