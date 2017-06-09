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

template void FoldedSpectrumScalapackOrtho<double> (int, int, int, int *, int *, double *, double *, double *, double *, double *, Scalapack *);

// Note V and B are dist matrices
template <typename KpointType>
void FoldedSpectrumScalapackOrtho(int n, int eig_start, int eig_stop, int *fs_eigcounts, int *fs_eigstart, KpointType *Vdist, KpointType *V, KpointType *B, KpointType *work1, KpointType *work2, Scalapack *MainSp)
{
    KpointType ZERO_t(0.0);
    KpointType ONE_t(1.0);
    KpointType alpha(1.0);
    KpointType beta(0.0);
    KpointType *Bgpu = NULL; 
    KpointType *Ggpu = NULL; 
    KpointType *Vgpu = NULL; 
    KpointType *Cgpu = NULL; 
#if GPU_ENABLED
    cublasStatus_t custat;
    KpointType *C = (KpointType *)GpuMallocHost(n * n * sizeof(KpointType));
    KpointType *G = (KpointType *)GpuMallocHost(n * n * sizeof(KpointType));
    Bgpu = (KpointType *)GpuMallocDevice(n * n * sizeof(KpointType));
    Ggpu = (KpointType *)GpuMallocDevice(n * n * sizeof(KpointType));
    Vgpu = (KpointType *)GpuMallocDevice(n * n * sizeof(KpointType));
    Cgpu = (KpointType *)GpuMallocDevice(n * n * sizeof(KpointType));
#else
    KpointType *C = work1;
    KpointType *G = work2;
#endif
    double *tarr = new double[n];
    int info = 0, ione = 1;
    char *trans_t="t", *trans_n="n", *cuplo = "l";

    // For mpi routines. Transfer twice as much data for complex orbitals
    int factor = 2;
    if(ct.is_gamma) factor = 1;


    int *m_f_desca = MainSp->GetDistDesca();
    int m_f_dist_length = MainSp->ComputeMdim(n) *  MainSp->ComputeNdim(n);

    // Folded spectrum scalapacks
    //Scalapack *FSp = MainSp->GetNextScalapack();

    static KpointType *m_distC;
    if(!m_distC) {
        int retval1 = MPI_Alloc_mem(m_f_dist_length * sizeof(double) * factor , MPI_INFO_NULL, &m_distC);
        if(retval1 != MPI_SUCCESS)
            rmg_error_handler (__FILE__, __LINE__, "Memory allocation failure in FoldedSpectrumScalapackOrtho");
 
    }

    // Overlaps
    RmgTimer *RT1 = new RmgTimer("4-Diagonalization: fs-Gram-overlaps");
    if(!B) {
#if GPU_ENABLED
        custat = cublasSetVector(n * n , sizeof(KpointType), V, 1, Vgpu, 1 );
        RmgCudaError(__FILE__, __LINE__, custat, "Problem transferreing C from system memory to gpu.");
        cublasDsyrk(ct.cublas_handle, CUBLAS_FILL_MODE_LOWER, CUBLAS_OP_T, n, n, &alpha, Vgpu, n, &beta, Cgpu, n);

#else
//        dsyrk (cuplo, trans_t, &n, &n, &alpha, V, &n, &beta, C, &n);
        pdsyrk_ (cuplo, trans_t, &n, &n, &alpha, Vdist, &ione, &ione, m_f_desca, &beta, m_distC, &ione, &ione, m_f_desca);
#endif
    }
    else {
        // transfer V and B to the GPU for the multiplication and leave the result there
//        RmgGemm(trans_n, trans_n, n, n, n, ONE_t, B, n, V, n, ZERO_t, G, n, Vgpu, Bgpu, Ggpu, true, true, false, false);
        RmgSymm("l", cuplo, n, n, ONE_t, B, n, V, n, ZERO_t, G, n, Bgpu, Vgpu, Ggpu, true, true, false, false);
        // Multiply G by V and leave result in Cgpu for the magma_dpotrf_gpu call coming up next
        RmgGemm(trans_t, trans_n, n, n, n, ONE_t, V, n, G, n, ZERO_t, C, n, Vgpu, Ggpu, Cgpu, false, false, false, false);
    }
    delete(RT1);


    // Cholesky factorization
    RT1 = new RmgTimer("4-Diagonalization: fs-Gram-cholesky");
#if GPU_ENABLED && MAGMA_LIBS
    magma_dpotrf_gpu(MagmaLower, n, Cgpu, n, &info);
    custat = cublasGetVector(n * n, sizeof( KpointType ), Cgpu, 1, C, 1 );
    RmgCudaError(__FILE__, __LINE__, custat, "Problem transferring C matrix from GPU to system memory.");
#elif GPU_ENABLED
    custat = cublasGetVector(n * n, sizeof( KpointType ), Cgpu, 1, C, 1 );
    RmgCudaError(__FILE__, __LINE__, custat, "Problem transferring C matrix from GPU to system memory.");
    dpotrf(cuplo, &n, C, &n, &info);
#else
    //dpotrf(cuplo, &n, C, &n, &info);
    pdpotrf_( cuplo, &n, m_distC, &ione, &ione, m_f_desca, &info );
    for(int i=0;i<n*n;i++)C[i]=0.0;
    MainSp->GatherMatrix(C, m_distC);

#endif
    MainSp->BcastRoot(C, factor * n * n, MPI_DOUBLE);
    delete(RT1);



    RT1 = new RmgTimer("4-Diagonalization: fs-Gram-update");
    // Get inverse of diagonal elements
    for(int ix = 0;ix < n;ix++) tarr[ix] = 1.0 / C[n * ix + ix];

//----------------------------------------------------------------
    for(int idx = 0;idx < n*n;idx++)G[idx] = ZERO_t;

    int idx;
    double *darr;
    int st, st1;
#pragma omp parallel private(idx,st,st1,darr)
{
    darr = new double[n];
#pragma omp barrier

#pragma omp for schedule(static, 1) nowait
    for(idx = eig_start;idx < eig_stop;idx++) {

        for (int st = 0; st < n; st++) darr[st] = V[st*n + idx];

        for (int st = 0; st < n; st++) {

            darr[st] *= tarr[st];

            for (int st1 = st+1; st1 < n; st1++) {
                darr[st1] -= C[st1 + n*st] * darr[st];
            }

        }

        for (st = 0; st < n; st++) G[st*n + idx] = darr[st];

    }

    delete [] darr;

} // end omp section

    delete(RT1);

    // The matrix transpose here lets us use an Allgatherv instead of an Allreduce which
    // greatly reduces the network bandwith required at the cost of doing local transposes.
    RT1 = new RmgTimer("4-Diagonalization: fs-Gram-allreduce");
    for(int i=0;i < n*n;i++)V[i] = ZERO_t;
#if 0
    for(int st1 = eig_start;st1 < eig_stop;st1++) {
        for(int st2 = 0;st2 < n;st2++) {
            V[st1*n + st2] = G[st1 + st2*n];
        }
    }
    MPI_Allreduce(MPI_IN_PLACE, G, n*n * factor, MPI_DOUBLE, MPI_SUM, MainSp->GetComm());
    MainSp->CopySquareMatrixToDistArray(G, Vdist, n, m_f_desca);
#else
    for(int st1 = eig_start;st1 < eig_stop;st1++) {
        for(int st2 = 0;st2 < n;st2++) {
            V[st1*n + st2] = G[st1 + st2*n];
        }
    }

    int eig_step = eig_stop - eig_start;
    MPI_Allgatherv(MPI_IN_PLACE, eig_step * n * factor, MPI_DOUBLE, V, fs_eigcounts, fs_eigstart, MPI_DOUBLE, MainSp->GetComm());
for(int st1 = 0;st1 < n;st1++) {
    for(int st2 = 0;st2 < n;st2++) {
        G[st1*n + st2] = V[st1 + st2*n];
    }
}
    MainSp->CopySquareMatrixToDistArray(G, Vdist, n, m_f_desca);

#endif

#if GPU_ENABLED
    GpuFreeDevice(Cgpu);
    GpuFreeDevice(Vgpu);
    GpuFreeDevice(Ggpu);
    GpuFreeDevice(Bgpu);
#endif
    delete(RT1);

    delete [] tarr;
#if GPU_ENABLED
    GpuFreeHost(G);
    GpuFreeHost(C);
#endif
}


