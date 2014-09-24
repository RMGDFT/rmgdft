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

template void FoldedSpectrumOrtho<double> (int, int, int, int *, int *, double *, double *);

template <typename KpointType>
void FoldedSpectrumOrtho(int n, int eig_start, int eig_stop, int *fs_eigcounts, int *fs_eigstart, KpointType *V, KpointType *B)
{
    RmgTimer RT0("Diagonalization: fs: Gram-Schmidt");
    KpointType ZERO_t(0.0);
    KpointType ONE_t(1.0);
    KpointType *NULLptr = NULL;
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
    Bgpu = (KpointType *)GpuMalloc(n * n * sizeof(KpointType));
    Ggpu = (KpointType *)GpuMalloc(n * n * sizeof(KpointType));
    Vgpu = (KpointType *)GpuMalloc(n * n * sizeof(KpointType));
    Cgpu = (KpointType *)GpuMalloc(n * n * sizeof(KpointType));
#else
    KpointType *C = new KpointType[n * n];
    KpointType *G = new KpointType[n * n];
#endif
    double *tarr = new double[n];
    int info = 0;

    char *trans_t="t", *trans_n="n", *cuplo = "l";

    // For mpi routines. Transfer twice as much data for complex orbitals
    int factor = 2;
    if(ct.is_gamma) factor = 1;

    // Overlaps
    RmgTimer *RT1 = new RmgTimer("Diagonalization: fs: overlaps");
    if(!B) {
#if GPU_ENABLED
        custat = cublasSetVector(n * n , sizeof(KpointType), V, 1, Vgpu, 1 );
        RmgCudaError(__FILE__, __LINE__, custat, "Problem transferreing C from system memory to gpu.");
        cublasDsyrk(ct.cublas_handle, CUBLAS_FILL_MODE_LOWER, CUBLAS_OP_T, n, n, &alpha, Vgpu, n, &beta, Cgpu, n);

#else
        dsyrk (cuplo, trans_t, &n, &n, &alpha, V, &n, &beta, C, &n);
#endif
    }
    else {
        // transfer V and B to the GPU for the multiplication and leave the result there
        RmgGemm(trans_t, trans_n, n, n, n, ONE_t, V, n, B, n, ZERO_t, G, n, Vgpu, Bgpu, Ggpu, true, true, false, false);
        // Multiply G by V and leave result in Cgpu for the magma_dpotrf_gpu call coming up next
        RmgGemm(trans_n, trans_n, n, n, n, ONE_t, G, n, V, n, ZERO_t, C, n, Ggpu, Vgpu, Cgpu, false, false, false, false);
    }
    delete(RT1);


    // Cholesky factorization
    RT1 = new RmgTimer("Diagonalization: fs: Gram-Schmidt cholesky");
#if GPU_ENABLED && MAGMA_LIBS
    magma_dpotrf_gpu(MagmaLower, n, Cgpu, n, &info);
    custat = cublasGetVector(n * n, sizeof( KpointType ), Cgpu, 1, C, 1 );
    RmgCudaError(__FILE__, __LINE__, custat, "Problem transferring C matrix from GPU to system memory.");
#elif GPU_ENABLED
    custat = cublasGetVector(n * n, sizeof( KpointType ), Cgpu, 1, C, 1 );
    RmgCudaError(__FILE__, __LINE__, custat, "Problem transferring C matrix from GPU to system memory.");
    dpotrf(cuplo, &n, C, &n, &info);
#else
    dpotrf(cuplo, &n, C, &n, &info);
#endif
    delete(RT1);


    RT1 = new RmgTimer("Diagonalization: fs: Gram-Schmidt update");
    // Get inverse of diagonal elements
    for(int ix = 0;ix < n;ix++) tarr[ix] = 1.0 / C[n * ix + ix];

//----------------------------------------------------------------
    for(int idx = 0;idx < n*n;idx++)G[idx] = ZERO_t;

    int idx, omp_tid;
    double *darr, *sarr;
    int st, st1;
#pragma omp parallel private(idx,st,st1,omp_tid,sarr)
{
    omp_tid = omp_get_thread_num();
    if(omp_tid == 0) darr = new double[n * omp_get_num_threads()];
#pragma omp barrier

#pragma omp for schedule(static, 1) nowait
    for(idx = eig_start;idx < eig_stop;idx++) {

        sarr = &darr[omp_tid*n];

        for (int st = 0; st < n; st++) sarr[st] = V[st*n + idx];

        for (int st = 0; st < n; st++) {

            sarr[st] *= tarr[st];

            for (int st1 = st+1; st1 < n; st1++) {
                sarr[st1] -= C[st1 + n*st] * sarr[st];
            }

        }

        for (st = 0; st < n; st++) G[st*n + idx] = sarr[st];

    }
} // end omp section

    delete [] darr;
    delete(RT1);

    // The matrix transpose here lets us use an Allgatherv instead of an Allreduce which
    // greatly reduces the network bandwith required at the cost of doing local transposes.
    RT1 = new RmgTimer("Diagonalization: fs: Gram-Schmidt allreduce3");
//    MPI_Allreduce(MPI_IN_PLACE, G, n*n * factor, MPI_DOUBLE, MPI_SUM, pct.grid_comm);

    for(int st1 = eig_start;st1 < eig_stop;st1++) {
        for(int st2 = 0;st2 < n;st2++) {
            V[st1*n + st2] = G[st1 + st2*n];
        }
    }

    int eig_step = eig_stop - eig_start;
    MPI_Allgatherv(MPI_IN_PLACE, eig_step * n * factor, MPI_DOUBLE, V, fs_eigcounts, fs_eigstart, MPI_DOUBLE, pct.grid_comm);

#if 0
    // Transpose the full matrix backwards. Maybe use OMP for this on the CPU?
#if GPU_ENABLED
    int ione = 1;
    custat = cublasSetVector(n * n , sizeof(KpointType), V, ione, Vgpu, ione );
    RmgCudaError(__FILE__, __LINE__, custat, "Problem transferring V matrix from system memory to GPU.");

    custat = cublasDgeam(ct.cublas_handle,  CUBLAS_OP_T,  CUBLAS_OP_N, n, n, &ONE_t, Vgpu, n, &ZERO_t, Bgpu, n, Ggpu, n);
    RmgCudaError(__FILE__, __LINE__, custat, "Matrix transpose with cublasDgeam failed.");

    custat = cublasGetVector(n * n, sizeof( KpointType ), Ggpu, 1, V, 1 );
    RmgCudaError(__FILE__, __LINE__, custat, "Problem transferring G matrix from GPU to system memory.");

#else
    for(int st1 = 0;st1 < n;st1++) {
        for(int st2 = 0;st2 < n;st2++) {
            G[st1*n + st2] = V[st1 + st2*n];
        }
    }
    for(int idx = 0;idx < n*n;idx++) V[idx] = G[idx];

#endif
#endif


#if GPU_ENABLED
    GpuFree(Cgpu);
    GpuFree(Vgpu);
    GpuFree(Ggpu);
    GpuFree(Bgpu);
#endif
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


