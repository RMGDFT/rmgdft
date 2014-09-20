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
    KpointType *C = new KpointType[n * n];
    KpointType *G = new KpointType[n * n];
    double *tarr = new double[n];
    int info = 0;

    char *trans_t="t", *trans_n="n", *cuplo = "l";

    // For mpi routines. Transfer twice as much data for complex orbitals
    int factor = 2;
    if(ct.is_gamma) factor = 1;

    // Overlaps
    RmgTimer *RT1 = new RmgTimer("Diagonalization: fs: overlaps");
    if(!B) {
        dsyrk (cuplo, trans_t, &n, &n, &alpha, V, &n, &beta, C, &n);
    }
    else {
        RmgGemm(trans_t, trans_n, n, n, n, ONE_t, V, n, B, n, ZERO_t, G, n, NULLptr, NULLptr, NULLptr, false, false, false, true);
        RmgGemm(trans_n, trans_n, n, n, n, ONE_t, G, n, V, n, ZERO_t, C, n, NULLptr, NULLptr, NULLptr, false, false, false, true);
    }
    delete(RT1);


    // Cholesky factorization
    RT1 = new RmgTimer("Diagonalization: fs: Gram-Schmidt cholesky");
    dpotrf(cuplo, &n, C, &n, &info);
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

    // Transpose the full matrix backwards. Maybe use OMP for this?
    for(int st1 = 0;st1 < n;st1++) {
        for(int st2 = 0;st2 < n;st2++) {
            G[st1*n + st2] = V[st1 + st2*n];
        }
    }

    delete(RT1);
    for(int idx = 0;idx < n*n;idx++) V[idx] = G[idx];

    delete [] tarr;
    delete [] G;
    delete [] C;
}
