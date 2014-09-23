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

#if GPU_ENABLED
    #include <cuda.h>
    #include <cuda_runtime_api.h>
    #include <cublas_v2.h>
    #if MAGMA_LIBS
        #include <magma.h>
    #endif
#endif

// Used to convert the generalized eigenvalue problem A*x = lambda*B*x into standard form
// B^(-1)*A*x = lambda*x
//
// Inputs are A and B matrices. Output is Z = B^(-1)*A. Parallelized over nodes, istart
// and istop delineate the rows of the matrix Z that this nodes is responsible for.
// fs_eigcounts and fs_eigstart are used to collect the data in the MPI_Allgatherv call
// at the end. A and B are both overwritten.
//
template void FoldedSpectrumGSE<double> (double *, double *, double *, int, int, int, int *, int *, int);
template <typename DataType>
void FoldedSpectrumGSE(DataType *A, DataType *B, DataType *Z, int n, int istart, int istop, int *fs_eigcounts, int *fs_eigstart, int iterations)
{
    RmgTimer RT0("Diagonalization: fs: GSE");
    DataType ZERO_t(0.0);
    DataType ONE_t(1.0);
    DataType *D = new DataType[n];

    DataType *NULLptr = NULL;
    int istep = istop - istart;

    // For mpi routines. Transfer twice as much data for complex orbitals
    int factor = 2;
    if(ct.is_gamma) factor = 1;

    char *trans_n = "n";
    char *trans_t = "t";



#if GPU_ENABLED

    RmgTimer *RT1 = new RmgTimer("Diagonalization: fs: GSE-setup");

    int ione = 1;
    cublasStatus_t custat;

    DataType *T1 = NULL;
    DataType *gpuT1 = (DataType *)GpuMalloc(n * n * sizeof(DataType));
    DataType *gpuT2 = (DataType *)GpuMalloc(n * n * sizeof(DataType));
    DataType *gpuA = (DataType *)GpuMalloc(istep * n * sizeof(DataType));
    DataType *gpuZ = (DataType *)GpuMalloc(istep * n * sizeof(DataType));
    DataType *gpuB = (DataType *)GpuMalloc(n * n * sizeof(DataType));
    DataType *gpuD = (DataType *)GpuMalloc(n * sizeof(DataType));

    // Set up D^(-1) and transfer it to the GPU
    for(int ix = 0;ix < n;ix++) D[ix] = 1.0 / B[ix*n + ix];
    custat = cublasSetVector(n , sizeof(DataType), D, ione, gpuD, ione );
    RmgCudaError(__FILE__, __LINE__, custat, "Problem transferreing D from system memory to gpu.");

    // Initial starting guess goes in Z and is just the identity
    //for(int ix = 0;ix < n*n;ix++) Z[ix] = ZERO_t;
    for(int st1 = istart;st1 < istop;st1++){
        for(int st2 = 0;st2 < n;st2++) Z[st1*n + st2] = ZERO_t;
    }
    for(int ix = istart;ix < istop;ix++) Z[ix*n + ix] = ONE_t;

    // T1 starts out as the identity but will hold (I - D-1 * B)
#if MAGMA_LIBS
    magmablas_dlaset(MagmaFull, n, n, ZERO_t, ONE_t, (double *)gpuT1, n);
#else
    T1 = (DataType *)GpuMallocHost(n * n * sizeof(DataType));
    for(int ix = 0;ix < n*n;ix++) T1[ix] = ZERO_t;
    for(int ix = 0;ix < n;ix++) T1[ix*n + ix] = 1.0;
    custat = cublasSetVector(n * n , sizeof(DataType), T1, ione, gpuT1, ione );
    RmgCudaError(__FILE__, __LINE__, custat, "Problem transferreing T1 from system memory to gpu.");
#endif

    // Transfer full B to gpu for Ddgmm call
    custat = cublasSetVector(n * n , sizeof(DataType), B, ione, gpuB, ione );
    RmgCudaError(__FILE__, __LINE__, custat, "Problem transferreing B from system memory to gpu.");

    // D-1 * B stored temporarily in gpuT2
    custat = cublasDdgmm(ct.cublas_handle, CUBLAS_SIDE_RIGHT, n, n, gpuB, n, gpuD, ione, gpuT2, n);

    // (I - D-1 * B) stored in T1
    DataType neg_rone(-1.0);
    custat = cublasDaxpy(ct.cublas_handle, n*n, &neg_rone, gpuT2, ione, gpuT1, ione);
    delete(RT1);


    int st1, st2;
    RT1 = new RmgTimer("Diagonalization: fs: GSE-Second term");
#pragma omp parallel private(st1, st2)
{
#pragma omp for schedule(static, 1) nowait
    // Compute D^(-1) * B * I and store in B
    for(st1 = istart;st1 < istop;st1++){
        for(st2 = 0;st2 < n;st2++){
              B[st1*n + st2] = D[st2] * A[st1*n + st2];
        }
    }
} // end omp parallel region
    delete(RT1);


    RT1 = new RmgTimer("Diagonalization: fs: GSE-First term");

    // Start the slice of A from istart to istop transferring to the GPU
    custat = cublasSetVector(istep * n , sizeof(DataType), &A[istart*n], ione, gpuA, ione );
    RmgCudaError(__FILE__, __LINE__, custat, "Problem transferreing A from system memory to gpu.");

    custat = cublasSetVector(istep * n , sizeof(DataType), &Z[istart*n], ione, gpuZ, ione );
    RmgCudaError(__FILE__, __LINE__, custat, "Problem transferreing T1 from system memory to gpu.");

    custat = cublasSetVector(istep * n , sizeof(DataType), &B[istart*n], ione, gpuB, ione );
    RmgCudaError(__FILE__, __LINE__, custat, "Problem transferreing T1 from system memory to gpu.");

    // outer loop over steps
    for(int step = 0;step < iterations;step++) {

        RmgGemm(trans_n, trans_n, n, istep, n, ONE_t, T1, n, &Z[istart*n], n, ZERO_t, &A[istart*n], n, gpuT1, gpuZ, gpuA, false, false, false, false);

        custat = cublasDgeam(ct.cublas_handle, CUBLAS_OP_N, CUBLAS_OP_N, n, istep, &ONE_t, gpuA, n, &ONE_t, gpuB, n, gpuZ, n);
        RmgCudaError(__FILE__, __LINE__, custat, "Problem executing cublasDgeam.");
    }

    custat = cublasGetVector(istep *n, sizeof( DataType ), gpuZ, 1, &Z[istart*n], 1 );
    RmgCudaError(__FILE__, __LINE__, custat, "Problem transferring Z matrix from GPU to system memory.");

    GpuFree(gpuD);
    GpuFree(gpuB);
    GpuFree(gpuZ);
    GpuFree(gpuA);
    GpuFree(gpuT2);
    GpuFree(gpuT1);

    delete(RT1);
#else

    RmgTimer *RT1 = new RmgTimer("Diagonalization: fs: GSE-setup");

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

    int st1, st2;
#pragma omp parallel private(st1, st2)
{
#pragma omp for schedule(static, 1) nowait
    // (I - D-1 * B)
    for(st1 = 0;st1 < n;st1++){
        for(st2 = 0;st2 < n;st2++){
            T1[st1*n + st2] -= D[st2] * B[st1*n + st2];
        }
    }
}

    delete(RT1);


    RT1 = new RmgTimer("Diagonalization: fs: GSE-Second term");
#pragma omp parallel private(st1, st2)
{
#pragma omp for schedule(static, 1) nowait
    // Compute D^(-1) * B * X and store in B
    for(int st1 = istart;st1 < istop;st1++){
        for(int st2 = 0;st2 < n;st2++){
            B[st1*n + st2] = D[st2] * A[st1*n + st2];
        }
    }
} // end omp parallel region
    delete(RT1);

    RT1 = new RmgTimer("Diagonalization: fs: GSE-First term");

    // outer loop over steps
    for(int step = 0;step < iterations;step++) {

        // Compute (I - D-1 * B) * Z(step) and store in A
        RmgGemm(trans_n, trans_n, n, istep, n, ONE_t, T1, n, &Z[istart*n], n, ZERO_t, &A[istart*n], n, NULLptr, NULLptr, NULLptr, false, false, false, true);

        // Finally generate Z(step+1) = (I - D-1 * B) * Z(step) + D^(-1) * B * X 
        //for(int ix=0;ix < n*n;ix++) Z[ix] = A[ix] + B[ix];
        for(int st1 = istart;st1 < istop;st1++){
            for(int st2 = 0;st2 < n;st2++){
                Z[st1*n + st2] =  A[st1*n + st2] +  B[st1*n + st2];
            }
        }

    }    

    delete(RT1);
#endif


    // Make sure everybody has a copy
    RT1 = new RmgTimer("Diagonalization: fs: GSE-Allgatherv");
    MPI_Allgatherv(MPI_IN_PLACE, istep * n * factor, MPI_DOUBLE, Z, fs_eigcounts, fs_eigstart, MPI_DOUBLE, pct.grid_comm);
    delete(RT1);

#if GPU_ENABLED
#if !MAGMA_LIBS
    GpuFreeHost(T1);
#endif
#else
    delete [] T1;
#endif
    delete [] D;
}
