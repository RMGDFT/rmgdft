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


template void FoldedSpectrumGSE<double> (double *, double *, double *, int, int, int, int *, int *, int);
template <typename DataType>
void FoldedSpectrumGSE(DataType *A, DataType *B, DataType *Z, int n, int istart, int istop, int *fs_eigcounts, int *fs_eigstart, int iterations)
{
    RmgTimer RT0("Diagonalization: fs: GSE");
    DataType ZERO_t(0.0);
    DataType ONE_t(1.0);
    DataType *D = new DataType[n];
    DataType *X = new DataType[n*n]();
    DataType *T1 = new DataType[n*n]();
    DataType *T2 = new DataType[n*n];
    DataType *NULLptr = NULL;
    int istep = istop - istart;

    // For mpi routines. Transfer twice as much data for complex orbitals
    int factor = 2;
    if(ct.is_gamma) factor = 1;

    char *trans_n = "n";
    char *trans_t = "n";

    RmgTimer *RT1 = new RmgTimer("Diagonalization: fs: GSE-setup");
    // Set up D^(-1)
    for(int ix = 0;ix < n;ix++) D[ix] = 1.0 / B[ix*n + ix];

    // Initial starting guess is just the identity
    for(int ix = 0;ix < n*n;ix++) Z[ix] = ZERO_t;
    for(int ix = istart;ix < istop;ix++) Z[ix*n + ix] = 1.0;

    // T1 starts out as the identity but will hold (I - D-1 * B)
    for(int ix = 0;ix < n;ix++) T1[ix*n + ix] = 1.0;

    // Create unitary matrix X
    for(int ix = 0;ix < n;ix++) X[ix*n + ix] = 1.0;

    // (I - D-1 * B)
    for(int st1 = 0;st1 < n;st1++){
        for(int st2 = 0;st2 < n;st2++){
            T1[st1*n + st2] -= D[st2] * B[st1*n + st2];
        }
    }

    // D-1 * A
    for(int st1 = 0;st1 < n;st1++){
        for(int st2 = 0;st2 < n;st2++){
            T2[st1*n + st2] = D[st2] * A[st1*n + st2];
        }
    }

    delete(RT1);

    RT1 = new RmgTimer("Diagonalization: fs: GSE-Second term");
    // Compute D^(-1) * B * X and store in B
    //RmgGemm(trans_n, trans_t, n, n, n, ONE_t, T2, n, X, n, ZERO_t, B, n, NULLptr, NULLptr, NULLptr, false, false, false, true);
    for(int st1 = istart;st1 < istop;st1++){
        for(int st2 = 0;st2 < n;st2++){
            B[st1*n + st2] = T2[st1*n + st2] * X[st1*n + st1];
        }
    }
    delete(RT1);

    // outer loop over steps
    for(int step = 0;step < iterations;step++) {

        RT1 = new RmgTimer("Diagonalization: fs: GSE-First term");
        // Compute (I - D-1 * B) * Z(step) and store in A
        RmgGemm(trans_n, trans_t, n, istep, n, ONE_t, T1, n, &Z[istart*n], n, ZERO_t, &A[istart*n], n, NULLptr, NULLptr, NULLptr, false, false, false, true);
        delete(RT1);

        RT1 = new RmgTimer("Diagonalization: fs: GSE-Sum");
        // Finally generate Z(step+1) = (I - D-1 * B) * Z(step) + D^(-1) * B * X 
        //for(int ix=0;ix < n*n;ix++) Z[ix] = A[ix] + B[ix];
        for(int st1 = istart;st1 < istop;st1++){
            for(int st2 = 0;st2 < n;st2++){
                Z[st1*n + st2] =  A[st1*n + st2] +  B[st1*n + st2];
            }
        }
        delete(RT1);

    }    

    // Make sure everybody has a copy
    RT1 = new RmgTimer("Diagonalization: fs: GSE-Allgatherv");
    MPI_Allgatherv(MPI_IN_PLACE, istep * n * factor, MPI_DOUBLE, Z, fs_eigcounts, fs_eigstart, MPI_DOUBLE, pct.grid_comm);
    delete(RT1);

    delete [] T2;
    delete [] T1;
    delete [] X;
    delete [] D;
}
