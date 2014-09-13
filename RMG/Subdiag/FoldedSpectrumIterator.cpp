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

template void FoldedSpectrumIterator<double> (double *, int, double *, int, double *, double, int);
template <typename DataType>
void FoldedSpectrumIterator(DataType *A, int n, double *eigs, int k, DataType *X, double alpha, int iterations)
{

    DataType ZERO_t(0.0);
    DataType ONE_t(1.0);
    DataType *Y = new DataType[n * k]();
    DataType *C = new DataType[n];
    DataType *NULLptr = NULL;

    char *trans_n = "n";


    // outer loop over steps
    for(int step = 0;step < iterations;step++) {

        // Generate A * X for entire block
        RmgGemm(trans_n, trans_n, n, k, n, ONE_t, A, n, X, n, ZERO_t, Y, n, NULLptr, NULLptr, NULLptr, false, false, false, true);

        // Subtract off lamda * I component
        for(int kcol = 0;kcol < k;kcol++) {
           for(int ix = 0;ix < n;ix++) {
               Y[kcol*n + ix] -= eigs[kcol] * X[kcol*n + ix];
           }
        }

        for(int ix = 0;ix < n*k;ix++) X[ix] += alpha * Y[ix];

    }    

    delete [] C;
    delete [] Y;
}
