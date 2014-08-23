#include <complex>
#include "const.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "rmg_alloc.h"
#include "rmg_error.h"
#include "Subdiag.h"


#if GPU_ENABLED
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cublas_v2.h>
#endif

#define         dgemm   dgemm_
#define         zgemm   zgemm_


/*
  These functions are used to simplify the subspace diagonalization code by hiding the
  details of data type and GPU utilization from the higher level routines.

*/

extern "C" {
void dgemm(const char *, const char *, int *, int *, int *, double *, double *, int *, double *, int *, double *, double *, int *);
void zgemm(const char *, const char *, int *, int *, int *, std::complex<double> *, std::complex<double> *, int *, std::complex<double> *, int *, std::complex<double> *, std::complex<double> *, int *);
}

template void SubdiagGemm<double>(char *, char *, int, int, int, double, double *, int, double *, int, double, double *, int);
template void SubdiagGemm<std::complex<double> >(char *, char *, int, int, int, std::complex<double>, std::complex<double> *, int, std::complex<double> *, int, std::complex<double>, std::complex<double> *, int);

template <typename DataType> void SubdiagGemm(char *transa, char *transb, int m, int n, int k, DataType alpha, DataType *A, int lda, DataType *B, int ldb, DataType beta, DataType *C, int ldc)
{

    if(typeid(DataType) == typeid(std::complex<double>)) {
        zgemm(transa, transb, &m, &n, &k, (std::complex<double> *)&alpha, (std::complex<double> *)A, &lda, (std::complex<double> *)B, &ldb, (std::complex<double> *)&beta, (std::complex<double> *)C, &ldc);
    }
    else {
        dgemm(transa, transb, &m, &n, &k, (double *)&alpha, (double *)A, &lda, (double *)B, &ldb, (double *)&beta, (double *)C, &ldc);
    }
}
