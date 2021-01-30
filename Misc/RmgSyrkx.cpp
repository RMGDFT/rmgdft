#include <complex>
#include <typeinfo>
#include <string.h>

#include "const.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "RmgGemm.h"
#include "GpuAlloc.h"
#include "ErrorFuncs.h"
#include "blas.h"
#include "transition.h"


#if CUDA_ENABLED
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cublas_v2.h>
#endif


/*
  These functions are used to hide the details of the matrix multiplication data types and GPU 
  utilization from the higher level routines.

  The first 13 arguments are the same as the standard dgemm args but with scalar quantities passed
  by value instead of by reference.

*/


template void RmgSyrkx<double>(char *, char *, int, int, double, double *, int, double *, int, 
                                  double, double *, int);

template void RmgSyrkx<std::complex<double> >(char *, char *, int, int, std::complex<double>, 
                      std::complex<double> *, int, std::complex<double> *, int, 
                      std::complex<double>, std::complex<double> *, int);


template <typename DataType> void RmgSyrkx(char *uplo, char *trans, int n, int k, 
                             DataType alpha, DataType *A, int lda, DataType *B, int ldb, DataType beta, 
                             DataType *C, int ldc)
{

#if BLAS_PROFILE
    if(typeid(DataType) == typeid(std::complex<double>))
    {
        if(pct.gridpe==0) printf("ZSYRK CALL n=%d k=%d\n",n,k);
    }
    else
    {
        if(pct.gridpe==0) printf("DSYRK CALL n=%d k=%d\n",n,k);
    }
#endif

#if CUDA_ENABLED

    cublasStatus_t custat;
    cublasOperation_t cu_trans = CUBLAS_OP_N;
    cublasFillMode_t fill_mode = CUBLAS_FILL_MODE_LOWER;

    if(!strcmp(uplo, "u")) fill_mode = CUBLAS_FILL_MODE_UPPER;
    if(!strcmp(uplo, "U")) fill_mode = CUBLAS_FILL_MODE_UPPER;

    if(!strcmp(trans, "t")) cu_trans = CUBLAS_OP_T;
    if(!strcmp(trans, "T")) cu_trans = CUBLAS_OP_T;
    if(!strcmp(trans, "c")) cu_trans = CUBLAS_OP_C;
    if(!strcmp(trans, "C")) cu_trans = CUBLAS_OP_C;

    DeviceSynchronize();
    if(typeid(DataType) == typeid(std::complex<double>)) {
        custat = cublasZsyrkx(ct.cublas_handle, fill_mode, cu_trans, n, k,
                            (cuDoubleComplex *)&alpha,
                            (cuDoubleComplex*)A, lda,
                            (cuDoubleComplex*)B, ldb,
                            (cuDoubleComplex*)&beta, (cuDoubleComplex*)C, ldc );
        ProcessGpublasError(custat);
        RmgGpuError(__FILE__, __LINE__, custat, "Problem executing cublasZsyrkx");
    }
    else {
        custat = cublasDsyrkx(ct.cublas_handle, fill_mode, cu_trans, n, k,
                            (double*)&alpha,
                            (double*)A, lda,
                            (double*)B, ldb,
                            (double*)&beta, (double*)C, ldc );
        ProcessGpublasError(custat);
        RmgGpuError(__FILE__, __LINE__, custat, "Problem executing cublasDsyrkx");
    }
    DeviceSynchronize();
    return;

#elif HIP_ENABLED
RmgGemm (trans, "N", n, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
return;
#if 0
    hipblasStatus_t hipstat;
    hipblasOperation_t hip_trans = HIPBLAS_OP_N;
    hipblasFillMode_t fill_mode = HIPBLAS_FILL_MODE_LOWER;

    if(!strcmp(uplo, "u")) fill_mode = HIPBLAS_FILL_MODE_UPPER;
    if(!strcmp(uplo, "U")) fill_mode = HIPBLAS_FILL_MODE_UPPER;

    if(!strcmp(trans, "t")) hip_trans = HIPBLAS_OP_T;
    if(!strcmp(trans, "T")) hip_trans = HIPBLAS_OP_T;
    if(!strcmp(trans, "c")) hip_trans = HIPBLAS_OP_C;
    if(!strcmp(trans, "C")) hip_trans = HIPBLAS_OP_C;

    int ka = n;
    if(!strcmp("n", trans)) ka = k;
    if(!strcmp("N", trans)) ka = k;

    int kb = n;
    if(!strcmp("n", trans)) kb = k;
    if(!strcmp("N", trans)) kb = k;

    size_t a_size = (size_t)lda * (size_t)ka;
    size_t b_size = (size_t)ldb * (size_t)kb;
    size_t c_size = (size_t)ldc * (size_t)n;

    if(typeid(DataType) == typeid(std::complex<double>)) {
        std::complex<double> *dA, *dB, *dC;
        gpuMalloc((void **)&dA, a_size * sizeof(std::complex<double>));
        gpuMalloc((void **)&dB, b_size * sizeof(std::complex<double>));
        gpuMalloc((void **)&dC, c_size * sizeof(std::complex<double>));
        hipMemcpyHtoD(dA, A, a_size * sizeof(std::complex<double>));
        hipMemcpyHtoD(dB, B, b_size * sizeof(std::complex<double>));
        if(std::abs(beta) != 0.0) hipMemcpyHtoD(dC, C, c_size * sizeof(std::complex<double>));
        hipstat = hipblasZsyrkx(ct.hipblas_handle, fill_mode, hip_trans, n, k,
                            (hipblasDoubleComplex *)&alpha,
                            (hipblasDoubleComplex*)dA, lda,
                            (hipblasDoubleComplex*)dB, ldb,
                            (hipblasDoubleComplex*)&beta, (hipblasDoubleComplex*)dC, ldc );
        hipMemcpyDtoH(dC, C, c_size * sizeof(std::complex<double>));
        gpuFree(dC);
        gpuFree(dB);
        gpuFree(dA);
        ProcessGpublasError(hipstat);
        RmgGpuError(__FILE__, __LINE__, hipstat, "Problem executing cublasZsyrkx");
    }
    else {
        double *dA, *dB, *dC;
        gpuMalloc((void **)&dA, a_size * sizeof(double));
        gpuMalloc((void **)&dB, b_size * sizeof(double));
        gpuMalloc((void **)&dC, c_size * sizeof(double));
        hipMemcpyHtoD(dA, A, a_size * sizeof(double));
        hipMemcpyHtoD(dB, B, b_size * sizeof(double));
        if(beta != 0.0) hipMemcpyHtoD(dC, C, c_size * sizeof(double));
        hipstat = hipblasDsyrkx(ct.hipblas_handle, fill_mode, hip_trans, n, k,
                            (double*)&alpha,
                            (double*)dA, lda,
                            (double*)dB, ldb,
                            (double*)&beta, (double*)dC, ldc );
        ProcessGpublasError(hipstat);
        hipMemcpyDtoH(C, dC, c_size * sizeof(double));
        gpuFree(dC);
        gpuFree(dB);
        gpuFree(dA);
        RmgGpuError(__FILE__, __LINE__, hipstat, "Problem executing hipblasDsyrkx");

    }
#endif

#else

    // No standard CPU version of syrkx so just use gemm call
    RmgGemm (trans, "N", n, n, k, alpha, A, lda, B, ldb, beta, C, ldc);

#endif
}

