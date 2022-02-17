#include <complex>
#include <typeinfo>
#include <string.h>

#include "const.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "RmgGemm.h"
#include "GpuAlloc.h"
#include "ErrorFuncs.h"
#include "transition.h"



#if CUDA_ENABLED
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cublas_v2.h>
#endif

#define         dsymm   dsymm_
#define         zsymm   zsymm_


extern "C" {
void dsymm(const char *, const char *, int *, int *, double *, double *, int *, double *, int *, double *, double *, int *);
void zsymm(const char *, const char *, int *, int *, std::complex<double> *, std::complex<double> *, int *, std::complex<double> *, int *, std::complex<double> *, std::complex<double> *, int *);
}




/*
  These functions are used to hide the details of the matrix multiplication data types and GPU 
  utilization from the higher level routines.

  The first 13 arguments are the same as the standard dgemm args but with scalar quantities passed
  by value instead of by reference. The last three arguments are used only when CUDA_ENABLED is true.
  In that case if

  [ABC]gpu == NULL   transfer [ABC] to gpu and perform matrix multiplication there. 

  If [ABC]gpu != NULL and Copyto[ABC]gpu is true then data needs to be transferred to the GPU.

  If CopyfromCgpu flag is true and Cgpu!=NULL then copy from Cgpu to C. Otherwise leave
  data in Cgpu for reuse.

*/


template void RmgSymm<double>(char *, char *, int, int, double, double *, int, double *, int, 
                                  double, double *, int);
template void RmgSymm<std::complex<double> >(char *, char *, int, int, std::complex<double>, 
                      std::complex<double> *, int, std::complex<double> *, int, 
                      std::complex<double>, std::complex<double> *, int);

template <typename DataType> void RmgSymm(char *side, char *uplo, int m, int n, DataType alpha,
                             DataType *A, int lda, DataType *B, int ldb, DataType beta, DataType *C, int ldc)
{

#if BLAS_PROFILE
    if(typeid(DataType) == typeid(std::complex<double>))
    {
        if(pct.gridpe==0) printf("ZSYMM CALL m=%d n=%d\n",m,n);
    }
    else
    {
        if(pct.gridpe==0) printf("DSYMM CALL m=%d n=%d\n",m,n);
    }
#endif

#if CUDA_ENABLED
    cublasStatus_t custat;
    cublasSideMode_t cu_side = CUBLAS_SIDE_LEFT;
    cublasFillMode_t cu_uplo = CUBLAS_FILL_MODE_LOWER;


    if(!strcmp(side, "l")) cu_side = CUBLAS_SIDE_LEFT;
    if(!strcmp(side, "L")) cu_side = CUBLAS_SIDE_LEFT;
    if(!strcmp(side, "r")) cu_side = CUBLAS_SIDE_RIGHT;
    if(!strcmp(side, "R")) cu_side = CUBLAS_SIDE_RIGHT;

    if(!strcmp(uplo, "l")) cu_uplo = CUBLAS_FILL_MODE_LOWER;
    if(!strcmp(uplo, "L")) cu_uplo = CUBLAS_FILL_MODE_LOWER;
    if(!strcmp(uplo, "u")) cu_uplo = CUBLAS_FILL_MODE_UPPER;
    if(!strcmp(uplo, "U")) cu_uplo = CUBLAS_FILL_MODE_UPPER;

    int ka = m;
    if(!strcmp("r", side)) ka = n;
    if(!strcmp("R", side)) ka = n;

    if(ct.use_cublasxt && (typeid(DataType) == typeid(std::complex<double>)))
    {
        custat = cublasXtZsymm(ct.cublasxt_handle, cu_side, cu_uplo, m, n,
                            (cuDoubleComplex *)&alpha,
                            (cuDoubleComplex*)A, lda,
                            (cuDoubleComplex*)B, ldb,
                            (cuDoubleComplex*)&beta, (cuDoubleComplex*)C, ldc );
        ProcessGpublasError(custat);
        RmgGpuError(__FILE__, __LINE__, custat, "Problem executing cublasZgemm");
        return;
    }
    if(ct.use_cublasxt && (typeid(DataType) == typeid(double)))
    {
        custat = cublasXtDsymm(ct.cublasxt_handle, cu_side, cu_uplo, m, n,
                            (double*)&alpha,
                            (double*)A, lda,
                            (double*)B, ldb,
                            (double*)&beta, (double*)C, ldc );
        ProcessGpublasError(custat);
        RmgGpuError(__FILE__, __LINE__, custat, "Problem executing cublasDgemm");
        return;
    }

    DeviceSynchronize();
    cudaPointerAttributes attr;
    cudaError_t cudaerr;
    cudaerr = cudaPointerGetAttributes(&attr, A);
    bool a_dev = false;
#if (CUDA_VERSION_MAJOR > 10)
    if(cudaerr == cudaSuccess && attr.type == cudaMemoryTypeDevice) a_dev = true;
    cudaerr = cudaPointerGetAttributes(&attr, B);
    bool b_dev = false;
    if(cudaerr == cudaSuccess && attr.type == cudaMemoryTypeDevice) b_dev = true;
    cudaerr = cudaPointerGetAttributes(&attr, C);
    bool c_dev = false;
    if(cudaerr == cudaSuccess && attr.type == cudaMemoryTypeDevice) c_dev = true;
#else
    if(cudaerr == cudaSuccess && attr.memoryType == cudaMemoryTypeDevice) a_dev = true;
    cudaerr = cudaPointerGetAttributes(&attr, B);
    bool b_dev = false;
    if(cudaerr == cudaSuccess && attr.memoryType == cudaMemoryTypeDevice) b_dev = true;
    cudaerr = cudaPointerGetAttributes(&attr, C);
    bool c_dev = false;
    if(cudaerr == cudaSuccess && attr.memoryType == cudaMemoryTypeDevice) c_dev = true;
#endif

    size_t a_size = (size_t)lda * (size_t)ka;
    size_t b_size = (size_t)ldb * (size_t)n;
    size_t c_size = (size_t)ldc * (size_t)n;

    if(typeid(DataType) == typeid(std::complex<double>)) {
        std::complex<double> *dA=(std::complex<double> *)A, *dB=(std::complex<double> *)B, *dC=(std::complex<double> *)C;
        if(!a_dev) gpuMalloc((void **)&dA, a_size * sizeof(std::complex<double>));
        if(!b_dev) gpuMalloc((void **)&dB, b_size * sizeof(std::complex<double>));
        if(!c_dev) gpuMalloc((void **)&dC, c_size * sizeof(std::complex<double>));
        if(!a_dev) cudaMemcpy(dA, A, a_size * sizeof(std::complex<double>), cudaMemcpyDefault);
        if(!b_dev) cudaMemcpy(dB, B, b_size * sizeof(std::complex<double>), cudaMemcpyDefault);
        if(!c_dev && std::abs(beta) != 0.0) cudaMemcpy(dC, C, c_size * sizeof(std::complex<double>), cudaMemcpyDefault);
        custat = cublasZsymm(ct.cublas_handle, cu_side, cu_uplo, m, n,
                            (cuDoubleComplex *)&alpha,
                            (cuDoubleComplex*)dA, lda,
                            (cuDoubleComplex*)dB, ldb,
                            (cuDoubleComplex*)&beta, (cuDoubleComplex*)dC, ldc );
        ProcessGpublasError(custat);
        RmgGpuError(__FILE__, __LINE__, custat, "Problem executing cublasZgemm");
        if(!c_dev) cudaMemcpy(C, dC, c_size * sizeof(std::complex<double>), cudaMemcpyDefault);
        if(!c_dev) gpuFree(dC);
        if(!b_dev) gpuFree(dB);
        if(!a_dev) gpuFree(dA);
    }
    else {
        double *dA=(double *)A, *dB=(double *)B, *dC=(double *)C;
        if(!a_dev) gpuMalloc((void **)&dA, a_size * sizeof(double));
        if(!b_dev) gpuMalloc((void **)&dB, b_size * sizeof(double));
        if(!c_dev) gpuMalloc((void **)&dC, c_size * sizeof(double));
        if(!a_dev) cudaMemcpy(dA, A, a_size * sizeof(double), cudaMemcpyDefault);
        if(!b_dev) cudaMemcpy(dB, B, b_size * sizeof(double), cudaMemcpyDefault);
        if(!c_dev && beta != 0.0) cudaMemcpy(dC, C, c_size * sizeof(double), cudaMemcpyDefault);
        custat = cublasDsymm(ct.cublas_handle, cu_side, cu_uplo, m, n,
                            (double*)&alpha,
                            (double*)dA, lda,
                            (double*)dB, ldb,
                            (double*)&beta, (double*)dC, ldc );
        ProcessGpublasError(custat);
        RmgGpuError(__FILE__, __LINE__, custat, "Problem executing cublasDgemm");
        if(!c_dev) cudaMemcpy(C, dC, c_size * sizeof(double), cudaMemcpyDefault);
        if(!c_dev) gpuFree(dC);
        if(!b_dev) gpuFree(dB);
        if(!a_dev) gpuFree(dA);
    }
    DeviceSynchronize();

#elif HIP_ENABLED

    hipblasStatus_t custat;
    hipblasSideMode_t cu_side = HIPBLAS_SIDE_LEFT;
    hipblasFillMode_t cu_uplo = HIPBLAS_FILL_MODE_LOWER;


    if(!strcmp(side, "l")) cu_side = HIPBLAS_SIDE_LEFT;
    if(!strcmp(side, "L")) cu_side = HIPBLAS_SIDE_LEFT;
    if(!strcmp(side, "r")) cu_side = HIPBLAS_SIDE_RIGHT;
    if(!strcmp(side, "R")) cu_side = HIPBLAS_SIDE_RIGHT;

    if(!strcmp(uplo, "l")) cu_uplo = HIPBLAS_FILL_MODE_LOWER;
    if(!strcmp(uplo, "L")) cu_uplo = HIPBLAS_FILL_MODE_LOWER;
    if(!strcmp(uplo, "u")) cu_uplo = HIPBLAS_FILL_MODE_UPPER;
    if(!strcmp(uplo, "U")) cu_uplo = HIPBLAS_FILL_MODE_UPPER;

    int ka = m;
    if(!strcmp("r", side)) ka = n;
    if(!strcmp("R", side)) ka = n;

    DeviceSynchronize();
    hipPointerAttribute_t attr;
    hipError_t cudaerr;
    cudaerr = hipPointerGetAttributes(&attr, A);
    bool a_dev = false;

    if(cudaerr == hipSuccess && attr.memoryType == hipMemoryTypeDevice) a_dev = true;
    cudaerr = hipPointerGetAttributes(&attr, B);
    bool b_dev = false;
    if(cudaerr == hipSuccess && attr.memoryType == hipMemoryTypeDevice) b_dev = true;
    cudaerr = hipPointerGetAttributes(&attr, C);
    bool c_dev = false;
    if(cudaerr == hipSuccess && attr.memoryType == hipMemoryTypeDevice) c_dev = true;

    size_t a_size = (size_t)lda * (size_t)ka;
    size_t b_size = (size_t)ldb * (size_t)n;
    size_t c_size = (size_t)ldc * (size_t)n;

    if(typeid(DataType) == typeid(std::complex<double>)) {
        std::complex<double> *dA=(std::complex<double> *)A, *dB=(std::complex<double> *)B, *dC=(std::complex<double> *)C;
        if(!a_dev) gpuMalloc((void **)&dA, a_size * sizeof(std::complex<double>));
        if(!b_dev) gpuMalloc((void **)&dB, b_size * sizeof(std::complex<double>));
        if(!c_dev) gpuMalloc((void **)&dC, c_size * sizeof(std::complex<double>));
        if(!a_dev) hipMemcpy(dA, A, a_size * sizeof(std::complex<double>), hipMemcpyDefault);
        if(!b_dev) hipMemcpy(dB, B, b_size * sizeof(std::complex<double>), hipMemcpyDefault);
        if(!c_dev && std::abs(beta) != 0.0) hipMemcpy(dC, C, c_size * sizeof(std::complex<double>), hipMemcpyDefault);
        custat = hipblasZsymm(ct.hipblas_handle, cu_side, cu_uplo, m, n,
                            (hipblasDoubleComplex *)&alpha,
                            (hipblasDoubleComplex*)dA, lda,
                            (hipblasDoubleComplex*)dB, ldb,
                            (hipblasDoubleComplex*)&beta, (hipblasDoubleComplex*)dC, ldc );
        ProcessGpublasError(custat);
        RmgGpuError(__FILE__, __LINE__, custat, "Problem executing hipblasZsymm");
        if(!c_dev) hipMemcpy(C, dC, c_size * sizeof(std::complex<double>), hipMemcpyDefault);
        if(!c_dev) gpuFree(dC);
        if(!b_dev) gpuFree(dB);
        if(!a_dev) gpuFree(dA);
    }
    else {
        double *dA=(double *)A, *dB=(double *)B, *dC=(double *)C;
        if(!a_dev) gpuMalloc((void **)&dA, a_size * sizeof(double));
        if(!b_dev) gpuMalloc((void **)&dB, b_size * sizeof(double));
        if(!c_dev) gpuMalloc((void **)&dC, c_size * sizeof(double));
        if(!a_dev) hipMemcpy(dA, A, a_size * sizeof(double), hipMemcpyDefault);
        if(!b_dev) hipMemcpy(dB, B, b_size * sizeof(double), hipMemcpyDefault);
        if(!c_dev && beta != 0.0) hipMemcpy(dC, C, c_size * sizeof(double), hipMemcpyDefault);
        custat = hipblasDsymm(ct.hipblas_handle, cu_side, cu_uplo, m, n,
                            (double*)&alpha,
                            (double*)dA, lda,
                            (double*)dB, ldb,
                            (double*)&beta, (double*)dC, ldc );
        ProcessGpublasError(custat);
        RmgGpuError(__FILE__, __LINE__, custat, "Problem executing hipblasDsymm");
        if(!c_dev) gpuMemcpy(C, dC, c_size * sizeof(double), gpuMemcpyDefault);
        if(!c_dev) gpuFree(dC);
        if(!b_dev) gpuFree(dB);
        if(!a_dev) gpuFree(dA);
    }
    DeviceSynchronize();
#else

    if(typeid(DataType) == typeid(std::complex<double>)) {
        zsymm(side, uplo, &m, &n, (std::complex<double> *)&alpha, (std::complex<double> *)A, &lda, 
             (std::complex<double> *)B, &ldb, (std::complex<double> *)&beta, (std::complex<double> *)C, &ldc);
    }
    else {
        dsymm(side, uplo, &m, &n, (double *)&alpha, (double *)A, &lda, 
        (double *)B, &ldb, (double *)&beta, (double *)C, &ldc);
    }

#endif
}

