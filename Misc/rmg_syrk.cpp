#include <complex>
#include <typeinfo>
#include <string.h>

#include "const.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "rmg_gemm.h"
#include "rmg_dev_allocate.h"
#include "GpuAlloc.h"

#include "transition.h"
#include "rmg_error.h"


#if CUDA_ENABLED
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cublas_v2.h>
#endif


/*
  These functions are used to hide the details of the matrix multiplication data types and GPU 
  utilization from the higher level routines.

*/

#define         dsyrk           RMG_FC_GLOBAL(dsyrk, DSYRK)
#define         zsyrk           RMG_FC_GLOBAL(zsyrk, ZSYRK)



#if SYCL_ENABLED
    #include <sycl/sycl.hpp>
    #include "oneapi/mkl/blas.hpp"
#else
extern "C" {
void dsyrk(const char *, const char *, int *, int *, double *, double *, int *, double *, double *, int *);
void zsyrk(const char *, const char *, int *, int *, std::complex<double> *, std::complex<double> *, int *, std::complex<double> *, std::complex<double> *, int *);
}
#endif

template void rmg::syrk<double>(char *, char *, int, int, double, double *, int,
                             double, double *, int);

template void rmg::syrk<std::complex<double> >(char *, char *, int, int, std::complex<double>, 
                    std::complex<double> *, int, std::complex<double>, std::complex<double> *, int);


template <typename DataType> void rmg::syrk(char *uplo, char *trans, int n, int k, 
                             DataType alpha, DataType *A, int lda, DataType beta, 
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

    cublasOperation_t cu_trans = CUBLAS_OP_N;
    cublasFillMode_t fill_mode = CUBLAS_FILL_MODE_LOWER;

    if(!strcmp(uplo, "u")) fill_mode = CUBLAS_FILL_MODE_UPPER;
    if(!strcmp(uplo, "U")) fill_mode = CUBLAS_FILL_MODE_UPPER;

    if(!strcmp(trans, "t")) cu_trans = CUBLAS_OP_T;
    if(!strcmp(trans, "T")) cu_trans = CUBLAS_OP_T;
    if(!strcmp(trans, "c")) cu_trans = CUBLAS_OP_C;
    if(!strcmp(trans, "C")) cu_trans = CUBLAS_OP_C;

    if(ct.use_cublasxt && (typeid(DataType) == typeid(std::complex<double>)))
    {
        rmg::error(cublasXtZsyrk(ct.cublasxt_handle, fill_mode, cu_trans, n, k,
                            (cuDoubleComplex *)&alpha,
                            (cuDoubleComplex*)A, lda,
                            (cuDoubleComplex*)&beta, (cuDoubleComplex*)C, ldc ));
        return;
    }
    if(ct.use_cublasxt && (typeid(DataType) == typeid(double)))
    {
        rmg::error(cublasXtDsyrk(ct.cublasxt_handle, fill_mode, cu_trans, (size_t)n, (size_t)k,
                            (double*)&alpha,
                            (double*)A, (size_t)lda,
                            (double*)&beta, (double*)C, (size_t)ldc ));
        return;
    }

    size_t a_size = (size_t)lda * (size_t)n;
    size_t c_size = (size_t)ldc * (size_t)n;

    cudaPointerAttributes attr;
    cudaError_t cudaerr;
    cudaerr = cudaPointerGetAttributes(&attr, A);
    bool a_dev = false;
#if (CUDA_VERSION_MAJOR > 10)
    if(cudaerr == cudaSuccess && attr.type == cudaMemoryTypeDevice) a_dev = true;
    cudaerr = cudaPointerGetAttributes(&attr, C);
    bool c_dev = false;
    if(cudaerr == cudaSuccess && attr.type == cudaMemoryTypeDevice) c_dev = true;
#else
    if(cudaerr == cudaSuccess && attr.type == cudaMemoryTypeDevice) a_dev = true;
    cudaerr = cudaPointerGetAttributes(&attr, C);
    bool c_dev = false;
    if(cudaerr == cudaSuccess && attr.type == cudaMemoryTypeDevice) c_dev = true;
#endif

    rmg::sync_device();
    if(typeid(DataType) == typeid(std::complex<double>)) {
        std::complex<double> *dA=(std::complex<double> *)A, *dC=(std::complex<double> *)C;
        if(!a_dev) rmg_device_pool->malloc(&dA, a_size);
        if(!c_dev) rmg_device_pool->malloc(&dC, c_size);
        if(!a_dev) cudaMemcpy(dA, A, a_size * sizeof(std::complex<double>), cudaMemcpyDefault);
        if(!c_dev && std::abs(beta) != 0.0) cudaMemcpy(dC, C, c_size * sizeof(std::complex<double>), cudaMemcpyDefault);
        rmg::error(cublasZsyrk(ct.cublas_handle, fill_mode, cu_trans, n, k,
                            (cuDoubleComplex *)&alpha,
                            (cuDoubleComplex*)dA, lda,
                            (cuDoubleComplex*)&beta, (cuDoubleComplex*)dC, ldc ));
        if(!c_dev) cudaMemcpy(C, dC, c_size * sizeof(std::complex<double>), cudaMemcpyDefault);
        if(!c_dev) rmg_device_pool->free(dC);
        if(!a_dev) rmg_device_pool->free(dA);
    }
    else {
        double *dA=(double *)A, *dC=(double *)C;
        if(!a_dev) rmg_device_pool->malloc(&dA, a_size);
        if(!c_dev) rmg_device_pool->malloc(&dC, c_size);
        if(!a_dev) cudaMemcpy(dA, A, a_size * sizeof(double), cudaMemcpyDefault);
        if(!c_dev && beta != 0.0) cudaMemcpy(dC, C, c_size * sizeof(double), cudaMemcpyDefault);
        rmg::error(cublasDsyrk(ct.cublas_handle, fill_mode, cu_trans, n, k,
                            (double*)&alpha,
                            (double*)dA, lda,
                            (double*)&beta, (double*)dC, ldc ));
        if(!c_dev) cudaMemcpy(C, dC, c_size * sizeof(double), cudaMemcpyDefault);
        if(!c_dev) rmg_device_pool->free(dC);
        if(!a_dev) rmg_device_pool->free(dA);
    }
    rmg::sync_device();
    return;

#elif HIP_ENABLED
    hipblasOperation_t hip_trans = HIPBLAS_OP_N;
    hipblasFillMode_t fill_mode = HIPBLAS_FILL_MODE_LOWER;

    if(!strcmp(uplo, "u")) fill_mode = HIPBLAS_FILL_MODE_UPPER;
    if(!strcmp(uplo, "U")) fill_mode = HIPBLAS_FILL_MODE_UPPER;

    if(!strcmp(trans, "t")) hip_trans = HIPBLAS_OP_T;
    if(!strcmp(trans, "T")) hip_trans = HIPBLAS_OP_T;
    if(!strcmp(trans, "c")) hip_trans = HIPBLAS_OP_C;
    if(!strcmp(trans, "C")) hip_trans = HIPBLAS_OP_C;

    size_t a_size = (size_t)lda * (size_t)n;
    size_t c_size = (size_t)ldc * (size_t)n;

    hipPointerAttribute_t attr;
    hipError_t hiperr;
    hiperr = hipPointerGetAttributes(&attr, A);
    bool a_dev = false;
    if(hiperr == hipSuccess && attr.type == hipMemoryTypeDevice) a_dev = true;
    hiperr = hipPointerGetAttributes(&attr, C);
    bool c_dev = false;
    if(hiperr == hipSuccess && attr.type == hipMemoryTypeDevice) c_dev = true;

    if(typeid(DataType) == typeid(std::complex<double>)) {
        std::complex<double> *dA=(std::complex<double> *)A, *dC=(std::complex<double> *)C;
        if(!a_dev) rmg_device_pool->malloc(&dA, a_size);
        if(!c_dev) rmg_device_pool->malloc(&dC, c_size);
        if(!a_dev) rmg::error(hipMemcpyHtoD(dA, A, a_size * sizeof(std::complex<double>)));
        if(!c_dev && std::abs(beta) != 0.0) rmg::error(hipMemcpyHtoD(dC, C, c_size * sizeof(std::complex<double>)));
        rmg::error(hipblasZsyrk(ct.hipblas_handle, fill_mode, hip_trans, n, k,
                            (hipDoubleComplex *)&alpha,
                            (hipDoubleComplex*)dA, lda,
                            (hipDoubleComplex*)&beta, (hipDoubleComplex*)dC, ldc ));
        if(!c_dev) rmg::error(hipMemcpyDtoH(dC, C, c_size * sizeof(std::complex<double>)));
        if(!c_dev) rmg_device_pool->free(dC);
        if(!a_dev) rmg_device_pool->free(dA);
    }
    else {
        double *dA=(double *)A, *dC=(double *)C;
        if(!a_dev) rmg_device_pool->malloc(&dA, a_size);
        if(!c_dev) rmg_device_pool->malloc(&dC, c_size);
        if(!a_dev) rmg::error(hipMemcpyHtoD(dA, A, a_size * sizeof(double)));
        if(!c_dev && beta != 0.0) rmg::error(hipMemcpyHtoD(dC, C, c_size * sizeof(double)));
        rmg::error(hipblasDsyrk(ct.hipblas_handle, fill_mode, hip_trans, n, k,
                            (double*)&alpha,
                            (double*)dA, lda,
                            (double*)&beta, (double*)dC, ldc ));
        if(!c_dev) rmg::error(hipMemcpyDtoH(C, dC, c_size * sizeof(double)));
        if(!c_dev) rmg_device_pool->free(dC);
        if(!a_dev) rmg_device_pool->free(dA);

    }
#elif SYCL_ENABLED

    // Determine upper/lower fill mode
    oneapi::mkl::uplo fill_mode = oneapi::mkl::uplo::lower;
    if(!strcmp(uplo, "u")) fill_mode = oneapi::mkl::uplo::upper;
    if(!strcmp(uplo, "U")) fill_mode = oneapi::mkl::uplo::upper;

    // Determine transpose operation
    oneapi::mkl::transpose sycl_transA = oneapi::mkl::transpose::nontrans;

    if(!strcmp(trans, "t")) sycl_transA = oneapi::mkl::transpose::trans;
    if(!strcmp(trans, "T")) sycl_transA = oneapi::mkl::transpose::trans;
    if(!strcmp(trans, "c")) sycl_transA = oneapi::mkl::transpose::conjtrans;
    if(!strcmp(trans, "C")) sycl_transA = oneapi::mkl::transpose::conjtrans;

    // SYRK: C = alpha * op(A) * op(A)^T + beta * C
    // No separate B matrix needed — SYRK uses A for both operands
    try {
        sycl::event syrk_event = oneapi::mkl::blas::column_major::syrk(
            ct.sycl_Q,
            fill_mode,
            sycl_transA,
            n, k,            
            alpha,
            A, lda,
            beta,
            C, ldc);
        // Wait for completion
        syrk_event.wait();
    }
    catch(sycl::exception const& e) {
        std::cout << "\t\tCaught synchronous SYCL exception during SYRK:\n"
                  << e.what() << std::endl << std::endl;
        rmg::error("Terminating");
    }


#else

    if(typeid(DataType) == typeid(std::complex<double>)) {
        zsyrk(uplo, trans, &n, &k, (std::complex<double> *)&alpha, (std::complex<double> *)A, &lda,
             (std::complex<double> *)&beta, (std::complex<double> *)C, &ldc);
    }
    else {
        dsyrk(uplo, trans, &n, &k, (double *)&alpha, (double *)A, &lda,
              (double *)&beta, (double *)C, &ldc);
    }

#endif
}

