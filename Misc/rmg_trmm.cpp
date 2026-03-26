#include <complex>
#include <typeinfo>
#include <string.h>

#include "const.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "rmg_gemm.h"
#include "GpuAlloc.h"

#include "transition.h"
#include "rmg_error.h"
#include "rmg_dev_allocate.h"


#if CUDA_ENABLED
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cublas_v2.h>
#endif


/*
  These functions are used to hide the details of the matrix multiplication data types and GPU 
  utilization from the higher level routines.

*/

#define         dtrmm           RMG_FC_GLOBAL(dtrmm, DTRMM)
#define         ztrmm           RMG_FC_GLOBAL(ztrmm, ZTRMM)



#if SYCL_ENABLED
    #include <sycl/sycl.hpp>
    #include "oneapi/mkl/blas.hpp"
#else
extern "C" {
void dtrmm(const char *side, const char *uplo, const char *transa, const char *diag,
        int *m, int *n, double *alpha, double *a, int *lda, double *b, int *ldb);
void ztrmm(const char *side, const char *uplo, const char *transa, const char *diag,
        int *m, int *n, std::complex<double> *alpha, std::complex<double> *a, int *lda, std::complex<double> *b, int *ldb);
}
#endif

template void rmg::trmm<double>(char *, char *, char *, char *, int, int, double, double *, int,
                             double *, int);

template void rmg::trmm<std::complex<double> >(char *, char *, char *, char *, int, int,
                     std::complex<double>, std::complex<double> *, int,
                     std::complex<double> *, int);


template <typename DataType> void rmg::trmm(char *side, char *uplo, char *trans, char *diag,
                             int m, int n, DataType alpha, DataType *A, int lda,
                             DataType *B, int ldb)
{

#if CUDA_ENABLED

    cublasOperation_t cu_trans = CUBLAS_OP_N;
    cublasSideMode_t cu_side = CUBLAS_SIDE_LEFT;
    cublasFillMode_t fill_mode = CUBLAS_FILL_MODE_LOWER;
    cublasDiagType_t diag_mode = CUBLAS_DIAG_NON_UNIT;

    if(!strcmp(uplo, "u")) fill_mode = CUBLAS_FILL_MODE_UPPER;
    if(!strcmp(uplo, "U")) fill_mode = CUBLAS_FILL_MODE_UPPER;

    if(!strcmp(trans, "t")) cu_trans = CUBLAS_OP_T;
    if(!strcmp(trans, "T")) cu_trans = CUBLAS_OP_T;
    if(!strcmp(trans, "c")) cu_trans = CUBLAS_OP_C;
    if(!strcmp(trans, "C")) cu_trans = CUBLAS_OP_C;

    if(!strcmp(side, "l")) cu_side = CUBLAS_SIDE_LEFT;
    if(!strcmp(side, "L")) cu_side = CUBLAS_SIDE_LEFT;
    if(!strcmp(side, "r")) cu_side = CUBLAS_SIDE_RIGHT;
    if(!strcmp(side, "R")) cu_side = CUBLAS_SIDE_RIGHT;

    if(!strcmp(diag, "u")) diag_mode = CUBLAS_DIAG_UNIT;
    if(!strcmp(diag, "U")) diag_mode = CUBLAS_DIAG_UNIT;

    if(ct.use_cublasxt && (typeid(DataType) == typeid(std::complex<double>)))
    {
        rmg::error(cublasXtZtrmm(ct.cublasxt_handle, cu_side, fill_mode, cu_trans, diag_mode,
                            m, n,
                            (cuDoubleComplex *)&alpha,
                            (cuDoubleComplex *)A, lda,
                            (cuDoubleComplex *)B, ldb,
                            (cuDoubleComplex *)B, ldb));
        return;
    }
    if(ct.use_cublasxt && (typeid(DataType) == typeid(double)))
    {
        rmg::error(cublasXtDtrmm(ct.cublasxt_handle, cu_side, fill_mode, cu_trans, diag_mode,
                            m, n,
                            (double *)&alpha,
                            (double *)A, lda,
                            (double *)B, ldb,
                            (double *)B, ldb));
        return;
    }

    size_t a_size = (size_t)lda * (size_t)n;
    size_t b_size = (size_t)ldb * (size_t)n;

    cudaPointerAttributes attr;
    cudaError_t cudaerr;
    cudaerr = cudaPointerGetAttributes(&attr, A);
    bool a_dev = false;
#if (CUDA_VERSION_MAJOR > 10)
    if(cudaerr == cudaSuccess && attr.type == cudaMemoryTypeDevice) a_dev = true;
    cudaerr = cudaPointerGetAttributes(&attr, B);
    bool b_dev = false;
    if(cudaerr == cudaSuccess && attr.type == cudaMemoryTypeDevice) b_dev = true;
#else
    if(cudaerr == cudaSuccess && attr.type == cudaMemoryTypeDevice) a_dev = true;
    cudaerr = cudaPointerGetAttributes(&attr, B);
    bool b_dev = false;
    if(cudaerr == cudaSuccess && attr.type == cudaMemoryTypeDevice) b_dev = true;
#endif

    rmg::sync_device();
    if(typeid(DataType) == typeid(std::complex<double>)) {
        std::complex<double> *dA=(std::complex<double> *)A, *dB=(std::complex<double> *)B;
        if(!a_dev) rmg_device_pool->malloc(&dA, a_size);
        if(!b_dev) rmg_device_pool->malloc(&dB, b_size);
        if(!a_dev) cudaMemcpy(dA, A, a_size * sizeof(std::complex<double>), cudaMemcpyDefault);
        if(!b_dev) cudaMemcpy(dB, B, b_size * sizeof(std::complex<double>), cudaMemcpyDefault);
        rmg::error(cublasZtrmm(ct.cublas_handle, cu_side, fill_mode, cu_trans, diag_mode,
                            m, n,
                            (cuDoubleComplex *)&alpha,
                            (cuDoubleComplex *)dA, lda,
                            (cuDoubleComplex *)dB, ldb,
                            (cuDoubleComplex *)dB, ldb));
        if(!b_dev) cudaMemcpy(B, dB, b_size * sizeof(std::complex<double>), cudaMemcpyDefault);
        if(!b_dev) rmg_device_pool->free(dB);
        if(!a_dev) rmg_device_pool->free(dA);
    }
    else {
        double *dA=(double *)A, *dB=(double *)B;
        if(!a_dev) rmg_device_pool->malloc(&dA, a_size);
        if(!b_dev) rmg_device_pool->malloc(&dB, b_size);
        if(!a_dev) cudaMemcpy(dA, A, a_size * sizeof(double), cudaMemcpyDefault);
        if(!b_dev) cudaMemcpy(dB, B, b_size * sizeof(double), cudaMemcpyDefault);
        rmg::error(cublasDtrmm(ct.cublas_handle, cu_side, fill_mode, cu_trans, diag_mode,
                            m, n,
                            (double *)&alpha,
                            (double *)dA, lda,
                            (double *)dB, ldb,
                            (double *)dB, ldb));
        if(!b_dev) cudaMemcpy(B, dB, b_size * sizeof(double), cudaMemcpyDefault);
        if(!b_dev) rmg_device_pool->free(dB);
        if(!a_dev) rmg_device_pool->free(dA);
    }
    rmg::sync_device();
    return;

#elif HIP_ENABLED

    hipblasFillMode_t hip_fill = HIPBLAS_FILL_MODE_UPPER;
    hipblasSideMode_t hip_side = HIPBLAS_SIDE_RIGHT;
    hipblasDiagType_t hip_diag = HIPBLAS_DIAG_NON_UNIT;
    hipblasOperation_t hip_trans = HIPBLAS_OP_N;

    if(!strcmp(uplo, "l")) hip_fill = HIPBLAS_FILL_MODE_LOWER;
    if(!strcmp(uplo, "L")) hip_fill = HIPBLAS_FILL_MODE_LOWER;
    if(!strcmp(side, "l")) hip_side = HIPBLAS_SIDE_LEFT;
    if(!strcmp(side, "L")) hip_side = HIPBLAS_SIDE_LEFT;
    if(!strcmp(diag, "u")) hip_diag = HIPBLAS_DIAG_UNIT;
    if(!strcmp(diag, "U")) hip_diag = HIPBLAS_DIAG_UNIT;
    if(!strcmp(trans, "t")) hip_trans = HIPBLAS_OP_T;
    if(!strcmp(trans, "T")) hip_trans = HIPBLAS_OP_T;
    if(!strcmp(trans, "c")) hip_trans = HIPBLAS_OP_C;
    if(!strcmp(trans, "C")) hip_trans = HIPBLAS_OP_C;

    size_t a_size = (size_t)lda * (size_t)n;
    size_t b_size = (size_t)ldb * (size_t)n;

    hipPointerAttribute_t attr;
    hipError_t hiperr;
    hiperr = hipPointerGetAttributes(&attr, A);
    bool a_dev = false;
    if(hiperr == hipSuccess && attr.type == hipMemoryTypeDevice) a_dev = true;
    hiperr = hipPointerGetAttributes(&attr, B);
    bool b_dev = false;
    if(hiperr == hipSuccess && attr.type == hipMemoryTypeDevice) b_dev = true;

    if(typeid(DataType) == typeid(std::complex<double>)) {
        std::complex<double> *dA=(std::complex<double> *)A, *dB=(std::complex<double> *)B;
        if(!a_dev) rmg_device_pool->malloc(&dA, a_size);
        if(!b_dev) rmg_device_pool->malloc(&dB, b_size);
        if(!a_dev) rmg::error(hipMemcpyHtoD(dA, A, a_size * sizeof(std::complex<double>)));
        if(!b_dev) rmg::error(hipMemcpyHtoD(dB, B, b_size * sizeof(std::complex<double>)));
        rmg::error(hipblasZtrmm(ct.hipblas_handle,
                            hip_side, hip_fill, hip_trans, hip_diag,
                            m, n,
                            (hipDoubleComplex *)&alpha,
                            (hipDoubleComplex *)dA, lda,
                            (hipDoubleComplex *)dB, ldb,
                            (hipDoubleComplex *)dB, ldb));
        if(!b_dev) rmg::error(hipMemcpyDtoH(dB, B, b_size * sizeof(std::complex<double>)));
        if(!b_dev) rmg_device_pool->free(dB);
        if(!a_dev) rmg_device_pool->free(dA);
    }
    else {
        double *dA=(double *)A, *dB=(double *)B;
        if(!a_dev) rmg_device_pool->malloc(&dA, a_size);
        if(!b_dev) rmg_device_pool->malloc(&dB, b_size);
        if(!a_dev) rmg::error(hipMemcpyHtoD(dA, A, a_size * sizeof(double)));
        if(!b_dev) rmg::error(hipMemcpyHtoD(dB, B, b_size * sizeof(double)));
        rmg::error(hipblasDtrmm(ct.hipblas_handle,
                            hip_side, hip_fill, hip_trans, hip_diag,
                            m, n,
                            (double*)&alpha,
                            (double*)dA, lda,
                            (double*)dB, ldb,
                            (double*)dB, ldb ));
        if(!b_dev) rmg::error(hipMemcpyDtoH(B, dB, b_size * sizeof(double)));
        if(!b_dev) rmg_device_pool->free(dB);
        if(!a_dev) rmg_device_pool->free(dA);

    }
#elif SYCL_ENABLED

        // Determine side (left or right)
    oneapi::mkl::side side_mode = oneapi::mkl::side::left;
    if (!strcmp(side, "r")) side_mode = oneapi::mkl::side::right;
    if (!strcmp(side, "R")) side_mode = oneapi::mkl::side::right;

    // Determine upper/lower fill mode
    oneapi::mkl::uplo fill_mode = oneapi::mkl::uplo::lower;
    if (!strcmp(uplo, "u")) fill_mode = oneapi::mkl::uplo::upper;
    if (!strcmp(uplo, "U")) fill_mode = oneapi::mkl::uplo::upper;

    // Determine transpose operation
    oneapi::mkl::transpose sycl_transA = oneapi::mkl::transpose::nontrans;
    if (!strcmp(trans, "t")) sycl_transA = oneapi::mkl::transpose::trans;
    if (!strcmp(trans, "T")) sycl_transA = oneapi::mkl::transpose::trans;
    if (!strcmp(trans, "c")) sycl_transA = oneapi::mkl::transpose::conjtrans; 
    if (!strcmp(trans, "C")) sycl_transA = oneapi::mkl::transpose::conjtrans;

    // Determine whether the matrix is unit triangular
    oneapi::mkl::diag diag_mode = oneapi::mkl::diag::nonunit;
    if (!strcmp(diag, "u")) diag_mode = oneapi::mkl::diag::unit;
    if (!strcmp(diag, "U")) diag_mode = oneapi::mkl::diag::unit;

    // TRMM: Triangular matrix-matrix multiplication
    // B = alpha * op(A) * B  or  B = alpha * B * op(A)
    try {
        sycl::event trmm_event = oneapi::mkl::blas::column_major::trmm(
            ct.sycl_Q,
            side_mode,
            fill_mode,
            sycl_transA,
            diag_mode,
            m, n,
            alpha,
            A, lda,
            B, ldb);
        // Wait for completion
        trmm_event.wait();
    }
    catch (sycl::exception const& e) {
        std::cout << "\t\tCaught synchronous SYCL exception during TRMM:\n"
                  << e.what() << std::endl << std::endl;
        rmg::error("Terminating");
    }


#else

    if(typeid(DataType) == typeid(std::complex<double>)) {
        ztrmm(side, uplo, trans, diag, &m, &n, (std::complex<double> *)&alpha,
             (std::complex<double> *)A, &lda, (std::complex<double> *)B, &ldb);
    }
    else {
        dtrmm(side, uplo, trans, diag, &m, &n, (double *)&alpha,
             (double *)A, &lda, (double *)B, &ldb);
    }

#endif
}

