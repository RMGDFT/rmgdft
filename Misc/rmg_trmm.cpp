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

#define         dtrmm           RMG_FC_GLOBAL(dtrmm, DTRMM)
#define         ztrmm           RMG_FC_GLOBAL(ztrmm, ZTRMM)



#if SYCL_ENABLED
    #include <CL/sycl.hpp>
    #include "oneapi/mkl/blas.hpp"
    #include "mkl.h"
#else
extern "C" {
void dtrmm(const char *side, const char *uplo, const char *transa, const char *diag,
        int *m, int *n, double *alpha, double *a, int *lda, double *b, int *ldb);
void ztrmm(const char *side, const char *uplo, const char *transa, const char *diag,
        int *m, int *n, std::complex<double> *alpha, std::complex<double> *a, int *lda, std::complex<double> *b, int *ldb);
}
#endif

template void rmg_trmm<double>(char *, char *, char *, char *, int, int, double, double *, int,
                             double *, int);

template void rmg_trmm<std::complex<double> >(char *, char *, char *, char *, int, int,
                     std::complex<double>, std::complex<double> *, int,
                     std::complex<double> *, int);


template <typename DataType> void rmg_trmm(char *side, char *uplo, char *trans, char *diag,
                             int m, int n, DataType alpha, DataType *A, int lda,
                             DataType *B, int ldb)
{

#if CUDA_ENABLED

    cublasStatus_t custat;
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
        custat = cublasXtZtrmm(ct.cublasxt_handle, cu_side, fill_mode, cu_trans, diag_mode,
                            m, n,
                            (cuDoubleComplex *)&alpha,
                            (cuDoubleComplex *)dA, lda,
                            (cuDoubleComplex *)dB, ldb,
                            (cuDoubleComplex *)dB, ldb);
        ProcessGpublasError(custat);
        RmgGpuError(__FILE__, __LINE__, custat, "Problem executing cublasXtZtrmm");
        return;
    }
    if(ct.use_cublasxt && (typeid(DataType) == typeid(double)))
    {
        custat = cublasXtDtrmm(ct.cublasxt_handle, cu_side, fill_mode, cu_trans, diag_mode,
                            m, n,
                            (double *)&alpha,
                            (double *)dA, lda,
                            (double *)dB, ldb,
                            (double *)dB, ldb);
        ProcessGpublasError(custat);
        RmgGpuError(__FILE__, __LINE__, custat, "Problem executing cublasXtDtrmm");
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

    DeviceSynchronize();
    if(typeid(DataType) == typeid(std::complex<double>)) {
        std::complex<double> *dA=(std::complex<double> *)A, *dB=(std::complex<double> *)B;
        if(!a_dev) gpuMalloc((void **)&dA, a_size * sizeof(std::complex<double>));
        if(!b_dev) gpuMalloc((void **)&dB, b_size * sizeof(std::complex<double>));
        if(!a_dev) cudaMemcpy(dA, A, a_size * sizeof(std::complex<double>), cudaMemcpyDefault);
        if(!b_dev) cudaMemcpy(dB, B, b_size * sizeof(std::complex<double>), cudaMemcpyDefault);
        custat = cublasZtrmm(ct.cublas_handle, cu_side, fill_mode, cu_trans, diag_mode,
                            m, n,
                            (cuDoubleComplex *)&alpha,
                            (cuDoubleComplex *)dA, lda,
                            (cuDoubleComplex *)dB, ldb,
                            (cuDoubleComplex *)dB, ldb);
        ProcessGpublasError(custat);
        RmgGpuError(__FILE__, __LINE__, custat, "Problem executing cublasZtrmm");
        if(!b_dev) cudaMemcpy(B, dB, b_size * sizeof(std::complex<double>), cudaMemcpyDefault);
        if(!b_dev) gpuFree(dB);
        if(!a_dev) gpuFree(dA);
    }
    else {
        double *dA=(double *)A, *dB=(double *)B;
        if(!a_dev) gpuMalloc((void **)&dA, a_size * sizeof(double));
        if(!b_dev) gpuMalloc((void **)&dB, b_size * sizeof(double));
        if(!a_dev) cudaMemcpy(dA, A, a_size * sizeof(double), cudaMemcpyDefault);
        if(!b_dev) cudaMemcpy(dB, B, b_size * sizeof(double), cudaMemcpyDefault);
        custat = cublasDtrmm(ct.cublas_handle, cu_side, fill_mode, cu_trans, diag_mode,
                            m, n,
                            (cuDoubleComplex *)&alpha,
                            (cuDoubleComplex *)dA, lda,
                            (cuDoubleComplex *)dB, ldb,
                            (cuDoubleComplex *)dB, ldb);
        ProcessGpublasError(custat);
        RmgGpuError(__FILE__, __LINE__, custat, "Problem executing cublasDtrmm");
        if(!b_dev) cudaMemcpy(B, dB, b_size * sizeof(double), cudaMemcpyDefault);
        if(!b_dev) gpuFree(dB);
        if(!a_dev) gpuFree(dA);
    }
    DeviceSynchronize();
    return;

#elif HIP_ENABLED
    rocblas_status rocstat;
    rocblas_fill fill_mode = rocblas_fill_lower;
    rocblas_diagonal diag_mode = rocblas_diagonal_non_unit;
    rocblas_side side_mode = rocblas_side_left;
    rocblas_operation roc_trans = rocblas_operation_none;

    if(!strcmp(uplo, "u")) fill_mode = rocblas_fill_upper;
    if(!strcmp(uplo, "U")) fill_mode = rocblas_fill_upper;

    if(!strcmp(trans, "t")) roc_trans = rocblas_operation_transpose;
    if(!strcmp(trans, "T")) roc_trans = rocblas_operation_transpose;
    if(!strcmp(trans, "c")) roc_trans = rocblas_operation_conjugate_transpose;
    if(!strcmp(trans, "C")) roc_trans = rocblas_operation_conjugate_transpose;

    if(!strcmp(side, "r")) side_mode = rocblas_side_right;
    if(!strcmp(side, "R")) side_mode = rocblas_side_right;

    if(!strcmp(diag, "u")) diag_mode = rocblas_diagonal_unit;
    if(!strcmp(diag, "U")) diag_mode = rocblas_diagonal_unit;


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
        if(!a_dev) gpuMalloc((void **)&dA, a_size * sizeof(std::complex<double>));
        if(!b_dev) gpuMalloc((void **)&dB, b_size * sizeof(std::complex<double>));
        if(!a_dev) hipMemcpyHtoD(dA, A, a_size * sizeof(std::complex<double>));
        if(!b_dev) hipMemcpyHtoD(dB, B, b_size * sizeof(std::complex<double>));
        rocstat = rocblas_ztrmm(ct.roc_handle, side_mode, fill_mode, roc_trans, diag_mode,
                            m, n,
                            (rocblas_double_complex *)&alpha,
                            (rocblas_double_complex *)dA, lda,
                            (rocblas_double_complex *)dB, ldb,
                            (rocblas_double_complex *)dB, ldb);
//        ProcessGpublasError(hipstat);
//        RmgGpuError(__FILE__, __LINE__, rocstat, "Problem executing cublasZsyrkx");
        if(!b_dev) hipMemcpyDtoH(dB, B, b_size * sizeof(std::complex<double>));
        if(!b_dev) gpuFree(dB);
        if(!a_dev) gpuFree(dA);
    }
    else {
        double *dA=(double *)A, *dB=(double *)B;
        if(!a_dev) gpuMalloc((void **)&dA, a_size * sizeof(double));
        if(!b_dev) gpuMalloc((void **)&dB, b_size * sizeof(double));
        if(!a_dev) hipMemcpyHtoD(dA, A, a_size * sizeof(double));
        if(!b_dev) hipMemcpyHtoD(dB, B, b_size * sizeof(double));
        rocstat = rocblas_dtrmm(ct.roc_handle, side_mode, fill_mode, roc_trans, diag_mode,
                            m, n,
                            (double*)&alpha,
                            (double*)dA, lda,
                            (double*)dB, ldb,
                            (double*)dB, ldb );
        //ProcessGpublasError(hipstat);
        //RmgGpuError(__FILE__, __LINE__, rocstat, "Problem executing hipblasDsyrkx");
        if(!b_dev) hipMemcpyDtoH(B, dB, b_size * sizeof(double));
        if(!b_dev) gpuFree(dB);
        if(!a_dev) gpuFree(dA);

    }
#elif SYCL_ENABLED

this should cause a compile error since as I have no access to a machine to test this on right now
    oneapi::mkl::uplo fill_mode = oneapi::mkl::uplo::lower;
    if(!strcmp(uplo, "u")) fill_mode = oneapi::mkl::uplo::upper;
    if(!strcmp(uplo, "U")) fill_mode = oneapi::mkl::uplo::upper;

    oneapi::mkl::transpose sycl_transA = oneapi::mkl::transpose::nontrans, sycl_transB;

    if(!strcmp(trans, "t")) sycl_transA = oneapi::mkl::transpose::trans;
    if(!strcmp(trans, "T")) sycl_transA = oneapi::mkl::transpose::trans;
    if(!strcmp(trans, "c")) sycl_transA = oneapi::mkl::transpose::conjtrans;
    if(!strcmp(trans, "C")) sycl_transA = oneapi::mkl::transpose::conjtrans;
    if(sycl_transA == oneapi::mkl::transpose::nontrans)
    {
	    if(!strcmp(trans, "t")) sycl_transB = oneapi::mkl::transpose::trans;
	    if(!strcmp(trans, "T")) sycl_transB = oneapi::mkl::transpose::trans;
	    if(!strcmp(trans, "c")) sycl_transB = oneapi::mkl::transpose::conjtrans;
	    if(!strcmp(trans, "C")) sycl_transB = oneapi::mkl::transpose::conjtrans;
    }
    else
    {
        sycl_transB = oneapi::mkl::transpose::nontrans;
    }

    size_t a_size = (size_t)lda * (size_t)n;
    size_t b_size = (size_t)ldb * (size_t)n;

    cl::sycl::buffer<DataType, 1> bufA((DataType *)A, a_size, {cl::sycl::property::buffer::use_host_ptr()});
    bufA.set_final_data(nullptr);
    cl::sycl::buffer<DataType, 1> bufB((DataType *)B, b_size, {cl::sycl::property::buffer::use_host_ptr()});
    if(A == B)
    {
        try {
            oneapi::mkl::blas::gemmt(ct.sycl_Q, fill_mode, sycl_transA, sycl_transB, n, k, alpha, 
                                    bufA, lda, bufA, ldb, beta, bufB, ldb);
        }
        catch(cl::sycl::exception const& e) {
            std::cout << "\t\tCaught synchronous SYCL exception during GEMMT:\n"
            << e.what() << std::endl << std::endl;
            rmg_error_handler (__FILE__, __LINE__, "Terminating");
        }
    }
    else
    {
        cl::sycl::buffer<DataType, 1> bufB((DataType *)B, b_size, {cl::sycl::property::buffer::use_host_ptr()});
        bufB.set_final_data(nullptr);
        try {
            oneapi::mkl::blas::gemmt(ct.sycl_Q, fill_mode, sycl_transA, sycl_transB, n, k, alpha, 
                                    bufA, lda, bufB, ldb, beta, bufB, ldb);
        }
        catch(cl::sycl::exception const& e) {
            std::cout << "\t\tCaught synchronous SYCL exception during GEMMT:\n"
            << e.what() << std::endl << std::endl;
            rmg_error_handler (__FILE__, __LINE__, "Terminating");
        }
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

