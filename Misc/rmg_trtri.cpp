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

#define         dtrtri          RMG_FC_GLOBAL(dtrtri, DTRTRI)
#define         ztrtri          RMG_FC_GLOBAL(ztrtri, ZTRTRI)

extern "C" {
void dtrtri(const char *uplo, const char *diag, int *n, double *a, int *lda, int *info );
void ztrtri(const char *uplo, const char *diag, int *n, std::complex<double> *a, int *lda, int *info );
};


/*
  These functions are used to hide the details of the xtrtri calls and GPU
  interface from the higher level routines.

*/

#if SYCL_ENABLED
    #include <CL/sycl.hpp>
    #include "oneapi/mkl/blas.hpp"
    #include "mkl.h"
#endif

template void rmg_trtri<double>(char *, char *, int, double *, int, int *);
template void rmg_trtri<std::complex<double>>(char *, char *, int, std::complex<double> *, int, int *);


template <typename DataType> void rmg_trtri(char *uplo, char *diag, int n, DataType *A, int lda, int *info)
{

#if CUDA_ENABLED
    int *dev_info;
    cusolverStatus_t custat;
    cublasFillMode_t fill_mode = CUBLAS_FILL_MODE_LOWER;
    cublasDiagType_t diag_mode = CUBLAS_DIAG_NON_UNIT;

    if(!strcmp(uplo, "u")) fill_mode = CUBLAS_FILL_MODE_UPPER;
    if(!strcmp(uplo, "U")) fill_mode = CUBLAS_FILL_MODE_UPPER;
    if(!strcmp(diag, "u")) diag_mode = CUBLAS_DIAG_UNIT;
    if(!strcmp(diag, "U")) diag_mode = CUBLAS_DIAG_UNIT;

    size_t a_size = (size_t)lda * (size_t)n;
    gpuMalloc((void **)&dev_info, sizeof(int));
    cudaPointerAttributes attr;
    cudaError_t cudaerr;
    cudaerr = cudaPointerGetAttributes(&attr, A);
    bool a_dev = false;
#if (CUDA_VERSION_MAJOR > 10)
    if(cudaerr == cudaSuccess && attr.type == cudaMemoryTypeDevice) a_dev = true;
#else
    if(cudaerr == cudaSuccess && attr.type == cudaMemoryTypeDevice) a_dev = true;
#endif

    DeviceSynchronize();
    if(typeid(DataType) == typeid(std::complex<double>)) {
        size_t dwork;
        size_t hwork;
        char *hbuf;
	void *dbuf;

        std::complex<double> *dA=(std::complex<double> *)A;
        if(!a_dev) gpuMalloc((void **)&dA, a_size * sizeof(std::complex<double>));
        if(!a_dev) cudaMemcpy(dA, A, a_size * sizeof(std::complex<double>), cudaMemcpyDefault);
        custat = cusolverDnXtrtri_bufferSize(ct.cusolver_handle, fill_mode, diag_mode,
            n, CUDA_C_64F, (void *)dA, lda, &dwork, &hwork);
        gpuMalloc((void **)&dbuf, dwork);
        hbuf = new char[hwork];
        custat = cusolverDnXtrtri(ct.cusolver_handle, fill_mode, diag_mode,
                 n, CUDA_C_64F, (void *)dA, lda, (void *)dbuf, dwork, (void *)hbuf, hwork, dev_info);

	if(custat != CUSOLVER_STATUS_SUCCESS)
            rmg_error_handler (__FILE__, __LINE__, " cusolverDnZtrtri failed.");

	delete [] hbuf;
	gpuFree(dbuf);
        if(!a_dev) cudaMemcpy(A, dA, a_size * sizeof(std::complex<double>), cudaMemcpyDefault);
        if(!a_dev) gpuFree(dA);
    }
    else {
        size_t dwork;
        size_t hwork;
        char *hbuf;
	void *dbuf;

        double *dA=(double *)A;
        if(!a_dev) gpuMalloc((void **)&dA, a_size * sizeof(double));
        if(!a_dev) cudaMemcpy(dA, A, a_size * sizeof(double), cudaMemcpyDefault);
        custat = cusolverDnXtrtri_bufferSize(ct.cusolver_handle, fill_mode, diag_mode,
            n, CUDA_R_64F, (void *)dA, lda, &dwork, &hwork);
        gpuMalloc((void **)&dbuf, dwork);
        hbuf = new char[hwork];
        custat = cusolverDnXtrtri(ct.cusolver_handle, fill_mode, diag_mode,
                 n, CUDA_R_64F, (void *)dA, lda, (void *)dbuf, dwork, (void *)hbuf, hwork, dev_info);

	if(custat != CUSOLVER_STATUS_SUCCESS)
            rmg_error_handler (__FILE__, __LINE__, " cusolverDnDtrtri failed.");

	delete [] hbuf;
	gpuFree(dbuf);
        if(!a_dev) cudaMemcpy(A, dA, a_size * sizeof(double), cudaMemcpyDefault);
        if(!a_dev) gpuFree(dA);
    }
    cudaMemcpy(info, dev_info, sizeof(int), cudaMemcpyDefault);
    DeviceSynchronize();
    gpuFree(dev_info);
    return;

#elif HIP_ENABLED
    rocblas_status rocstat;
    rocblas_fill fill_mode = rocblas_fill_lower;
    rocblas_diagonal diag_mode = rocblas_diagonal_non_unit;

    if(!strcmp(uplo, "u")) fill_mode = rocblas_fill_upper;
    if(!strcmp(uplo, "U")) fill_mode = rocblas_fill_upper;
    if(!strcmp(diag, "u")) diag_mode = rocblas_diagonal_unit;
    if(!strcmp(diag, "U")) diag_mode = rocblas_diagonal_unit;

    size_t a_size = (size_t)lda * (size_t)n;

    hipPointerAttribute_t attr;
    hipError_t hiperr;
    hiperr = hipPointerGetAttributes(&attr, A);
    bool a_dev = false;
    int *dev_info;
    gpuMalloc((void **)&dev_info, sizeof(int));
    if(hiperr == hipSuccess && attr.type == hipMemoryTypeDevice) a_dev = true;

    if(typeid(DataType) == typeid(std::complex<double>)) {
        std::complex<double> *dA=(std::complex<double> *)A;
        if(!a_dev) gpuMalloc((void **)&dA, a_size * sizeof(std::complex<double>));
        if(!a_dev) hipMemcpyHtoD(dA, A, a_size * sizeof(std::complex<double>));
        rocstat = rocsolver_ztrtri(ct.roc_handle, fill_mode, diag_mode, n,
                                   (rocblas_double_complex *)dA, lda, dev_info);
        if (rocstat != rocblas_status_success) 
            rmg_error_handler(__FILE__, __LINE__, "Problem executing rocsolver_ztrtri");
        if(!a_dev) hipMemcpyDtoH(A, dA, a_size * sizeof(std::complex<double>));
        if(!a_dev) gpuFree(dA);
    }
    else {
        double *dA=(double *)A;
        if(!a_dev) gpuMalloc((void **)&dA, a_size * sizeof(double));
        if(!a_dev) hipMemcpyHtoD(dA, A, a_size * sizeof(double));
        rocstat = rocsolver_dtrtri(ct.roc_handle, fill_mode, diag_mode, n, dA, lda, dev_info);
        if (rocstat != rocblas_status_success) 
            rmg_error_handler(__FILE__, __LINE__, "Problem executing rocsolver_dtrtri");
        if(!a_dev) hipMemcpyDtoH(A, dA, a_size * sizeof(double));
        if(!a_dev) gpuFree(dA);
    }
    hipMemcpyDtoH(info, dev_info, sizeof(int));
    gpuFree(dev_info);
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
    size_t c_size = (size_t)ldc * (size_t)n;

    cl::sycl::buffer<DataType, 1> bufA((DataType *)A, a_size, {cl::sycl::property::buffer::use_host_ptr()});
    bufA.set_final_data(nullptr);
    cl::sycl::buffer<DataType, 1> bufC((DataType *)C, c_size, {cl::sycl::property::buffer::use_host_ptr()});
    if(A == B)
    {
        try {
            oneapi::mkl::blas::gemmt(ct.sycl_Q, fill_mode, sycl_transA, sycl_transB, n, k, alpha, 
                                    bufA, lda, bufA, ldb, beta, bufC, ldc);
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
                                    bufA, lda, bufB, ldb, beta, bufC, ldc);
        }
        catch(cl::sycl::exception const& e) {
            std::cout << "\t\tCaught synchronous SYCL exception during GEMMT:\n"
            << e.what() << std::endl << std::endl;
            rmg_error_handler (__FILE__, __LINE__, "Terminating");
        }
    }
#else

    if(typeid(DataType) == typeid(std::complex<double>)) {
        ztrtri(uplo, diag, &n, (std::complex<double> *)A, &n, info);
    }
    else {
        dtrtri(uplo, diag, &n, (double *)A, &n, info);

    }

#endif
}

