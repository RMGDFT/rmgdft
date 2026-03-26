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


#if CUDA_ENABLED
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cublas_v2.h>
#endif


/*
  These functions are used to hide the details of the choleski decomposition and GPU
  interface from the higher level routines.

*/

#if SYCL_ENABLED
    #include <sycl/sycl.hpp>
    #include "oneapi/mkl/blas.hpp"
    #include <cstring>        // For strcmp
    #include <oneapi/mkl/lapack.hpp>
#endif

template void rmg::potrf<double>(char *, int, double *, int, int *);
template void rmg::potrf<std::complex<double>>(char *, int, std::complex<double> *, int, int *);


template <typename DataType> void rmg::potrf(char *uplo, int n, DataType *A, int lda, int *info)
{

#if CUDA_ENABLED
    int *dev_info, lwork;
    cusolverStatus_t custat;
    cublasFillMode_t fill_mode = CUBLAS_FILL_MODE_LOWER;

    if(!strcmp(uplo, "u")) fill_mode = CUBLAS_FILL_MODE_UPPER;
    if(!strcmp(uplo, "U")) fill_mode = CUBLAS_FILL_MODE_UPPER;

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

    rmg::sync_device();
    if(typeid(DataType) == typeid(std::complex<double>)) {
        std::complex<double> *dA=(std::complex<double> *)A, *work;
        if(!a_dev) gpuMalloc((void **)&dA, a_size * sizeof(std::complex<double>));
        if(!a_dev) cudaMemcpy(dA, A, a_size * sizeof(std::complex<double>), cudaMemcpyDefault);
        custat = cusolverDnZpotrf_bufferSize(ct.cusolver_handle, fill_mode, n,
			          (cuDoubleComplex *)dA, lda, &lwork);
        gpuMalloc((void **)&work, lwork * sizeof(std::complex<double>));
        custat = cusolverDnZpotrf(ct.cusolver_handle, fill_mode, n, (cuDoubleComplex *)dA, lda,
			          (cuDoubleComplex *)work, lwork, dev_info);
	if(custat != CUSOLVER_STATUS_SUCCESS)
            rmg::error(" cusolverDnZpotrf failed.");
	gpuFree(work);
        if(!a_dev) cudaMemcpy(A, dA, a_size * sizeof(std::complex<double>), cudaMemcpyDefault);
        if(!a_dev) gpuFree(dA);
    }
    else {
        double *dA=(double *)A, *work;
        if(!a_dev) gpuMalloc((void **)&dA, a_size * sizeof(double));
        if(!a_dev) cudaMemcpy(dA, A, a_size * sizeof(double), cudaMemcpyDefault);
        custat = cusolverDnDpotrf_bufferSize(ct.cusolver_handle, fill_mode, n, dA, lda, &lwork);
        gpuMalloc((void **)&work, lwork * sizeof(double));
        custat = cusolverDnDpotrf(ct.cusolver_handle, fill_mode, n, dA, lda, work, lwork, dev_info);
	if(custat != CUSOLVER_STATUS_SUCCESS)
            rmg::error(" cusolverDnDpotrf failed.");
	gpuFree(work);
        if(!a_dev) cudaMemcpy(A, dA, a_size * sizeof(double), cudaMemcpyDefault);
        if(!a_dev) gpuFree(dA);
    }
    cudaMemcpy(info, dev_info, sizeof(int), cudaMemcpyDefault);
    rmg::sync_device();
    gpuFree(dev_info);
    return;

#elif HIP_ENABLED
    rocblas_status rocstat;
    rocblas_fill fill_mode = rocblas_fill_lower;

    if(!strcmp(uplo, "u")) fill_mode = rocblas_fill_upper;
    if(!strcmp(uplo, "U")) fill_mode = rocblas_fill_upper;

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
        if(!a_dev) rmg::error(hipMemcpyHtoD(dA, A, a_size * sizeof(std::complex<double>)));
        rocstat = rocsolver_zpotrf(ct.roc_handle, fill_mode, n, (rocblas_double_complex *)dA, lda, dev_info);
        if (rocstat != rocblas_status_success) 
            rmg::error("Problem executing rocsolver_zpotrf");
        if(!a_dev) rmg::error(hipMemcpyDtoH(A, dA, a_size * sizeof(std::complex<double>)));
        if(!a_dev) gpuFree(dA);
    }
    else {
        double *dA=(double *)A;
        if(!a_dev) gpuMalloc((void **)&dA, a_size * sizeof(double));
        if(!a_dev) rmg::error(hipMemcpyHtoD(dA, A, a_size * sizeof(double)));
        rocstat = rocsolver_dpotrf(ct.roc_handle, fill_mode, n, dA, lda, dev_info);
        if (rocstat != rocblas_status_success) 
            rmg::error("Problem executing rocsolver_dpotrf");
        if(!a_dev) rmg::error(hipMemcpyDtoH(A, dA, a_size * sizeof(double)));
        if(!a_dev) gpuFree(dA);
    }
    rmg::error(hipMemcpyDtoH(info, dev_info, sizeof(int)));
    gpuFree(dev_info);
#elif SYCL_ENABLED

    // Determine upper/lower fill mode
    oneapi::mkl::uplo fill_mode = oneapi::mkl::uplo::lower;
    if (!strcmp(uplo, "u")) fill_mode = oneapi::mkl::uplo::upper;
    if (!strcmp(uplo, "U")) fill_mode = oneapi::mkl::uplo::upper;

    // Allocate scratchpad memory with the correct type
    std::int64_t scratchpad_size = oneapi::mkl::lapack::potrf_scratchpad_size<DataType>(ct.sycl_Q, fill_mode, n, lda);
    DataType *scratchpad = sycl::malloc_device<DataType>(scratchpad_size, ct.sycl_Q);

    // POTRF: Cholesky factorization
    // Computes L or U such that A = L * L^T or A = U^T * U
    try {
        // Use USM pointer-based API; returns a sycl::event
        sycl::event potrf_event = oneapi::mkl::lapack::potrf(
            ct.sycl_Q, fill_mode, n, A, lda, scratchpad, scratchpad_size);
        // Wait for completion
        potrf_event.wait();
    }
    catch (sycl::exception const& e) {
        std::cout << "\t\tCaught synchronous SYCL exception during POTRF:\n"
                  << e.what() << std::endl << std::endl;
        rmg::error("Terminating");
    }

    // Free scratchpad memory
    sycl::free(scratchpad, ct.sycl_Q);


#else

    if(typeid(DataType) == typeid(std::complex<double>)) {
        zpotrf(uplo, &n, (std::complex<double> *)A, &n, info);
    }
    else {
        dpotrf(uplo, &n, (double *)A, &n, info);
    }

#endif
}

