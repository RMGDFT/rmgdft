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

#define         dsyrk           RMG_FC_GLOBAL(dsyrk, DSYRK)
#define         zsyrk           RMG_FC_GLOBAL(zsyrk, ZSYRK)



#if SYCL_ENABLED
    #include <CL/sycl.hpp>
    #include "oneapi/mkl/blas.hpp"
    #include "mkl.h"
#else
extern "C" {
void dsyrk(const char *, const char *, int *, int *, double *, double *, int *, double *, double *, int *);
void zsyrk(const char *, const char *, int *, int *, std::complex<double> *, std::complex<double> *, int *, std::complex<double> *, std::complex<double> *, int *);
}
#endif

template void RmgSyrk<double>(char *, char *, int, int, double, double *, int,
                             double, double *, int);

template void RmgSyrk<std::complex<double> >(char *, char *, int, int, std::complex<double>, 
                    std::complex<double> *, int, std::complex<double>, std::complex<double> *, int);


template <typename DataType> void RmgSyrk(char *uplo, char *trans, int n, int k, 
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

    cublasStatus_t custat;
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
        custat = cublasXtZsyrk(ct.cublasxt_handle, fill_mode, cu_trans, n, k,
                            (cuDoubleComplex *)&alpha,
                            (cuDoubleComplex*)A, lda,
                            (cuDoubleComplex*)&beta, (cuDoubleComplex*)C, ldc );
        ProcessGpublasError(custat);
        RmgGpuError(__FILE__, __LINE__, custat, "Problem executing cublasXtZsyrkx");
        return;
    }
    if(ct.use_cublasxt && (typeid(DataType) == typeid(double)))
    {
        custat = cublasXtDsyrk(ct.cublasxt_handle, fill_mode, cu_trans, (size_t)n, (size_t)k,
                            (double*)&alpha,
                            (double*)A, (size_t)lda,
                            (double*)&beta, (double*)C, (size_t)ldc );
        ProcessGpublasError(custat);
        RmgGpuError(__FILE__, __LINE__, custat, "Problem executing cublasXtDsyrkx");
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

    DeviceSynchronize();
    if(typeid(DataType) == typeid(std::complex<double>)) {
        std::complex<double> *dA=(std::complex<double> *)A, *dC=(std::complex<double> *)C;
        if(!a_dev) gpuMalloc((void **)&dA, a_size * sizeof(std::complex<double>));
        if(!c_dev) gpuMalloc((void **)&dC, c_size * sizeof(std::complex<double>));
        if(!a_dev) cudaMemcpy(dA, A, a_size * sizeof(std::complex<double>), cudaMemcpyDefault);
        if(!c_dev && std::abs(beta) != 0.0) cudaMemcpy(dC, C, c_size * sizeof(std::complex<double>), cudaMemcpyDefault);
        custat = cublasZsyrk(ct.cublas_handle, fill_mode, cu_trans, n, k,
                            (cuDoubleComplex *)&alpha,
                            (cuDoubleComplex*)dA, lda,
                            (cuDoubleComplex*)&beta, (cuDoubleComplex*)dC, ldc );
        ProcessGpublasError(custat);
        RmgGpuError(__FILE__, __LINE__, custat, "Problem executing cublasZsyrkx");
        if(!c_dev) cudaMemcpy(C, dC, c_size * sizeof(std::complex<double>), cudaMemcpyDefault);
        if(!c_dev) gpuFree(dC);
        if(!a_dev) gpuFree(dA);
    }
    else {
        double *dA=(double *)A, *dC=(double *)C;
        if(!a_dev) gpuMalloc((void **)&dA, a_size * sizeof(double));
        if(!c_dev) gpuMalloc((void **)&dC, c_size * sizeof(double));
        if(!a_dev) cudaMemcpy(dA, A, a_size * sizeof(double), cudaMemcpyDefault);
        if(!c_dev && beta != 0.0) cudaMemcpy(dC, C, c_size * sizeof(double), cudaMemcpyDefault);
        custat = cublasDsyrk(ct.cublas_handle, fill_mode, cu_trans, n, k,
                            (double*)&alpha,
                            (double*)dA, lda,
                            (double*)&beta, (double*)dC, ldc );
        ProcessGpublasError(custat);
        RmgGpuError(__FILE__, __LINE__, custat, "Problem executing cublasDsyrkx");
        if(!c_dev) cudaMemcpy(C, dC, c_size * sizeof(double), cudaMemcpyDefault);
        if(!c_dev) gpuFree(dC);
        if(!a_dev) gpuFree(dA);
    }
    DeviceSynchronize();
    return;

#elif HIP_ENABLED
    hipblasStatus_t hipstat;
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
        if(!a_dev) gpuMalloc((void **)&dA, a_size * sizeof(std::complex<double>));
        if(!c_dev) gpuMalloc((void **)&dC, c_size * sizeof(std::complex<double>));
        if(!a_dev) hipMemcpyHtoD(dA, A, a_size * sizeof(std::complex<double>));
        if(!c_dev && std::abs(beta) != 0.0) hipMemcpyHtoD(dC, C, c_size * sizeof(std::complex<double>));
        hipstat = hipblasZsyrk(ct.hipblas_handle, fill_mode, hip_trans, n, k,
                            (hipblasDoubleComplex *)&alpha,
                            (hipblasDoubleComplex*)dA, lda,
                            (hipblasDoubleComplex*)&beta, (hipblasDoubleComplex*)dC, ldc );
        ProcessGpublasError(hipstat);
        RmgGpuError(__FILE__, __LINE__, hipstat, "Problem executing cublasZsyrkx");
        if(!c_dev) hipMemcpyDtoH(dC, C, c_size * sizeof(std::complex<double>));
        if(!c_dev) gpuFree(dC);
        if(!a_dev) gpuFree(dA);
    }
    else {
        double *dA=(double *)A, *dC=(double *)C;
        if(!a_dev) gpuMalloc((void **)&dA, a_size * sizeof(double));
        if(!c_dev) gpuMalloc((void **)&dC, c_size * sizeof(double));
        if(!a_dev) hipMemcpyHtoD(dA, A, a_size * sizeof(double));
        if(!c_dev && beta != 0.0) hipMemcpyHtoD(dC, C, c_size * sizeof(double));
        hipstat = hipblasDsyrk(ct.hipblas_handle, fill_mode, hip_trans, n, k,
                            (double*)&alpha,
                            (double*)dA, lda,
                            (double*)&beta, (double*)dC, ldc );
        ProcessGpublasError(hipstat);
        RmgGpuError(__FILE__, __LINE__, hipstat, "Problem executing hipblasDsyrkx");
        if(!c_dev) hipMemcpyDtoH(C, dC, c_size * sizeof(double));
        if(!c_dev) gpuFree(dC);
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
        zsyrk(uplo, trans, &n, &k, (std::complex<double> *)&alpha, (std::complex<double> *)A, &lda,
             (std::complex<double> *)&beta, (std::complex<double> *)C, &ldc);
    }
    else {
        dsyrk(uplo, trans, &n, &k, (double *)&alpha, (double *)A, &lda,
              (double *)&beta, (double *)C, &ldc);
    }

#endif
}

