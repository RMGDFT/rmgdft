#include <complex>
#include <typeinfo>
#include <string.h>

#include "const.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "rmg_gemm.h"
#include "GpuAlloc.h"
#include "rmg_dev_allocate.h"

#include "RmgTimer.h"
#include "transition.h"


#if CUDA_ENABLED
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cublas_v2.h>
#endif

#define         dgemm           RMG_FC_GLOBAL(dgemm, DGEMM)
#define         zgemm           RMG_FC_GLOBAL(zgemm, ZGEMM)

#if SYCL_ENABLED
    #include <CL/sycl.hpp>
    #include "oneapi/mkl/blas.hpp"
    #include "mkl.h"
#else

extern "C" {
void dgemm(const char *, const char *, int *, int *, int *, double *, double *, int *, double *, int *, double *, double *, int *);
void zgemm(const char *, const char *, int *, int *, int *, std::complex<double> *, std::complex<double> *, int *, std::complex<double> *, int *, std::complex<double> *, std::complex<double> *, int *);
}
#endif

/*
  These functions are used to hide the details of the matrix multiplication data types and GPU 
  utilization from the higher level routines.

  The first 13 arguments are the same as the standard dgemm args but with scalar quantities passed
  by value instead of by reference.

*/


template void rmg::gemm_strided_batched<double>(char *, char *, int, int, int, double, double *, int, size_t, double *, int,  size_t,
                                  double, double *, int,  size_t, int);

template void rmg::gemm_strided_batched<std::complex<double> >(char *, char *, int, int, int, std::complex<double>, 
                      std::complex<double> *, int,  size_t, std::complex<double> *, int,  size_t,
                      std::complex<double>, std::complex<double> *, int,  size_t, int);


template <typename DataType> void rmg::gemm_strided_batched(char *transa, char *transb, int m, int n, int k, 
                             DataType alpha, DataType *A, int lda,  size_t strideA, DataType *B, int ldb,  size_t strideB, DataType beta, 
                             DataType *C, int ldc,  size_t strideC, int batchCount)
{

#if BLAS_PROFILE
    if(typeid(DataType) == typeid(std::complex<double>))
    {
        if(pct.gridpe==0) printf("ZGEMM CALL m=%d n=%d k=%d, batch=%d\n",m,n,k, batchCount);
    }
    else
    {
        if(pct.gridpe==0) printf("DGEMM CALL m=%d n=%d k=%d batch=%d\n",m,n,k, batchCount);
    }
#endif

#if CUDA_ENABLED

    cublasOperation_t cu_transA = CUBLAS_OP_N, cu_transB = CUBLAS_OP_N;

    if(!strcmp(transa, "t")) cu_transA = CUBLAS_OP_T;
    if(!strcmp(transa, "T")) cu_transA = CUBLAS_OP_T;
    if(!strcmp(transa, "c")) cu_transA = CUBLAS_OP_C;
    if(!strcmp(transa, "C")) cu_transA = CUBLAS_OP_C;

    if(!strcmp(transb, "t")) cu_transB = CUBLAS_OP_T;
    if(!strcmp(transb, "T")) cu_transB = CUBLAS_OP_T;
    if(!strcmp(transb, "c")) cu_transB = CUBLAS_OP_C;
    if(!strcmp(transb, "C")) cu_transB = CUBLAS_OP_C;

    int ka = m;
    if(!strcmp("n", transa)) ka = k;
    if(!strcmp("N", transa)) ka = k;

    int kb = k;
    if(!strcmp("n", transb)) kb = n;
    if(!strcmp("N", transb)) kb = n;

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
    if(cudaerr == cudaSuccess && attr.type == cudaMemoryTypeDevice) a_dev = true;
    cudaerr = cudaPointerGetAttributes(&attr, B);
    bool b_dev = false;
    if(cudaerr == cudaSuccess && attr.type == cudaMemoryTypeDevice) b_dev = true;
    cudaerr = cudaPointerGetAttributes(&attr, C);
    bool c_dev = false;
    if(cudaerr == cudaSuccess && attr.type == cudaMemoryTypeDevice) c_dev = true;
#endif

    size_t a_size = (size_t)lda * (size_t)ka * (size_t)batchCount;
    size_t b_size = (size_t)ldb * (size_t)kb * (size_t)batchCount;
    size_t c_size = (size_t)ldc * (size_t)n * (size_t)batchCount;

    rmg::sync_device();
    if(typeid(DataType) == typeid(std::complex<double>)) {
        std::complex<double> *dA=(std::complex<double> *)A, *dB=(std::complex<double> *)B, *dC=(std::complex<double> *)C;
        if(!a_dev) rmg_device_pool->malloc(&dA, a_size);
        if(!b_dev) rmg_device_pool->malloc(&dB, b_size);
        if(!c_dev) rmg_device_pool->malloc(&dC, c_size);
        if(!a_dev) gpuMemcpy(dA, A, a_size * sizeof(std::complex<double>), gpuMemcpyDefault);
        if(!b_dev) gpuMemcpy(dB, B, b_size * sizeof(std::complex<double>), gpuMemcpyDefault);
        if(!c_dev && std::abs(beta) != 0.0) gpuMemcpy(dC, C, c_size * sizeof(std::complex<double>), gpuMemcpyDefault);
        rmg::error(cublasZgemmStridedBatched(ct.cublas_handle, cu_transA, cu_transB, m, n, k,
                            (cuDoubleComplex *)&alpha,
                            (cuDoubleComplex*)dA, lda, strideA,
                            (cuDoubleComplex*)dB, ldb, strideB,
                            (cuDoubleComplex*)&beta, (cuDoubleComplex*)dC, ldc, strideC, batchCount ));
        if(!c_dev) gpuMemcpy(C, dC, c_size * sizeof(std::complex<double>), gpuMemcpyDefault);
        if(!c_dev) rmg_device_pool->free(dC);
        if(!b_dev) rmg_device_pool->free(dB);
        if(!a_dev) rmg_device_pool->free(dA);
    }
    else {
        double *dA=(double *)A, *dB=(double *)B, *dC=(double *)C;
        if(!a_dev) rmg_device_pool->malloc(&dA, a_size);
        if(!b_dev) rmg_device_pool->malloc(&dB, b_size);
        if(!c_dev) rmg_device_pool->malloc(&dC, c_size);
        if(!a_dev) gpuMemcpy(dA, A, a_size * sizeof(double), gpuMemcpyDefault);
        if(!b_dev) gpuMemcpy(dB, B, b_size * sizeof(double), gpuMemcpyDefault);
        if(!c_dev && beta != 0.0) gpuMemcpy(dC, C, c_size * sizeof(double), gpuMemcpyDefault);
        rmg::error(cublasDgemmStridedBatched(ct.cublas_handle, cu_transA, cu_transB, m, n, k,
                            (double*)&alpha,
                            (double*)dA, lda, strideA,
                            (double*)dB, ldb, strideB,
                            (double*)&beta, (double*)dC, ldc, strideC, batchCount ));
        if(!c_dev) gpuMemcpy(C, dC, c_size * sizeof(double), gpuMemcpyDefault);
        if(!c_dev) rmg_device_pool->free(dC);
        if(!b_dev) rmg_device_pool->free(dB);
        if(!a_dev) rmg_device_pool->free(dA);
    }
    rmg::sync_device();
    return;

#elif HIP_ENABLED
    hipPointerAttribute_t attr;
    hipError_t hiperr;
    hiperr = hipPointerGetAttributes(&attr, A);
    bool a_dev = false;
    if(hiperr == hipSuccess && attr.type == hipMemoryTypeDevice) a_dev = true;
    hiperr = hipPointerGetAttributes(&attr, B);
    bool b_dev = false;
    if(hiperr == hipSuccess && attr.type == hipMemoryTypeDevice) b_dev = true;
    hiperr = hipPointerGetAttributes(&attr, C);
    bool c_dev = false;
    if(hiperr == hipSuccess && attr.type == hipMemoryTypeDevice) c_dev = true;
    hipblasOperation_t hip_transA = HIPBLAS_OP_N, hip_transB = HIPBLAS_OP_N;

    if(!strcmp(transa, "t")) hip_transA = HIPBLAS_OP_T;
    if(!strcmp(transa, "T")) hip_transA = HIPBLAS_OP_T;
    if(!strcmp(transa, "c")) hip_transA = HIPBLAS_OP_C;
    if(!strcmp(transa, "C")) hip_transA = HIPBLAS_OP_C;

    if(!strcmp(transb, "t")) hip_transB = HIPBLAS_OP_T;
    if(!strcmp(transb, "T")) hip_transB = HIPBLAS_OP_T;
    if(!strcmp(transb, "c")) hip_transB = HIPBLAS_OP_C;
    if(!strcmp(transb, "C")) hip_transB = HIPBLAS_OP_C;

    int ka = m;
    if(!strcmp("n", transa)) ka = k;
    if(!strcmp("N", transa)) ka = k;

    int kb = k;
    if(!strcmp("n", transb)) kb = n;
    if(!strcmp("N", transb)) kb = n;

    size_t a_size = (size_t)lda * (size_t)ka * (size_t)batchCount;
    size_t b_size = (size_t)ldb * (size_t)kb * (size_t)batchCount;
    size_t c_size = (size_t)ldc * (size_t)n * (size_t)batchCount;

    if(typeid(DataType) == typeid(std::complex<double>)) {
        std::complex<double> *dA=(std::complex<double> *)A, *dB=(std::complex<double> *)B, *dC=(std::complex<double> *)C;
        if(!a_dev) rmg_device_pool->malloc(&dA, a_size);
        if(!b_dev) rmg_device_pool->malloc(&dB, b_size);
        if(!c_dev) rmg_device_pool->malloc(&dC, c_size);
        if(!a_dev) rmg::error(hipMemcpyHtoD(dA, A, a_size * sizeof(std::complex<double>)));
        if(!b_dev) rmg::error(hipMemcpyHtoD(dB, B, b_size * sizeof(std::complex<double>)));
        if(!c_dev && std::abs(beta) != 0.0) rmg::error(hipMemcpyHtoD(dC, C, c_size * sizeof(std::complex<double>)));
        rmg::error(hipblasZgemmStridedBatched(ct.hipblas_handle, hip_transA, hip_transB, m, n, k,
                            (hipDoubleComplex *)&alpha,
                            (hipDoubleComplex*)dA, lda, strideA,
                            (hipDoubleComplex*)dB, ldb, strideB,
                            (hipDoubleComplex*)&beta, (hipDoubleComplex*)dC, ldc, strideC, batchCount ));
        if(!c_dev) rmg::error(hipMemcpyDtoH(C, dC, c_size * sizeof(std::complex<double>)));
        if(!c_dev) rmg_device_pool->free(dC);
        if(!b_dev) rmg_device_pool->free(dB);
        if(!a_dev) rmg_device_pool->free(dA);
    }
    else {
        double *dA=(double *)A, *dB=(double *)B, *dC=(double *)C;
        if(!a_dev) rmg_device_pool->malloc(&dA, a_size);
        if(!b_dev) rmg_device_pool->malloc(&dB, b_size);
        if(!c_dev) rmg_device_pool->malloc(&dC, c_size);
        if(!a_dev) rmg::error(hipMemcpyHtoD(dA, A, a_size * sizeof(double)));
        if(!b_dev) rmg::error(hipMemcpyHtoD(dB, B, b_size * sizeof(double)));
        if(!c_dev && beta != 0.0) rmg::error(hipMemcpyHtoD(dC, C, c_size * sizeof(double)));
        rmg::error(hipblasDgemmStridedBatched(ct.hipblas_handle, hip_transA, hip_transB, m, n, k,
                            (double*)&alpha,
                            (double*)dA, lda, strideA,
                            (double*)dB, ldb, strideB,
                            (double*)&beta, (double*)dC, ldc, strideC,batchCount));
        if(!c_dev) rmg::error(hipMemcpyDtoH(C, dC, c_size * sizeof(double)));
        if(!c_dev) rmg_device_pool->free(dC);
        if(!b_dev) rmg_device_pool->free(dB);
        if(!a_dev) rmg_device_pool->free(dA);
    }

#elif SYCL_ENABLED

    oneapi::mkl::transpose sycl_transA = oneapi::mkl::transpose::nontrans;
    oneapi::mkl::transpose sycl_transB = oneapi::mkl::transpose::nontrans;

    if(!strcmp(transa, "t")) sycl_transA = oneapi::mkl::transpose::trans;
    if(!strcmp(transa, "T")) sycl_transA = oneapi::mkl::transpose::trans;
    if(!strcmp(transa, "c")) sycl_transA = oneapi::mkl::transpose::conjtrans;
    if(!strcmp(transa, "C")) sycl_transA = oneapi::mkl::transpose::conjtrans;

    if(!strcmp(transb, "t")) sycl_transB = oneapi::mkl::transpose::trans;
    if(!strcmp(transb, "T")) sycl_transB = oneapi::mkl::transpose::trans;
    if(!strcmp(transb, "c")) sycl_transB = oneapi::mkl::transpose::conjtrans;
    if(!strcmp(transb, "C")) sycl_transB = oneapi::mkl::transpose::conjtrans;

    int ka = m;
    if(!strcmp("n", transa)) ka = k;
    if(!strcmp("N", transa)) ka = k;

    int kb = k;
    if(!strcmp("n", transb)) kb = n;
    if(!strcmp("N", transb)) kb = n;

    size_t a_size = (size_t)lda * (size_t)ka * (size_t)batchCount;
    size_t b_size = (size_t)ldb * (size_t)kb * (size_t)batchCount;
    size_t c_size = (size_t)ldc * (size_t)n * (size_t)batchCount;


    cl::sycl::buffer<DataType, 1> bufA((DataType *)A, a_size, {cl::sycl::property::buffer::use_host_ptr()});
    cl::sycl::buffer<DataType, 1> bufB((DataType *)B, b_size, {cl::sycl::property::buffer::use_host_ptr()});
    cl::sycl::buffer<DataType, 1> bufC((DataType *)C, c_size, {cl::sycl::property::buffer::use_host_ptr()});
    try {
        oneapi::mkl::blas::gemm_batch(ct.sycl_Q, sycl_transA, sycl_transB, m, n, k, alpha,
                                bufA, lda, strideA, bufB, ldb, strideB,
                                beta, bufC, ldc, strideC, batchCount);
    }
    catch(cl::sycl::exception const& e) {
        std::cout << "\t\tCaught synchronous SYCL exception during GEMM:\n"
        << e.what() << std::endl << std::endl;
    }

#else

    
//    RmgTimer *RT = new RmgTimer("gemmmmm ");
    for(int p = 0; p < batchCount; p++)
    {
        if(typeid(DataType) == typeid(std::complex<double>))
        {
            if(ct.use_alt_zgemm)
                MyZgemm(transa, transb, m, n, k, (std::complex<double> *)(&alpha), (std::complex<double> *)&A[p*strideA], lda, 
                        (std::complex<double> *)&B[p*strideB], ldb, (std::complex<double> *)(&beta), (std::complex<double> *)&C[p*strideC], ldc);
            else
                zgemm(transa, transb, &m, &n, &k, (std::complex<double> *)&alpha, (std::complex<double> *)&A[p*strideA], &lda,
                        (std::complex<double> *)&B[p*strideB], &ldb, (std::complex<double> *)&beta, (std::complex<double> *)&C[p*strideC], &ldc);
        }
        else {
            dgemm(transa, transb, &m, &n, &k, (double *)&alpha, (double *)&A[p*strideA], &lda, 
                    (double *)&B[p*strideB], &ldb, (double *)&beta, (double *)&C[p*strideC], &ldc);
        }
    }
    //    delete RT;

#endif
}
