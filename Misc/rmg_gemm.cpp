/*  
 *  
 * Copyright 2026 The RMG Project Developers. See the COPYRIGHT file 
 * at the top-level directory of this distribution or in the current
 * directory.
 *      
 * This file is part of RMG. 
 * RMG is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * any later version.
 *
 * RMG is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *      
 * You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
*/

#if HIP_ENABLED
#include "tiled_mm.hpp"
#endif
#include <complex>
#include <typeinfo>
#include <string.h>

#include "const.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "rmg_gemm.h"
#include "rmg_dev_allocate.h"
#include "GpuAlloc.h"

#include "RmgTimer.h"
#include "transition.h"
#include "rmg_error.h"


#if CUDA_ENABLED
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cublas_v2.h>
#endif


#define         dgemm           RMG_FC_GLOBAL(dgemm, DGEMM)
#define         zgemm           RMG_FC_GLOBAL(zgemm, ZGEMM)
#define         cgemm           RMG_FC_GLOBAL(cgemm, CGEMM)
#define         sgemm           RMG_FC_GLOBAL(sgemm, SGEMM)


#if SYCL_ENABLED
    #include <sycl/sycl.hpp>
    #include "oneapi/mkl/blas.hpp"
	template <typename T>
	class DeviceTemp
	{
	public:
	    DeviceTemp(sycl::queue& q, size_t count)
		: queue(q)
	    {
		if (count > 0)
		    ptr = sycl::malloc_device<T>(count, queue);
	    }

	    ~DeviceTemp()
	    {
		if (ptr)
		    sycl::free(ptr, queue);
	    }
	    
	    T* get() const { return ptr; }
	    
	    DeviceTemp(const DeviceTemp&) = delete;
	    DeviceTemp& operator=(const DeviceTemp&) = delete;

	private:
	    sycl::queue& queue;
	    T* ptr = nullptr;
	};
	// Determine whether a pointer is accessible to the SYCL device.
	template <typename T>
	bool is_device_accessible(T* ptr, sycl::queue& q)
	{
	    auto type = sycl::get_pointer_type(ptr, q.get_context());

	    return type == sycl::usm::alloc::device ||
		   type == sycl::usm::alloc::shared;
	}
// Generic GEMM implementation.
template <typename DataType>
void sycl_gemm(
    sycl::queue& q,
    const char* transa,
    const char* transb,
    int m, int n, int k,
    DataType alpha,
    DataType* A, int lda,
    DataType* B, int ldb,
    DataType beta,
    DataType* C, int ldc)
{
    // Determine transpose operations.
    auto get_transpose = [](const char* trans)
    {
        char c = static_cast<char>(
            std::tolower(static_cast<unsigned char>(trans[0])));

        switch (c)
        {
            case 't':
                return oneapi::mkl::transpose::trans;

            case 'c':
                return oneapi::mkl::transpose::conjtrans;

            default:
                return oneapi::mkl::transpose::nontrans;
        }
    };

    auto transA = get_transpose(transa);
    auto transB = get_transpose(transb);

    // Determine matrix storage sizes (column-major).
    size_t ka = (transA == oneapi::mkl::transpose::nontrans) ? k : m;
    size_t kb = (transB == oneapi::mkl::transpose::nontrans) ? n : k;

    size_t a_size = static_cast<size_t>(lda) * ka;
    size_t b_size = static_cast<size_t>(ldb) * kb;
    size_t c_size = static_cast<size_t>(ldc) * n;

    // Determine whether pointers are already device-accessible.
    bool a_dev = is_device_accessible(A, q);
    bool b_dev = is_device_accessible(B, q);
    bool c_dev = is_device_accessible(C, q);

    // Allocate temporary device storage where necessary.
    DeviceTemp<DataType> tempA(q, a_dev ? 0 : a_size);
    DeviceTemp<DataType> tempB(q, b_dev ? 0 : b_size);
    DeviceTemp<DataType> tempC(q, c_dev ? 0 : c_size);

    DataType* dA = a_dev ? A : tempA.get();
    DataType* dB = b_dev ? B : tempB.get();
    DataType* dC = c_dev ? C : tempC.get();

    // Copy input matrices to device.
    if (!a_dev)
        q.memcpy(dA, A, a_size * sizeof(DataType));

    if (!b_dev)
        q.memcpy(dB, B, b_size * sizeof(DataType));

    // C only needs initialization when beta is nonzero.
    if (!c_dev && beta != DataType{})
        q.memcpy(dC, C, c_size * sizeof(DataType));

    // Ensure input transfers are complete before GEMM.
    q.wait_and_throw();

    // Execute oneMKL GEMM.
    oneapi::mkl::blas::column_major::gemm(
        q,
        transA,
        transB,
        m, n, k,
        alpha,
        dA, lda,
        dB, ldb,
        beta,
        dC, ldc
    );

    // Ensure GEMM is complete before copying results or freeing memory.
    q.wait_and_throw();

    // Copy result back to host if C was not device-accessible.
    if (!c_dev)
        q.memcpy(C, dC, c_size * sizeof(DataType)).wait_and_throw();
}

#else
extern "C" {
void dgemm(const char *, const char *, int *, int *, int *, double *, double *, int *, double *, int *, double *, double *, int *);
void sgemm(const char *, const char *, int *, int *, int *, float *, float *, int *, float *, int *, float *, float *, int *);
void zgemm(const char *, const char *, int *, int *, int *, std::complex<double> *, std::complex<double> *, int *, std::complex<double> *, int *, std::complex<double> *, std::complex<double> *, int *);
void cgemm(const char *, const char *, int *, int *, int *, std::complex<float> *, std::complex<float> *, int *, std::complex<float> *, int *, std::complex<float> *, std::complex<float> *, int *);
}
#endif


/*
  These functions are used to hide the details of the matrix multiplication data types and GPU 
  utilization from the higher level routines.

  The first 13 arguments are the same as the standard dgemm args but with scalar quantities passed
  by value instead of by reference.

*/


template void rmg::gemm<double>(char *, char *, int, int, int, double, double *, int, double *, int, 
                                  double, double *, int);

template void rmg::gemm<float>(char *, char *, int, int, int, float, float *, int, float *, int, 
                                  float, float *, int);

template void rmg::gemm<std::complex<double> >(char *, char *, int, int, int, std::complex<double>, 
                      std::complex<double> *, int, std::complex<double> *, int, 
                      std::complex<double>, std::complex<double> *, int);

template void rmg::gemm<std::complex<float> >(char *, char *, int, int, int, std::complex<float>, 
                      std::complex<float> *, int, std::complex<float> *, int, 
                      std::complex<float>, std::complex<float> *, int);


template <typename DataType> void rmg::gemm(char *transa, char *transb, int m, int n, int k, 
                             DataType alpha, DataType *A, int lda, DataType *B, int ldb, DataType beta, 
                             DataType *C, int ldc)
{

#if BLAS_PROFILE
    if(typeid(DataType) == typeid(std::complex<double>))
    {
        if(pct.gridpe==0) printf("ZGEMM CALL m=%d n=%d k=%d\n",m,n,k);
    }
    else
    {
        if(pct.gridpe==0) printf("DGEMM CALL m=%d n=%d k=%d\n",m,n,k);
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

    if(ct.use_cublasxt && (typeid(DataType) == typeid(std::complex<double>)))
    {
        rmg::error(cublasXtZgemm(ct.cublasxt_handle, cu_transA, cu_transB, m, n, k,
                            (cuDoubleComplex *)&alpha,
                            (cuDoubleComplex*)A, lda,
                            (cuDoubleComplex*)B, ldb,
                            (cuDoubleComplex*)&beta, (cuDoubleComplex*)C, ldc ));
        return;
    }
    if(ct.use_cublasxt && (typeid(DataType) == typeid(std::complex<float>)))
    {
        rmg::error(cublasXtCgemm(ct.cublasxt_handle, cu_transA, cu_transB, m, n, k,
                            (cuComplex *)&alpha,
                            (cuComplex*)A, lda,
                            (cuComplex*)B, ldb,
                            (cuComplex*)&beta, (cuComplex*)C, ldc ));
        return;
    }
    if(ct.use_cublasxt && (typeid(DataType) == typeid(float)))
    {
        rmg::error(cublasXtSgemm(ct.cublasxt_handle, cu_transA, cu_transB, m, n, k,
                            (float*)&alpha,
                            (float*)A, lda,
                            (float*)B, ldb,
                            (float*)&beta, (float*)C, ldc ));
        return;
    }
    if(ct.use_cublasxt && (typeid(DataType) == typeid(double)))
    {
        rmg::error(cublasXtDgemm(ct.cublasxt_handle, cu_transA, cu_transB, m, n, k,
                            (double*)&alpha,
                            (double*)A, lda,
                            (double*)B, ldb,
                            (double*)&beta, (double*)C, ldc ));
        return;
    }

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

    size_t a_size = (size_t)lda * (size_t)ka;
    size_t b_size = (size_t)ldb * (size_t)kb;
    size_t c_size = (size_t)ldc * (size_t)n;

    rmg::sync_device();
    if(typeid(DataType) == typeid(std::complex<double>)) {
        std::complex<double> *dA=(std::complex<double> *)A, *dB=(std::complex<double> *)B, *dC=(std::complex<double> *)C;
        if(!a_dev) rmg_device_pool->malloc(&dA, a_size);
        if(!b_dev) rmg_device_pool->malloc(&dB, b_size);
        if(!c_dev) rmg_device_pool->malloc(&dC, c_size);
        if(!a_dev) gpuMemcpy(dA, A, a_size * sizeof(std::complex<double>), gpuMemcpyDefault);
        if(!b_dev) gpuMemcpy(dB, B, b_size * sizeof(std::complex<double>), gpuMemcpyDefault);
        if(!c_dev && std::abs(beta) != 0.0) gpuMemcpy(dC, C, c_size * sizeof(std::complex<double>), gpuMemcpyDefault);
        rmg::error(cublasZgemm(ct.cublas_handle, cu_transA, cu_transB, m, n, k,
                            (cuDoubleComplex *)&alpha,
                            (cuDoubleComplex*)dA, lda,
                            (cuDoubleComplex*)dB, ldb,
                            (cuDoubleComplex*)&beta, (cuDoubleComplex*)dC, ldc ));
        if(!c_dev) gpuMemcpy(C, dC, c_size * sizeof(std::complex<double>), gpuMemcpyDefault);
        if(!c_dev) rmg_device_pool->free(dC);
        if(!b_dev) rmg_device_pool->free(dB);
        if(!a_dev) rmg_device_pool->free(dA);
    }
    else if(typeid(DataType) == typeid(std::complex<float>)) {
        std::complex<float> *dA=(std::complex<float> *)A, *dB=(std::complex<float> *)B, *dC=(std::complex<float> *)C;
        if(!a_dev) rmg_device_pool->malloc(&dA, a_size);
        if(!b_dev) rmg_device_pool->malloc(&dB, b_size);
        if(!c_dev) rmg_device_pool->malloc(&dC, c_size);
        if(!a_dev) gpuMemcpy(dA, A, a_size * sizeof(std::complex<float>), gpuMemcpyDefault);
        if(!b_dev) gpuMemcpy(dB, B, b_size * sizeof(std::complex<float>), gpuMemcpyDefault);
        if(!c_dev && std::abs(beta) != 0.0) gpuMemcpy(dC, C, c_size * sizeof(std::complex<float>), gpuMemcpyDefault);
        rmg::error(cublasCgemm(ct.cublas_handle, cu_transA, cu_transB, m, n, k,
                            (cuComplex *)&alpha,
                            (cuComplex*)dA, lda,
                            (cuComplex*)dB, ldb,
                            (cuComplex*)&beta, (cuComplex*)dC, ldc ));
        if(!c_dev) gpuMemcpy(C, dC, c_size * sizeof(std::complex<float>), gpuMemcpyDefault);
        if(!c_dev) rmg_device_pool->free(dC);
        if(!b_dev) rmg_device_pool->free(dB);
        if(!a_dev) rmg_device_pool->free(dA);
    }
    else if(typeid(DataType) == typeid(float)) {
        float *dA=(float *)A, *dB=(float *)B, *dC=(float *)C;
        if(!a_dev) rmg_device_pool->malloc(&dA, a_size);
        if(!b_dev) rmg_device_pool->malloc(&dB, b_size);
        if(!c_dev) rmg_device_pool->malloc(&dC, c_size);
        if(!a_dev) gpuMemcpy(dA, A, a_size * sizeof(float), gpuMemcpyDefault);
        if(!b_dev) gpuMemcpy(dB, B, b_size * sizeof(float), gpuMemcpyDefault);
        if(!c_dev && std::abs(beta) != 0.0) gpuMemcpy(dC, C, c_size * sizeof(float), gpuMemcpyDefault);
        rmg::error(cublasSgemm(ct.cublas_handle, cu_transA, cu_transB, m, n, k,
                            (float*)&alpha,
                            (float*)dA, lda,
                            (float*)dB, ldb,
                            (float*)&beta, (float*)dC, ldc ));
        if(!c_dev) gpuMemcpy(C, dC, c_size * sizeof(float), gpuMemcpyDefault);
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
        if(!c_dev && std::abs(beta) != 0.0) gpuMemcpy(dC, C, c_size * sizeof(double), gpuMemcpyDefault);
        rmg::error(cublasDgemm(ct.cublas_handle, cu_transA, cu_transB, m, n, k,
                            (double*)&alpha,
                            (double*)dA, lda,
                            (double*)dB, ldb,
                            (double*)&beta, (double*)dC, ldc ));
        if(!c_dev) gpuMemcpy(C, dC, c_size * sizeof(double), gpuMemcpyDefault);
        if(!c_dev) rmg_device_pool->free(dC);
        if(!b_dev) rmg_device_pool->free(dB);
        if(!a_dev) rmg_device_pool->free(dA);

    }
    rmg::sync_device();
    return;

#elif HIP_ENABLED
    // For Tiled-MM
    if(ct.use_cublasxt)
    {
       static std::unique_ptr<gpu::mm_handle<DataType>> ctx;

       char trans_a = 'N', trans_b = 'N';
       if(!strcmp(transa, "n")) trans_a = 'N';
       if(!strcmp(transb, "n")) trans_b = 'N';
       if(!strcmp(transa, "t")) trans_a = 'T';
       if(!strcmp(transb, "t")) trans_b = 'T';
       if(!strcmp(transa, "c")) trans_a = 'C';
       if(!strcmp(transb, "c")) trans_b = 'C';
       if(!strcmp(transa, "C")) trans_a = 'C';
       if(!strcmp(transb, "C")) trans_b = 'C';

       // can probably tune these
       if(!ctx) ctx = gpu::make_context<DataType>(2, 5000, 5000, 5000);
       gpu::gemm<DataType>(*ctx,
          trans_a, trans_b,
          m, n, k,
          alpha,
          A, lda, 
          B, ldb, 
          beta,
          C, ldc, false, true);
          return;
    }

    //hipDeviceSynchronize();
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

    size_t a_size = (size_t)lda * (size_t)ka;
    size_t b_size = (size_t)ldb * (size_t)kb;
    size_t c_size = (size_t)ldc * (size_t)n;

    if(typeid(DataType) == typeid(std::complex<double>)) {
        std::complex<double> *dA=(std::complex<double> *)A, *dB=(std::complex<double> *)B, *dC=(std::complex<double> *)C;
        if(!a_dev) rmg_device_pool->malloc(&dA, a_size);
        if(!b_dev) rmg_device_pool->malloc(&dB, b_size);
        if(!c_dev) rmg_device_pool->malloc(&dC, c_size);
        if(!a_dev) rmg::error(hipMemcpyHtoD(dA, A, a_size * sizeof(std::complex<double>)));
        if(!b_dev) rmg::error(hipMemcpyHtoD(dB, B, b_size * sizeof(std::complex<double>)));
        if(!c_dev && std::abs(beta) != 0.0) rmg::error(hipMemcpyHtoD(dC, C, c_size * sizeof(std::complex<double>)));
        rmg::error(hipblasZgemm(ct.hipblas_handle, hip_transA, hip_transB, m, n, k,
                            (hipDoubleComplex *)&alpha,
                            (hipDoubleComplex*)dA, lda,
                            (hipDoubleComplex*)dB, ldb,
                            (hipDoubleComplex*)&beta, (hipDoubleComplex*)dC, ldc ));
        if(!c_dev) rmg::error(hipMemcpyDtoH(C, dC, c_size * sizeof(std::complex<double>)));
        if(!c_dev) rmg_device_pool->free(dC);
        if(!b_dev) rmg_device_pool->free(dB);
        if(!a_dev) rmg_device_pool->free(dA);
    }
    else if(typeid(DataType) == typeid(std::complex<float>)) {
        std::complex<float> *dA=(std::complex<float> *)A, *dB=(std::complex<float> *)B, *dC=(std::complex<float> *)C;
        if(!a_dev) rmg_device_pool->malloc(&dA, a_size);
        if(!b_dev) rmg_device_pool->malloc(&dB, b_size);
        if(!c_dev) rmg_device_pool->malloc(&dC, c_size);
        if(!a_dev) rmg::error(hipMemcpyHtoD(dA, A, a_size * sizeof(std::complex<float>)));
        if(!b_dev) rmg::error(hipMemcpyHtoD(dB, B, b_size * sizeof(std::complex<float>)));
        if(!c_dev && std::abs(beta) != 0.0) rmg::error(hipMemcpyHtoD(dC, C, c_size * sizeof(std::complex<float>)));
        rmg::error(hipblasCgemm(ct.hipblas_handle, hip_transA, hip_transB, m, n, k,
                            (hipFloatComplex *)&alpha,
                            (hipFloatComplex*)dA, lda,
                            (hipFloatComplex*)dB, ldb,
                            (hipFloatComplex*)&beta, (hipFloatComplex*)dC, ldc ));
        if(!c_dev) rmg::error(hipMemcpyDtoH(C, dC, c_size * sizeof(std::complex<float>)));
        if(!c_dev) rmg_device_pool->free(dC);
        if(!b_dev) rmg_device_pool->free(dB);
        if(!a_dev) rmg_device_pool->free(dA);
    }
    else if(typeid(DataType) == typeid(float)) {
        float *dA=(float *)A, *dB=(float *)B, *dC=(float *)C;
        if(!a_dev) rmg_device_pool->malloc(&dA, a_size);
        if(!b_dev) rmg_device_pool->malloc(&dB, b_size);
        if(!c_dev) rmg_device_pool->malloc(&dC, c_size);
        if(!a_dev) rmg::error(hipMemcpyHtoD(dA, A, a_size * sizeof(float)));
        if(!b_dev) rmg::error(hipMemcpyHtoD(dB, B, b_size * sizeof(float)));
        if(!c_dev && std::abs(beta) != 0.0) rmg::error(hipMemcpyHtoD(dC, C, c_size * sizeof(float)));
        rmg::error(hipblasSgemm(ct.hipblas_handle, hip_transA, hip_transB, m, n, k,
                            (float*)&alpha,
                            (float*)dA, lda,
                            (float*)dB, ldb,
                            (float*)&beta, (float*)dC, ldc ));
        if(!c_dev) rmg::error(hipMemcpyDtoH(C, dC, c_size * sizeof(float)));
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
        if(!c_dev && std::abs(beta) != 0.0) rmg::error(hipMemcpyHtoD(dC, C, c_size * sizeof(double)));
        rmg::error(hipblasDgemm(ct.hipblas_handle, hip_transA, hip_transB, m, n, k,
                            (double*)&alpha,
                            (double*)dA, lda,
                            (double*)dB, ldb,
                            (double*)&beta, (double*)dC, ldc ));
        if(!c_dev) rmg::error(hipMemcpyDtoH(C, dC, c_size * sizeof(double)));
        if(!c_dev) rmg_device_pool->free(dC);
        if(!b_dev) rmg_device_pool->free(dB);
        if(!a_dev) rmg_device_pool->free(dA);

    }

    //hipDeviceSynchronize();
#elif SYCL_ENABLED

sycl_gemm(
    ct.sycl_Q,
    transa,
    transb,
    m, n, k,
    alpha,
    A, lda,
    B, ldb,
    beta,
    C, ldc
);
#else

    
//    RmgTimer *RT = new RmgTimer("gemmmmm ");
    if(typeid(DataType) == typeid(std::complex<double>))
    {
        if(ct.use_alt_zgemm)
            MyZgemm(transa, transb, m, n, k, (std::complex<double> *)(&alpha), (std::complex<double> *)A, lda, 
             (std::complex<double> *)B, ldb, (std::complex<double> *)(&beta), (std::complex<double> *)C, ldc);
        else
            zgemm(transa, transb, &m, &n, &k, (std::complex<double> *)&alpha, (std::complex<double> *)A, &lda,
            (std::complex<double> *)B, &ldb, (std::complex<double> *)&beta, (std::complex<double> *)C, &ldc);
    }
    else if(typeid(DataType) == typeid(std::complex<float>)) {
            cgemm(transa, transb, &m, &n, &k, (std::complex<float> *)&alpha, (std::complex<float> *)A, &lda,
             (std::complex<float> *)B, &ldb, (std::complex<float> *)(&beta), (std::complex<float> *)C, &ldc);
    }
    else if(typeid(DataType) == typeid(float)) {
        sgemm(transa, transb, &m, &n, &k, (float *)&alpha, (float *)A, &lda, 
        (float *)B, &ldb, (float *)&beta, (float *)C, &ldc);
    }
    else {
        dgemm(transa, transb, &m, &n, &k, (double *)&alpha, (double *)A, &lda, 
        (double *)B, &ldb, (double *)&beta, (double *)C, &ldc);
    }
//    delete RT;

#endif
}
