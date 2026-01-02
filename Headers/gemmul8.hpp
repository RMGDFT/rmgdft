#pragma once
#define __NVCC__ 1
#if defined(__NVCC__)
    #include <cuComplex.h>
    #include <cublas_v2.h>
    #include <cuda_runtime.h>
#endif
#if defined(__HIPCC__)
    #include <hip/hip_complex.h>
    #include <hip/hip_runtime.h>
    #include <hipblas/hipblas.h>
#endif
#include <vector>

namespace gemmul8 {

/***
 * workSize returns the required workspace size in bytes.
 */
template <bool is_Complex = false> size_t workSize(
    const size_t m,                         // Number of rows of C
    const size_t n,                         // Number of columns of C
    const size_t k,                         // Inner dimension <= 2^17
    const unsigned num_moduli,              // #moduli, 2 <= num_moduli <= 20
    const bool enable_skip_scalA = false,   // [option] Reserve extra space for A to allow skip_scalA
    const bool enable_skip_scalB = false,   // [option] Reserve extra space for B to allow skip_scalB
    size_t *workSizeA            = nullptr, // [option] Output: workspace size used for A8i and sftA
    size_t *workSizeB            = nullptr  // [option] Output: workspace size used for B8i and sftB
);

/***
 * GEMM emulation using INT8 Tensor Cores
 */
#if defined(__NVCC__)
template <typename T> std::vector<double> gemm(
    cublasHandle_t handle,                  // Handle to the cuBLAS library context
    const cublasOperation_t op_A,           // CUBLAS_OP_N or CUBLAS_OP_T
    const cublasOperation_t op_B,           // CUBLAS_OP_N or CUBLAS_OP_T
    const size_t m,                         // Number of rows of C
    const size_t n,                         // Number of columns of C
    const size_t k,                         // Inner dimension <= 2^17
    const T *alpha,                         // Scaling factor for op(A)*op(B)
    const T *const A,                       // 1-D device array of dimensions lda*k (CUBLAS_OP_N) or lda*m (CUBLAS_OP_T)
    const size_t lda,                       // Leading dimension of A
    const T *const B,                       // 1-D device array of dimensions ldb*n (CUBLAS_OP_N) or ldb*k (CUBLAS_OP_T)
    const size_t ldb,                       // Leading dimension of B
    const T *beta,                          // Scaling factor for C
    T *const C,                             // 1-D device array of dimensions ldc*n
    const size_t ldc,                       // Leading dimension of C
    const unsigned num_moduli,              // #moduli, 2 <= num_moduli <= 20
    const bool fastmode,                    // false (accurate mode) or true (fast mode)
    void *const work,                       // Preallocated workspace
    void *const workA            = nullptr, // [optional] Separate workspace for A (if nullptr, uses work)
    void *const workB            = nullptr, // [optional] Separate workspace for B (if nullptr, uses work)
    const bool enable_skip_scalA = false,   // [optional] Enables scaling-skip mechanism for A
    const bool enable_skip_scalB = false,   // [optional] Enables scaling-skip mechanism for B
    const bool skip_scalA        = false,   // [optional] If true, skip preprocessing for A
    const bool skip_scalB        = false    // [optional] If true, skip preprocessing for B
);
#endif

#if defined(__HIPCC__)
template <typename T> std::vector<double> gemm(
    hipblasHandle_t handle,                 // Handle to the cuBLAS library context
    const hipblasOperation_t op_A,          // CUBLAS_OP_N or CUBLAS_OP_T
    const hipblasOperation_t op_B,          // CUBLAS_OP_N or CUBLAS_OP_T
    const size_t m,                         // Number of rows of C
    const size_t n,                         // Number of columns of C
    const size_t k,                         // Inner dimension <= 2^17
    const T *alpha,                         // Scaling factor for op(A)*op(B)
    const T *const A,                       // 1-D device array of dimensions lda*k (CUBLAS_OP_N) or lda*m (CUBLAS_OP_T)
    const size_t lda,                       // Leading dimension of A
    const T *const B,                       // 1-D device array of dimensions ldb*n (CUBLAS_OP_N) or ldb*k (CUBLAS_OP_T)
    const size_t ldb,                       // Leading dimension of B
    const T *beta,                          // Scaling factor for C
    T *const C,                             // 1-D device array of dimensions ldc*n
    const size_t ldc,                       // Leading dimension of C
    const unsigned num_moduli,              // #moduli, 2 <= num_moduli <= 20
    const bool fastmode,                    // false (accurate mode) or true (fast mode)
    void *const work,                       // Preallocated workspace
    void *const workA            = nullptr, // [optional] Separate workspace for A (if nullptr, uses work)
    void *const workB            = nullptr, // [optional] Separate workspace for B (if nullptr, uses work)
    const bool enable_skip_scalA = false,   // [optional] Enables scaling-skip mechanism for A
    const bool enable_skip_scalB = false,   // [optional] Enables scaling-skip mechanism for B
    const bool skip_scalA        = false,   // [optional] If true, skip preprocessing for A
    const bool skip_scalB        = false    // [optional] If true, skip preprocessing for B
);
#endif

} // namespace gemmul8
