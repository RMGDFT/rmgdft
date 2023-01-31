#include "tiled_mm.hpp"

#include "util.hpp"
#include "device_stream.hpp"
#include "device_event.hpp"
#include "gpu_blas_handle.hpp"
#include "tile_coord.hpp"
#include "gpu_context.hpp"
#include "device_buffer.hpp"
#include "tiled_matrix.hpp"
#include "tile_coord.hpp"
#include "gpu_blas_api.hpp"
#include "gpu_runtime_api.hpp"

// #include <omp.h>
// #include <cublasXt.h>
// #include <libsci_acc.h>

#include <vector>
#include <cstring>
#include <iostream>
#include <cmath>
#include <cstdio>
#include <complex>
#include <tuple>

namespace gpu {

using zfloat = std::complex<float>;
using zdouble = std::complex<double>;

int get_first(char trans, int m, int n) {
    return trans == 'N' ? m : n;
}

int get_second(char trans, int m, int n) {
    return trans == 'N' ? n : m;
}

gpu::tile_coord swap_if_transposed(char trans, int m, int n) {
    return trans == 'N' ? gpu::tile_coord(m, n) : gpu::tile_coord(n, m);
}


template<typename Scalar>
void copy_tile_to_device_async(tiled_matrix<Scalar>& tiled_mat, Scalar* d_buffer,
        tile_coord tile, device_stream& stream) {
    Scalar* from = tiled_mat.tile_data(tile);
    Scalar* to = d_buffer;

    tile_dim tile_dims = tiled_mat.tile_dimensions(tile);
    // std::cout << "host->device" << std::endl;

    auto status=
    runtime_api::memcpy_2d_async(to, tile_dims.rows() * sizeof(Scalar),
            from, tiled_mat.leading_dim() * sizeof(Scalar),
            tile_dims.rows() * sizeof(Scalar), tile_dims.cols(),
            runtime_api::flag::MemcpyHostToDevice, stream.stream());
    check_runtime_status(status);
}

template<typename Scalar>
void copy_tile_to_device_async(tiled_matrix<Scalar>& tiled_mat, device_buffer<Scalar>& d_buffer,
        tile_coord tile, gpu_context& ctx, int stream_id) {
    copy_tile_to_device_async(tiled_mat, d_buffer.stream_buffer(stream_id), tile, ctx.get_device_stream(stream_id));
}

template<typename Scalar>
void strided_copy_tile_to_device_async(
        tiled_matrix<Scalar>& tiled_mat, 
        Scalar* d_buffer,
        tile_coord tile, 
        int device_stride, 
        device_stream& stream) {
    Scalar* from = tiled_mat.tile_data(tile);
    Scalar* to = d_buffer;

    tile_dim tile_dims = tiled_mat.tile_dimensions(tile);
    // std::cout << "host->device" << std::endl;

    auto status=
    runtime_api::memcpy_2d_async(to, device_stride * sizeof(Scalar),
            from, tiled_mat.leading_dim() * sizeof(Scalar),
            tile_dims.rows() * sizeof(Scalar), tile_dims.cols(),
            runtime_api::flag::MemcpyHostToDevice, stream.stream());
    check_runtime_status(status);
}

template<typename Scalar>
void strided_copy_tile_to_device_async(
        tiled_matrix<Scalar>& tiled_mat, 
        Scalar* d_buffer,
        tile_coord tile, 
        int stride_device, 
        gpu_context& ctx, 
        int stream_id) {
    strided_copy_tile_to_device_async(tiled_mat, d_buffer, tile, stride_device, ctx.get_device_stream(stream_id));
}

template<typename Scalar>
void copy_tile_to_host_async(tiled_matrix<Scalar>& tiled_mat,
                             Scalar* d_buffer,
                             tile_coord tile,
                             device_stream& stream) {
    Scalar* from  = d_buffer;
    Scalar* to = tiled_mat.tile_data(tile);

    tile_dim tile_dims = tiled_mat.tile_dimensions(tile);

    // std::cout << "device->host" << std::endl;
    auto status=
    runtime_api::memcpy_2d_async(to, tiled_mat.leading_dim() * sizeof(Scalar),
            from, tile_dims.rows() * sizeof(Scalar),
            tile_dims.rows() * sizeof(Scalar), tile_dims.cols(),
            runtime_api::flag::MemcpyDeviceToHost, stream.stream());
    check_runtime_status(status);
}

template<typename Scalar>
void copy_tile_to_host_async(tiled_matrix<Scalar>& tiled_mat, device_buffer<Scalar>& d_buffer,
        tile_coord tile, gpu_context& ctx, int stream_id) {
    copy_tile_to_host_async(tiled_mat, d_buffer.stream_buffer(stream_id), tile, ctx.get_device_stream(stream_id));
}


template<typename Scalar>
std::tuple<int, int, int> get_num_tiles(tiled_matrix<Scalar>& a, tiled_matrix<Scalar>& b, tiled_matrix<Scalar>& c,
		                        int m, int n, int k,
					char trans_a, char trans_b) {
    int tile_size_m, tile_size_n, tile_size_k;
    auto tile_c = c.tile_dimensions();
    tile_size_m = tile_c.rows();
    tile_size_n = tile_c.cols();

    auto a_tile = a.tile_dimensions();
    int temp_m = a_tile.rows();
    int temp_k = a_tile.cols();

    tile_size_k = trans_a == 'N' ? temp_k : temp_m;

    int n_tiles_m = (int) std::ceil(1.0 * m / tile_size_m);
    int n_tiles_n = (int) std::ceil(1.0 * n / tile_size_n);
    int n_tiles_k = (int) std::ceil(1.0 * k / tile_size_k);

    return {n_tiles_m, n_tiles_n, n_tiles_k};
}


template<typename Scalar>
std::tuple<int, int, int> get_tile_sizes(tiled_matrix<Scalar>& a,
        tiled_matrix<Scalar>& b, tiled_matrix<Scalar>& c,
        int m_tile_id, int n_tile_id, int k_tile_id,
	char trans_a, char trans_b) {
    int tile_size_m, tile_size_n, tile_size_k;
    auto tile_c = c.tile_dimensions();
    tile_size_m = tile_c.rows();
    tile_size_n = tile_c.cols();

    auto a_tile = a.tile_dimensions();
    int temp_m = a_tile.rows();
    int temp_k = a_tile.cols();

    tile_size_k = trans_a == 'N' ? temp_k : temp_m;
    return {tile_size_m, tile_size_n, tile_size_k};
}


blas_api::OperationType get_blas_operation(char trans) {
    blas_api::OperationType op = trans == 'T'
	                            ? 
				    blas_api::operation::Transpose 
			            : 
				    (trans == 'C'
			                ? 
					blas_api::operation::ConjugateTranspose 
			                : 
					blas_api::operation::None);
    return op;
}

blas_api::StatusType cublas_gemm_wrapper(blas_api::HandleType handle,
		                   char trans_a, char trans_b,
                                   int m, int n, int k,
                                   const float* alpha,
                                   const float* a,
                                   const float* b,
                                   const float* beta,
                                   float* c, 
                                   int lld_c) {
    blas_api::OperationType op_a = get_blas_operation(trans_a);
    blas_api::OperationType op_b = get_blas_operation(trans_b);

    int ld_a = get_first(trans_a, m , k);
    int ld_b = get_first(trans_b, k , n);

    return blas_api::sgemm(handle, op_a, op_b, m, n, k,
                         alpha, a, ld_a, b, ld_b, beta, c, lld_c);
}

blas_api::StatusType cublas_gemm_wrapper(blas_api::HandleType handle,
		                   char trans_a, char trans_b,
                                   int m, int n, int k,
                                   const double* alpha,
                                   const double* a,
                                   const double* b,
                                   const double* beta,
                                   double* c,
                                   int lld_c) {
    blas_api::OperationType op_a = get_blas_operation(trans_a);
    blas_api::OperationType op_b = get_blas_operation(trans_b);

    int ld_a = get_first(trans_a, m , k);
    int ld_b = get_first(trans_b, k , n);

    return blas_api::dgemm(handle, op_a, op_b, m, n, k,
                         alpha, a, ld_a, b, ld_b, beta, c, lld_c);
}

// Note: Converting from std::complex to cuComplex and cuDoubleComple
//       works because they are binary compatible.
//
//       http://icl.cs.utk.edu/magma/forum/viewtopic.php?f=2&t=902
//
blas_api::StatusType cublas_gemm_wrapper(blas_api::HandleType handle,
		                   char trans_a, char trans_b,
                                   int m, int n, int k,
                                   const zfloat* alpha,
                                   const zfloat* a,
                                   const zfloat* b,
                                   const zfloat* beta,
                                   zfloat* c,
                                   int lld_c) {
    blas_api::OperationType op_a = get_blas_operation(trans_a);
    blas_api::OperationType op_b = get_blas_operation(trans_b);

    int ld_a = get_first(trans_a, m , k);
    int ld_b = get_first(trans_b, k , n);

    return blas_api::cgemm(handle, op_a, op_b, m, n, k,
                         reinterpret_cast<const blas_api::ComplexFloatType*>(alpha),
                         reinterpret_cast<const blas_api::ComplexFloatType*>(a), ld_a,
                         reinterpret_cast<const blas_api::ComplexFloatType*>(b), ld_b,
                         reinterpret_cast<const blas_api::ComplexFloatType*>(beta),
                         reinterpret_cast<blas_api::ComplexFloatType*>(c), lld_c);
}

blas_api::StatusType cublas_gemm_wrapper(blas_api::HandleType handle,
		                   char trans_a, char trans_b,
                                   int m, int n, int k,
                                   const zdouble* alpha,
                                   const zdouble* a,
                                   const zdouble* b,
                                   const zdouble* beta,
                                   zdouble* c,
                                   int lld_c) {
    blas_api::OperationType op_a = get_blas_operation(trans_a);
    blas_api::OperationType op_b = get_blas_operation(trans_b);

    int ld_a = get_first(trans_a, m , k);
    int ld_b = get_first(trans_b, k , n);

    return blas_api::zgemm(handle, op_a, op_b, m, n, k,
                         reinterpret_cast<const blas_api::ComplexDoubleType*>(alpha),
                         reinterpret_cast<const blas_api::ComplexDoubleType*>(a), ld_a,
                         reinterpret_cast<const blas_api::ComplexDoubleType*>(b), ld_b,
                         reinterpret_cast<const blas_api::ComplexDoubleType*>(beta),
                         reinterpret_cast<blas_api::ComplexDoubleType*>(c), lld_c);
}

template<typename Scalar>
void round_robin(tiled_matrix<Scalar>& a_host, tiled_matrix<Scalar>& b_host, tiled_matrix<Scalar>& c_host,
        device_buffer<Scalar>& a_device,
        device_buffer<Scalar>& b_device,
        device_buffer<Scalar>& c_device,
        char trans_a, char trans_b,
	int m, int n, int k, 
	Scalar alpha, Scalar beta, mm_handle<Scalar>& handle) {

    int n_tiles_m, n_tiles_n, n_tiles_k;
    std::tie(n_tiles_m, n_tiles_n, n_tiles_k) = get_num_tiles(a_host, b_host, c_host, 
		                                              m, n, k, 
							      trans_a, trans_b);

    int n_streams = std::min(handle.get_num_streams(), n_tiles_m * n_tiles_n);
    auto& gpu_ctx = handle.get_gpu_context();

    auto& result_stream = gpu_ctx.get_result_stream();

    std::vector<device_event> c_computed_on_device(n_streams);
    std::vector<device_event> c_copied_to_host(n_streams);

    for (int i = 0; i < n_tiles_m * n_tiles_n; i += n_streams) {
        for (int k_tile_id = 0; k_tile_id < n_tiles_k; ++k_tile_id) {
            for (int round = 0; round < 2; ++round) {
                int current_i = i;

                for (int stream_id = 0; stream_id < n_streams
                        && current_i < n_tiles_m * n_tiles_n; ++stream_id) {

                    int m_tile_id = current_i / n_tiles_n;
                    int n_tile_id = current_i % n_tiles_n;

                    int actual_size_m, actual_size_n, actual_size_k;
                    std::tie(actual_size_m, actual_size_n, actual_size_k) =
                        get_tile_sizes(a_host, b_host, c_host,
                                m_tile_id, n_tile_id, k_tile_id,
				trans_a, trans_b);

                    Scalar new_beta = k_tile_id == 0 ? beta : Scalar{1};

                    auto& current_stream = gpu_ctx.get_device_stream(stream_id);

                    if (round == 0) {
                        // copy A tile
                        copy_tile_to_device_async(a_host, a_device,
                                swap_if_transposed(trans_a, m_tile_id, k_tile_id),
                                gpu_ctx, stream_id);

                        // copy B tile
                        copy_tile_to_device_async(b_host, b_device,
                                swap_if_transposed(trans_b, k_tile_id, n_tile_id),
                                gpu_ctx, stream_id);

                        // copy C tile if this is the first partial result and beta > 0
                        if (k_tile_id == 0 && std::abs(beta) > 0) {
                            current_stream.wait_on_event(c_copied_to_host[stream_id]);
                            copy_tile_to_device_async(c_host, c_device,
				    {m_tile_id, n_tile_id},
                                    gpu_ctx, stream_id);
                        }
                    } else {
                        // perform dgemm
                        // cublasSetStream(get_blas_handle(stream_id), streams[stream_id].stream());
                        // std::cout << "performing dgemm" << std::endl;
                        if (k_tile_id == 0) {
                            current_stream.wait_on_event(c_copied_to_host[stream_id]);
                        }

                        auto status = cublas_gemm_wrapper(
                                gpu_ctx.get_blas_handle(stream_id),
				trans_a, trans_b,
                                actual_size_m, actual_size_n, actual_size_k,
                                &alpha,
                                a_device.stream_buffer(stream_id),
                                b_device.stream_buffer(stream_id),
                                &new_beta,
                                c_device.stream_buffer(stream_id), actual_size_m);
                        check_blas_status(status);

                        if (k_tile_id == n_tiles_k - 1) {
                            c_computed_on_device[stream_id] = current_stream.enqueue_event();
                            // copy result back to host
                            result_stream.wait_on_event(c_computed_on_device[stream_id]);
                            copy_tile_to_host_async(c_host, c_device.stream_buffer(stream_id),
			            {m_tile_id, n_tile_id},
                                    result_stream);
                            c_copied_to_host[stream_id] = result_stream.enqueue_event();
                        }
                    }
                    current_i++;
                }
            }
        }
    }
}

template<typename Scalar>
void round_robin_without_copy_c(tiled_matrix<Scalar>& a_host, tiled_matrix<Scalar>& b_host, tiled_matrix<Scalar>& c_host,
        device_buffer<Scalar>& a_device,
        device_buffer<Scalar>& b_device,
        tiled_matrix<Scalar>& tiled_c_device,
	char trans_a, char trans_b,
        int m, int n, int k, Scalar alpha, Scalar beta, mm_handle<Scalar>& handle) {

    int n_tiles_m, n_tiles_n, n_tiles_k;
    std::tie(n_tiles_m, n_tiles_n, n_tiles_k) = get_num_tiles(a_host, b_host, c_host, 
		                                              m, n, k, 
							      trans_a, trans_b);

    int n_streams = std::min(handle.get_num_streams(), n_tiles_m * n_tiles_n);
    auto& gpu_ctx = handle.get_gpu_context();

    std::vector<device_event> copied_to_device(n_streams);

    for (int i = 0; i < n_tiles_m * n_tiles_n; i += n_streams) {
        for (int k_tile_id = 0; k_tile_id < n_tiles_k; ++k_tile_id) {
            for (int round = 0; round < 2; ++round) {
                int current_i = i;

                for (int stream_id = 0; stream_id < n_streams
                        && current_i < n_tiles_m * n_tiles_n; ++stream_id) {

                    int m_tile_id = current_i / n_tiles_n;
                    int n_tile_id = current_i % n_tiles_n;

                    int actual_size_m, actual_size_n, actual_size_k;
                    std::tie(actual_size_m, actual_size_n, actual_size_k) =
                        get_tile_sizes(a_host, b_host, c_host,
                                m_tile_id, n_tile_id, k_tile_id,
				trans_a, trans_b);

                    Scalar new_beta = k_tile_id == 0 ? beta : Scalar{1};

                    auto& current_stream = gpu_ctx.get_device_stream(stream_id);

                    auto c_device_ptr = tiled_c_device.tile_data({m_tile_id, n_tile_id});

                    if (round == 0) {
                        int prev_stream_id = (stream_id + n_streams - 1) % n_streams;
                        current_stream.wait_on_event(copied_to_device[prev_stream_id]);

                        // copy A tile
                        copy_tile_to_device_async(a_host, a_device,
			        swap_if_transposed(trans_a, m_tile_id, k_tile_id),
                                gpu_ctx, stream_id);

                        // copy B tile
                        copy_tile_to_device_async(b_host, b_device,
			        swap_if_transposed(trans_b, k_tile_id, n_tile_id),
                                gpu_ctx, stream_id);

                        // copy C tile if this is the first partial result and beta > 0
                        if (k_tile_id == 0 && std::abs(beta) > 0) {
                            // current_stream.wait_on_event(c_copied_to_host[stream_id]);
                            strided_copy_tile_to_device_async(
                                    c_host, 
                                    c_device_ptr,
                                    {m_tile_id, n_tile_id},
                                    m,
                                    gpu_ctx, 
                                    stream_id);
                        }

                        copied_to_device[stream_id] = current_stream.enqueue_event();

                    } else {
                        // perform dgemm
                        auto status = cublas_gemm_wrapper(
                                gpu_ctx.get_blas_handle(stream_id),
				trans_a, trans_b,
                                actual_size_m, actual_size_n, actual_size_k,
                                &alpha,
                                a_device.stream_buffer(stream_id),
                                b_device.stream_buffer(stream_id),
                                &new_beta,
                                c_device_ptr, m);
                        check_blas_status(status);
                    }
                    current_i++;
                }
            }
        }
    }
}

/*
void gpu_dgemm_(mm_handle& m_handle, double* a, double* b, double* c,
        int m, int n, int k,
        double alpha, double beta) {

    int tile_size_m, tile_size_n, tile_size_k;
    std::tie(tile_size_m, tile_size_n, tile_size_k) = m_handle.get_tile_sizes();

    // perform dgemm
    dgemm('N', 'N', m, n, k, alpha, a, m, b, k, beta, c, m);
}

void gpu_dgemm_(mm_handle& m_handle, double* a, double* b, double* c,
        int m, int n, int k,
        double alpha, double beta) {
    cublasXtHandle_t handle;
    cublasXtCreate(&handle);
    int devices[1] = {0};
    cublasXtDeviceSelect(handle, 1, devices);

    int tile_size_m, tile_size_n, tile_size_k;
    std::tie(tile_size_m, tile_size_n, tile_size_k) = m_handle.get_tile_sizes();

    cublasXtSetBlockDim(handle, tile_size_m);

    // cublasXtSetCpuRoutine(handle, CUBLASXT_GEMM, CUBLASXT_DOUBLE, (void*)(&dgemm_));
    // perform dgemm
    cublasXtDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N,
        m, n, k, &alpha, a, m, b, k, &beta, c, m);

    cudaDeviceSynchronize();

    if (handle)
        cublasXtDestroy(handle);
}
*/
template<typename Scalar>
void gemm(mm_handle<Scalar>& handle, 
	  char transa, char transb,
          int m, int n, int k,
	  Scalar alpha,
	  Scalar* a, int ld_a,
	  Scalar* b, int ld_b, 
	  Scalar beta,
	  Scalar* c, int ld_c,
          bool pin_host_buffers, bool copy_c_back) {

    char trans_a = std::toupper(transa);
    char trans_b = std::toupper(transb);

    // sumatrix size to multiply
    int a_subm = trans_a == 'N' ? m : k;
    int a_subn = trans_a == 'N' ? k : m;

    int b_subm = trans_b == 'N' ? k : n;
    int b_subn = trans_b == 'N' ? n : k;

    int c_subm = m;
    int c_subn = n;

    // check the leading dimensions
    if (ld_a < a_subm) {
        std::runtime_error("[ERROR] Leading dimensions for matrix A must be >= the number of rows.");
    }
    // check the leading dimensions
    if (ld_b < b_subm) {
        std::runtime_error("[ERROR] Leading dimensions for matrix B must be >= the number of rows.");
    }
    // check the leading dimensions
    if (ld_c < c_subm) {
        std::runtime_error("[ERROR] Leading dimensions for matrix C must be >= the number of rows.");
    }

    if (pin_host_buffers) {
        // pin matrix A
        // auto start = std::chrono::steady_clock::now();
        auto status = gpu::runtime_api::host_register(
                a,
                ld_a * a_subn * sizeof(Scalar),
                gpu::runtime_api::flag::HostRegisterDefault);
        gpu::check_runtime_status(status);
        // pin matrix B
        status = gpu::runtime_api::host_register(
                b,
                ld_b * b_subn * sizeof(Scalar),
                gpu::runtime_api::flag::HostRegisterDefault);
        gpu::check_runtime_status(status);
        // pin matrix C
        if (copy_c_back == true || std::abs(beta) >0) {
            status = gpu::runtime_api::host_register(
                    c,
                    ld_c * c_subn * sizeof(Scalar),
                    gpu::runtime_api::flag::HostRegisterDefault);
            gpu::check_runtime_status(status);
        }
        // auto end = std::chrono::steady_clock::now();
        // auto time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        // std::cout << "Pinning time = " << time << std::endl;
    }

    int tile_size_m, tile_size_n, tile_size_k;
    std::tie(tile_size_m, tile_size_n, tile_size_k) = handle.optimal_tile_sizes(m, n, k);

    // set these tile sizes to GPU matrices
    handle.set_tile_sizes(tile_size_m, tile_size_n, tile_size_k);
    if (!copy_c_back) {
        handle.set_full_sizes(m, n, k);
    }

    // tile sizes
    int a_tile_subm = trans_a == 'N' ? tile_size_m : tile_size_k;
    int a_tile_subn = trans_a == 'N' ? tile_size_k : tile_size_m;

    int b_tile_subm = trans_b == 'N' ? tile_size_k : tile_size_n;
    int b_tile_subn = trans_b == 'N' ? tile_size_n : tile_size_k;

    int c_tile_subm = tile_size_m;
    int c_tile_subn = tile_size_n;

    // set these tile sizes to CPU matrices
    tiled_matrix<Scalar> a_host(a, a_subm, a_subn, ld_a, {a_tile_subm, a_tile_subn});
    tiled_matrix<Scalar> b_host(b, b_subm, b_subn, ld_b, {b_tile_subm, b_tile_subn});
    tiled_matrix<Scalar> c_host(c, c_subm, c_subn, ld_c, {c_tile_subm, c_tile_subn});

    // get GPU matrices
    device_buffer<Scalar>& a_device = handle.get_device_buffer_a();
    device_buffer<Scalar>& b_device = handle.get_device_buffer_b();
    device_buffer<Scalar>& c_device = handle.get_device_buffer_c();

    // used only when copy_c_back == false
    device_vector<Scalar>& full_c_device = handle.get_full_device_buffer_c();
    tiled_matrix<Scalar> tiled_c_device(full_c_device.data(), c_subm, c_subn, ld_c, {c_tile_subm, c_tile_subn});

    if (copy_c_back) {
        round_robin(a_host, b_host, c_host,
                    a_device, b_device, c_device,
		    trans_a, trans_b,
                    m, n, k, alpha, beta, handle);
    } else {
        round_robin_without_copy_c(
                    a_host, b_host, c_host,
                    a_device, b_device, tiled_c_device,
		    trans_a, trans_b,
                    m, n, k, alpha, beta, handle);
    }

    auto status =
    runtime_api::device_synchronize();
    check_runtime_status(status);

    if (pin_host_buffers) {
        // start = std::chrono::steady_clock::now();
        // unpin matrix A
        status = gpu::runtime_api::host_unregister(a);
        gpu::check_runtime_status(status);
        // unpin matrix B
        status = gpu::runtime_api::host_unregister(b);
        gpu::check_runtime_status(status);
        // unpin matrix C
        if (copy_c_back == true || std::abs(beta) >0) {
            status = gpu::runtime_api::host_unregister(c);
            gpu::check_runtime_status(status);
        }

        // end = std::chrono::steady_clock::now();
        // time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        // std::cout << "Unpinning time = " << time << std::endl;
    }
}

template void gemm<float>(
	mm_handle<float>& handle,
	char transa, char transb,
	int m, int n, int k,
	float alpha,
	float* a, int ld_a,
	float* b, int ld_b, 
	float beta,
	float* c, int ld_c,
	bool pin_host_buffers, bool copy_c_back);

template void gemm<double>(
	mm_handle<double>& handle,
	char transa, char transb,
        int m, int n, int k,
	double alpha,
	double* a, int ld_a,
	double* b, int ld_b, 
	double beta,
	double* c, int ld_c,
        bool pin_host_buffers, bool copy_c_back);

template void gemm<zfloat>(
	mm_handle<zfloat>& handle,
	char transa, char transb,
        int m, int n, int k,
	zfloat alpha,
	zfloat* a, int ld_a,
	zfloat* b, int ld_b, 
	zfloat beta,
	zfloat* c, int ld_c,
        bool pin_host_buffers, bool copy_c_back);

template void gemm<zdouble>(
	mm_handle<zdouble>& handle,
	char transa, char transb,
        int m, int n, int k,
	zdouble alpha,
	zdouble* a, int ld_a,
	zdouble* b, int ld_b, 
	zdouble beta,
	zdouble* c, int ld_c,
        bool pin_host_buffers, bool copy_c_back);

}
