#pragma once
#include "gpu_context.hpp"
#include "device_buffer.hpp"
#include "device_stream.hpp"
#include "gpu_blas_handle.hpp"
#include <tuple>
#include <memory>

namespace gpu{

template <typename Scalar>
class mm_handle {
public:
    // the constructor just reserves the device vectors, without resizing
    mm_handle(int streams, int max_tile_m, int max_tile_n, int max_tile_k);
    ~mm_handle();

    mm_handle(mm_handle&&) = delete;
    mm_handle(const mm_handle&) = delete;
    mm_handle& operator=(const mm_handle&& other) = delete;

    void set_num_streams(int streams);
    int get_num_streams();

    gpu_context& get_gpu_context();

    // these resize the device vectors
    void set_tile_sizes(int tile_size_m, int tile_size_n, int tile_size_k);
    void set_tile_sizes(int tile_size);
    void set_full_sizes(int m, int n, int k);
    // returns the tile sizes that are actually used for 
    // given problem size: m, n, k
    std::tuple<int, int, int> optimal_tile_sizes(int m, int n, int k);
    // returns the allocated tile sizes
    std::tuple<int, int, int> get_max_tile_sizes();

    void set_streams_and_tiles(int streams, int tile_size_m, int tile_size_n, int tile_size_k);

    device_buffer<Scalar>& get_device_buffer_a();
    device_buffer<Scalar>& get_device_buffer_b();
    device_buffer<Scalar>& get_device_buffer_c();
    device_vector<Scalar>& get_full_device_buffer_c();

private:
    int n_streams = 2;
    int max_tile_size_m = 5000;
    int max_tile_size_n = 5000;
    int max_tile_size_k = 5000;

    gpu_context ctx;

    // device_vector<Scalar> memory_pool;

    device_buffer<Scalar> a_buff;
    device_buffer<Scalar> b_buff;
    device_buffer<Scalar> c_buff; // this holds just a single tile of C for each stream
    device_vector<Scalar> full_c_buff; // this holds the full C matrix
};

template <typename Scalar>
std::unique_ptr<mm_handle<Scalar>> make_context(int streams, 
                                                int max_tile_m, 
                                                int max_tile_n, 
                                                int max_tile_k) {
    return std::make_unique<mm_handle<Scalar>>(streams, 
                                               max_tile_m, 
                                               max_tile_n, 
                                               max_tile_k);
}
template <typename Scalar>
std::unique_ptr<mm_handle<Scalar>> make_context() {
    return std::make_unique<mm_handle<Scalar>>(2, 
                                               5000, 
                                               5000, 
                                               5000);
}
}
