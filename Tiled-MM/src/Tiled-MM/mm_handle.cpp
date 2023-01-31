#include "mm_handle.hpp"
#include "gpu_runtime_api.hpp"

#include<complex>
#include <cmath>
#include <cassert>

namespace gpu {

template <typename Scalar>
mm_handle<Scalar>::mm_handle(int streams, int max_tile_m, int max_tile_n, int max_tile_k): 
        n_streams(streams),
        max_tile_size_m(max_tile_m),
        max_tile_size_n(max_tile_n),
        max_tile_size_k(max_tile_k),
        ctx(streams) {

    // here we don't set the device because:
    //     - if the user has already set it, 
    //       we don't want to oevrwrite it
    //     - if the user didn't set it,
    //       it's anyway 0 by default
    // runtime_api::set_device(0);
    // status = runtime_api::get_device(&old_device);

    a_buff = device_buffer<Scalar>(n_streams);
    b_buff = device_buffer<Scalar>(n_streams);
    c_buff = device_buffer<Scalar>(n_streams);
}

template <typename Scalar>
mm_handle<Scalar>::~mm_handle() {
    // std::cout << "freeing mm_handle" << std::endl;
}

template <typename Scalar>
void mm_handle<Scalar>::set_num_streams(int streams) {
    n_streams = streams;
    ctx.set_num_streams(streams);

    // this will resize the vectors
    a_buff.set_num_streams(streams);
    b_buff.set_num_streams(streams);
    c_buff.set_num_streams(streams);
}

template <typename Scalar>
int mm_handle<Scalar>::get_num_streams() {
    return n_streams;
}

template <typename Scalar>
gpu_context& mm_handle<Scalar>::get_gpu_context() {
    return ctx;
}

template <typename Scalar>
void mm_handle<Scalar>::set_tile_sizes(int tile_m, int tile_n, int tile_k) {
    assert(tile_m <= max_tile_size_m);
    assert(tile_n <= max_tile_size_n);
    assert(tile_k <= max_tile_size_k);

    a_buff.set_tile_sizes({tile_m, tile_k});
    b_buff.set_tile_sizes({tile_k, tile_n});
    c_buff.set_tile_sizes({tile_m, tile_n});
}

template <typename Scalar>
void mm_handle<Scalar>::set_tile_sizes(int tile_size) {
    set_tile_sizes(tile_size, tile_size, tile_size);
}

template <typename Scalar>
void mm_handle<Scalar>::set_full_sizes(int m, int n, int k) {
    assert(m > 0);
    assert(n > 0);
    assert(k > 0);

    full_c_buff.resize(m * n);
}

// tries to: either grow the current tile size if necessary and possible
// or finds a possibly smaller tile size that 
// perfectly divides the dimension.
// If the found tile size is too small, then 
// stick to the max_tile_size.
// The divisibility is not a requirement of the algorithm
// but can improve the performance.
int optimal_tile_size(int dim, int curr_tile_size, int max_tile_size) {
    // if we can still grow the tile size, then grow if necessary
    if (curr_tile_size < max_tile_size && dim <= max_tile_size) {
        return dim;
    }

    // otherwise, we ignore the curr_tile_size and consider only max_tile_size
    int tile = 1;
    int limit = std::min(max_tile_size, dim);
    // iterate over divisors of dim
    for (int i = 1; i <= limit; ++i) {
        if (dim % i == 0) {
            tile = i;
        }
    }

    // if the divisible tile size is not smaller than half of the 
    // max_tile_size then use it, otherwise, stick to the max_tile_size
    if (std::abs(max_tile_size - tile) <= max_tile_size/2) 
        return tile;
    return limit;
}

template <typename Scalar>
std::tuple<int, int, int> mm_handle<Scalar>::get_max_tile_sizes() {
    return {max_tile_size_m, max_tile_size_n, max_tile_size_k};
}

template <typename Scalar>
std::tuple<int, int, int> mm_handle<Scalar>::optimal_tile_sizes(
        int m, int n, int k) {
    int tile_m = optimal_tile_size(
            m, // matrix dimension
            c_buff.get_tile_sizes().rows(), // current tile size
            max_tile_size_m); // maximum allowed tile size
    int tile_n = optimal_tile_size(
            n,
            c_buff.get_tile_sizes().cols(),
            max_tile_size_n);
    int tile_k = optimal_tile_size(
            k, 
            a_buff.get_tile_sizes().cols(), // current tile size
            max_tile_size_k);
    return {tile_m, tile_n, tile_k};
}

template <typename Scalar>
void mm_handle<Scalar>::set_streams_and_tiles(int streams, int tile_m, int tile_n, int tile_k) {
    assert(tile_m <= max_tile_size_m);
    assert(tile_n <= max_tile_size_n);
    assert(tile_k <= max_tile_size_k);
    n_streams = streams;
    ctx.set_num_streams(n_streams);
    a_buff.set_streams_and_tiles(streams, {tile_m, tile_k});
    b_buff.set_streams_and_tiles(streams, {tile_k, tile_n});
    c_buff.set_streams_and_tiles(streams, {tile_m, tile_n});
}

template <typename Scalar>
device_buffer<Scalar>& mm_handle<Scalar>::get_device_buffer_a() {
    return a_buff;
}

template <typename Scalar>
device_buffer<Scalar>& mm_handle<Scalar>::get_device_buffer_b() {
    return b_buff;
}

template <typename Scalar>
device_buffer<Scalar>& mm_handle<Scalar>::get_device_buffer_c() {
    return c_buff;
}

template <typename Scalar>
device_vector<Scalar>& mm_handle<Scalar>::get_full_device_buffer_c() {
    return full_c_buff;
}

template class mm_handle<float>;
template class mm_handle<double>;
template class mm_handle<std::complex<float>>;
template class mm_handle<std::complex<double>>;

}
