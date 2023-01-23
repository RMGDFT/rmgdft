#pragma once
#include <iostream>
#include <string>
#include <cmath>
#include "gpu_runtime_api.hpp"
#include "gpu_blas_api.hpp"

namespace gpu {

///////////////////////////////////////////////////////////////////////////////
// CUDA error checking
///////////////////////////////////////////////////////////////////////////////
static void check_runtime_status(runtime_api::StatusType status) {
    if(status !=  runtime_api::status::Success) {
        std::cerr << "error: GPU API call : "
        << runtime_api::get_error_string(status) << std::endl;
        throw(std::runtime_error("GPU ERROR"));
    }
}

static void check_blas_status(blas_api::StatusType status) {
    if(status != blas_api::status::Success) {
        auto error = blas_api::status::get_string(status); 
        std::cerr << "error: BLAS API call: " << error << std::endl;
        throw(std::runtime_error("GPU ERROR"));
    }
}

static void check_last_device_kernel(std::string const& errstr) {
    auto status = runtime_api::get_last_error();
    if(status != runtime_api::status::Success) {
        std::cout << "error: GPU kernel launch : " << errstr << " : "
        << runtime_api::get_error_string(status) << std::endl;
        throw(std::runtime_error("GPU ERROR"));
    }
}

///////////////////////////////////////////////////////////////////////////////
// TOTAL AVAILABLE MEMORY ON GPU 
///////////////////////////////////////////////////////////////////////////////
inline
std::size_t gpu_allocated_memory() {
    runtime_api::device_synchronize();
    auto status = runtime_api::get_last_error();
    check_runtime_status(status);
    std::size_t free;
    std::size_t total;
    status = runtime_api::mem_get_info(&free, &total);
    return status == runtime_api::status::Success ? total-free : -1;
}

///////////////////////////////////////////////////////////////////////////////
// allocating memory
///////////////////////////////////////////////////////////////////////////////

// allocate space on GPU for n instances of type T
template <typename T>
T* malloc_device(size_t n) {
    void* p;
    auto status = runtime_api::malloc(&p, n*sizeof(T));
    check_runtime_status(status);
    return (T*)p;
}

template <typename T>
T* malloc_pinned(size_t N, T value=T()) {
    T* ptr;
    auto status = runtime_api::host_alloc((void**)&ptr, N*sizeof(T), 0);
    check_runtime_status(status);
    std::fill(ptr, ptr+N, value);
    return ptr;
}

///////////////////////////////////////////////////////////////////////////////
// copying memory
///////////////////////////////////////////////////////////////////////////////

// copy n*T from host to device
template <typename T>
void copy_to_device(const T* from, T* to, size_t n) {
    runtime_api::memcpy(to, from, n*sizeof(T), runtime_api::flag::MemcpyHostToDevice);
}

// copy n*T from device to host
template <typename T>
void copy_to_host(const T* from, T* to, size_t n) {
    runtime_api::memcpy(to, from, n*sizeof(T), runtime_api::flag::MemcpyDeviceToHost);
}

// copy n*T from host to device
// If a cuda stream is passed as the final argument the copy will be performed
// asynchronously in the specified stream, otherwise it will be serialized in
// the default (NULL) stream
template <typename T>
void copy_to_device_async(const T* from, T* to, size_t n, runtime_api::StreamType stream=NULL) {
    //cudaDeviceSynchronize();
    // auto status = cudaGetLastError();
    // if(status != cudaSuccess) {
    //    std::cout << "error: CUDA kernel launch:"
    //    << cudaGetErrorString(status) << std::endl;
    //    throw(std::runtime_error("CUDA ERROR"));
    //}

    auto status = runtime_api::memcpy_async(to, from, n * sizeof(T),
                                            runtime_api::flag::MemcpyHostToDevice, stream);
    check_runtime_status(status);
}

// copy n*T from device to host
// If a cuda stream is passed as the final argument the copy will be performed
// asynchronously in the specified stream, otherwise it will be serialized in
// the default (NULL) stream
template <typename T>
void copy_to_host_async(const T* from, T* to, size_t n, runtime_api::StreamType stream=NULL) {
    auto status = runtime_api::memcpy_async(to, from, n * sizeof(T),
                                            runtime_api::flag::MemcpyDeviceToHost, stream);
    check_runtime_status(status);
}
}
