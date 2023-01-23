#pragma once

#include "util.hpp"
#include "gpu_blas_api.hpp"
#include "gpu_runtime_api.hpp"

namespace gpu {
// wrapper around cublasHandle
class gpu_blas_handle {
public:
    gpu_blas_handle() {
        // runtime_api::set_device(0);
        auto status = 
        blas_api::create(&handle_);
        check_blas_status(status);
        valid_ = true;
    }

    ~gpu_blas_handle() {
        if (valid_) {
            // std::cout << "freeing cublas handle" << std::endl;
            auto status =
            blas_api::destroy(handle_);
            check_blas_status(status);
        }
    }

    // move constructor
    gpu_blas_handle(gpu_blas_handle&& other) {
        handle_ = other.handle_;
        valid_ = other.valid_;
        other.valid_ = false;
    }

    // move-assignment operator
    gpu_blas_handle& operator=(gpu_blas_handle&& other) {
        if (this != &other) {
            if (valid_) {
                auto status = 
                blas_api::destroy(handle_);
                check_blas_status(status);
            }
            handle_ = other.handle_;
            valid_ = other.valid_;
            other.valid_ = false;
        }
        return *this;
    }

    // copy-constructor disabled
    gpu_blas_handle(gpu_blas_handle&) = delete;
    // copy-operator disabled
    gpu_blas_handle& operator=(gpu_blas_handle&) = delete;

    // return the unerlying cublas handle
    blas_api::HandleType handle() const {
        return handle_;
    }

private:
    bool valid_ = false;
    blas_api::HandleType handle_;
};
}
