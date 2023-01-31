#pragma once
#include <cmath>
#include "util.hpp"
#include "gpu_runtime_api.hpp"

namespace gpu {

// ****************************************************************** 
// class definition
// ****************************************************************** 
template <typename T>
class device_vector {
public:
    device_vector() = default;

    device_vector(std::size_t n);

    // copy-constructor are not supported
    device_vector(device_vector& other) = delete;
    // copy-operator disabled
    device_vector& operator=(device_vector&) = delete;

    // assignment operators are supported
    device_vector& operator=(device_vector&& other);

    T* data();
    std::size_t size();
    std::size_t capacity();

    void resize(std::size_t size);

    ~device_vector();

private:
    T* data_;
    std::size_t size_ = (std::size_t) 0;
    std::size_t capacity_ = (std::size_t) 0;
};

// ****************************************************************** 
// class implementation
// ****************************************************************** 
template <typename T>
device_vector<T>::device_vector(std::size_t n): 
    size_(n),
    capacity_((std::size_t) std::ceil(1.2 * n)) {
    data_ = malloc_device<T>(capacity_);
}

// assignment operators are supported
template <typename T>
device_vector<T>& device_vector<T>::operator=(device_vector<T>&& other) {
    if (this != &other) {
        if (this->capacity() > 0) {
            auto status = runtime_api::free(this->data_);
            check_runtime_status(status);
        }
        this->data_ = other.data_;
        this->size_ = other.size_;
        this->capacity_ = other.capacity_;
        other.size_ = 0;
        other.capacity_ = 0;
    }
    return *this;
}

template <typename T>
T* device_vector<T>::data() {
    return data_;
}

template <typename T>
std::size_t device_vector<T>::size() {
    return size_;
}

template <typename T>
std::size_t device_vector<T>::capacity() {
    return capacity_;
}

template <typename T>
device_vector<T>::~device_vector() {
    if (capacity() > 0) {
        auto status = runtime_api::free(data_);
        check_runtime_status(status);
        capacity_ = 0;
        size_ = 0;
    }
}

template <typename T>
void device_vector<T>::resize(std::size_t size) {
    if (size > 0) {
        if (size > capacity()) {
            if (capacity() > 0) {
                auto status = runtime_api::free(data_);
                check_runtime_status(status);
            }
            size_ = size;
            capacity_ = (std::size_t) std::ceil(1.2 * size);
            data_ = malloc_device<T>(capacity_);
        } else {
            size_ = size;
        }
    }
}
}
