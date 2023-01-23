#pragma once
#include "gpu_blas_api.hpp"
#include "gpu_blas_handle.hpp"
#include "device_stream.hpp"
#include "device_event.hpp"
#include "gpu_runtime_api.hpp"
#include <vector>
#include <stdexcept>
#include "util.hpp"

namespace gpu {

class gpu_context {
public:
    gpu_context(int streams);

    gpu_context(gpu_context&&) = delete;
    gpu_context(gpu_context&) = delete;
    // context& operator=(context&& other) = delete;

    blas_api::HandleType get_blas_handle(int stream_id) const;

    runtime_api::StreamType get_stream(int stream_id) const;
    device_stream& get_device_stream(int stream_id);

    device_event enqueue_event(int stream_id) const;

    int get_num_streams() const;
    void set_num_streams(int streams);

    device_stream& get_result_stream();

    ~gpu_context();

private:
    int n_streams;
    std::vector<gpu_blas_handle> handles;
    std::vector<device_stream> streams;
    device_stream result_stream;
};

}
