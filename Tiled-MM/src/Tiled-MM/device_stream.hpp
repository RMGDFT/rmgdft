#pragma once

#include "util.hpp"
#include "gpu_runtime_api.hpp"
#include "device_event.hpp"

namespace gpu {
// wrapper around device streams
class device_stream {
  public:
    device_stream():
        stream_(new_stream()),
        valid_(true)
    {};

    ~device_stream() {
        if (valid_) {
            // std::cout << "freeing device_stream" << std::endl;
            auto status = runtime_api::stream_destroy(stream_);
            check_runtime_status(status);
            // std::cout << "freed device stream" << std::endl;
        }
    }

    // return the device stream
    runtime_api::StreamType stream() const {
        return stream_;
    }

    // move-constructor
    device_stream(device_stream&& other) {
        stream_ = other.stream_;
        valid_ = other.valid_;
        other.valid_ = false;
    }

    // move-assignment operator
    device_stream& operator=(device_stream&& other) {
        if (this != &other) {
            if (valid_) {
                auto status = runtime_api::stream_destroy(stream_);
                check_runtime_status(status);
            }
            stream_ = other.stream_;
            valid_ = other.valid_;
            other.valid_ = false;
        }
        return *this;
    }

    // copy-constructor disabled
    device_stream(device_stream&) = delete;
    // copy-operator disabled
    device_stream& operator=(device_stream&) = delete;

    // insert event into stream
    // returns immediately
    device_event enqueue_event() const {
        device_event e;

        auto status = runtime_api::event_record(e.get(), stream_);
        check_runtime_status(status);

        return e;
    }

    // make all future work on stream wait until event has completed.
    // returns immediately, not waiting for event to complete
    void wait_on_event(device_event &e) const {
        auto status = runtime_api::stream_wait_event(stream_, e.get(), 0);
        check_runtime_status(status);
    }

  private:
    runtime_api::StreamType new_stream() {
        runtime_api::StreamType s;

        auto status = runtime_api::stream_create_with_flags(&s, runtime_api::flag::StreamNonBlocking);
        check_runtime_status(status);

        return s;
    }

    runtime_api::StreamType stream_;
    bool valid_;
};
}
