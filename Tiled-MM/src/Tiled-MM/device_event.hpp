#pragma once

#include "util.hpp"
#include "gpu_runtime_api.hpp"

namespace gpu {

// wrapper around device events
class device_event {
  public:

    device_event() {
        auto status = runtime_api::event_create_with_flags(&event_, runtime_api::flag::EventDisableTiming);
        check_runtime_status(status);
        valid_ = true;
    }

    ~device_event() {
        // note that EventDestroy can be called on an device_event before is has
        // been reached in a stream, and the GPU runtime will defer clean up
        // of the device_event until it has been completed.
        if (valid_) {
            auto status = runtime_api::event_destroy(event_);
            check_runtime_status(status);
        }
    }

    // move constructor
    device_event(device_event&& other) {
        event_ = other.event_;
        valid_ = other.valid_;
        other.valid_ = false;
    }

    // move-assignment operator
    device_event& operator=(device_event&& other) {
        if (this != &other) {
            if (valid_) {
                auto status = runtime_api::event_destroy(event_);
                check_runtime_status(status);
            }
            event_ = other.event_;
            valid_ = other.valid_;
            other.valid_ = false;
        }
        return *this;
    }

    // copy-constructor disabled
    device_event(device_event&) = delete;
    // copy-operator disabled
    device_event& operator=(device_event&) = delete;

    // return the underlying device_event handle
    runtime_api::EventType& get() {
        return event_;
    }

    // force host execution to wait for device_event completion
    void wait() {
        auto status = runtime_api::event_synchronize(event_);
        check_runtime_status(status);
    }

    // returns time in seconds taken between this device device_event and another device device_event
    double time_since(device_event& other) {
        float time_taken = 0.0f;

        auto status = runtime_api::event_elapsed_time(&time_taken, other.get(), event_);
        check_runtime_status(status);
        return double(time_taken/1.e3);
    }

  private:
    bool valid_ = false;
    runtime_api::EventType event_;
};
}

