#include <cstddef>
#include <utility>
#include <string>
#include "rmg_error.h"
#include "GpuAlloc.h"

template <typename T, int target>
class rmg_vector
{
public:
    using value_type = T;

    rmg_vector() = default;

    explicit rmg_vector(std::size_t n)
        : size_(n)
    {
        allocate();
    }

    // non-copyable (like unique_ptr)
    rmg_vector(const rmg_vector&) = delete;
    rmg_vector& operator=(const rmg_vector&) = delete;

    ~rmg_vector()
    {
        release();
    }

    T* data() noexcept { return data_; }
    const T* data() const noexcept { return data_; }

    std::size_t size() const noexcept { return size_; }

    // Host-only access (safe only for MallocHost!)
    T& operator[](std::size_t i) noexcept { return data_[i]; }
    const T& operator[](std::size_t i) const noexcept { return data_[i]; }

private:
    void allocate()
    {
        if(target == 0)
        {
            gpuMalloc((void **)&data_, sizeof(T) * size_);
        }
        else if (target == 1)
        {
            data_ = (T *)RmgMallocHost(sizeof(T) * size_);
        }
        else
        {
            rmg::error("rmg_vector accept device or host only");
        }

    }

    void release()
    {

        if(target == 0)
        {
            gpuFree(data_);
        }
        else if (target == 1)
        {
            RmgFreeHost(data_);
        }
        else
        {
            rmg::error("rmg_vector accept device or host only");
        }

        data_ = nullptr;
        size_ = 0;
    }

private:
    T* data_ = nullptr;
    std::size_t size_ = 0;
};

