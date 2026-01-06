
#include "blas_driver.h"
#include "GpuAlloc.h"


#if HIP_ENABLED
#include <hip/hip_runtime.h>
namespace rmg
{
    void sync_device(void)
    {
        hipDeviceSynchronize();
    }
} // end namespace rmg

#elif SYCL_ENABLED
#include <complex>
#include <typeinfo>
#include <string.h>
#include "const.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "rmg_control.h"

namespace rmg
{
    void sync_device(void)
    {
        ct.sycl_Q.wait();
    }
} // end namespace rmg

#elif CUDA_ENABLED
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cuda_runtime.h>
#include <driver_types.h>

namespace rmg
{
    void sync_device(void)
    {
        cudaDeviceSynchronize();
    }
}

#else
namespace rmg
{
    void sync_device(void)
    {
    }
}
#endif

