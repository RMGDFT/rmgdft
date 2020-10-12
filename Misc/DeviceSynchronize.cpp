

#if HIP_ENABLED
#include <hip/hip_runtime.h>
void DeviceSynchronize(void)
{
    hipDeviceSynchronize();
}
#elif CUDA_ENABLED
#include <cuda.h>
void DeviceSynchronize(void)
{
    DeviceSynchronize();
}
#else
void DeviceSynchronize(void)
{
}
#endif
