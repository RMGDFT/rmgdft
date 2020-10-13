

#if HIP_ENABLED
#include <hip/hip_runtime.h>
void DeviceSynchronize(void)
{
    hipDeviceSynchronize();
}
hipError_t gpuStreamSynchronize (hipStream_t stream)
{
    return hipStreamSynchronize (stream);
}
#elif CUDA_ENABLED
#include <cuda.h>
void DeviceSynchronize(void)
{
    DeviceSynchronize();
}
cudaError_t gpuStreamSynchronize (cudaStream_t stream)
{
    return cudaStreamSynchronize (stream);
}
#else
void DeviceSynchronize(void)
{
}
#endif
