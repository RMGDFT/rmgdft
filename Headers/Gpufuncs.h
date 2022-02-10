
#ifndef GPU_FUNCS_H
#define GPU_FUNCS_H 1

#include <complex>
template <typename T>
struct fdparms_o8 {
    T a0;
    T gpt1x;
    T gpt2x;
    T gpt3x;
    T gpt4x;

    T gmt1x;
    T gmt2x;
    T gmt3x;
    T gmt4x;

    T gpt1y;
    T gpt2y;
    T gpt3y;
    T gpt4y;

    T gmt1y;
    T gmt2y;
    T gmt3y;
    T gmt4y;

    T gpt1z;
    T gpt2z;
    T gpt3z;
    T gpt4z;

    T gmt1z;
    T gmt2z;
    T gmt3z;
    T gmt4z;

    T gpt1xh;
    T gpt2xh;
    T gpt3xh;
    T gpt4xh;

    T gmt1xh;
    T gmt2xh;
    T gmt3xh;
    T gmt4xh;
};

int getThreadId(void);

#if CUDA_ENABLED
#include <cublas_v2.h>

void init_cuda_fd(int max_threads, size_t bufsize);
void GpuFill(double *dptr, int n, double fillval);
void GpuNegate(double *dx, int incx, double *dy, int incy, int n);
void gramsch_update_psi(double *V,
                        double *C,
                        int N,
                        int eig_start,
                        int eig_stop,
                        cublasHandle_t cublasH);
void GpuEleMul(double *dx, std::complex<double> *dy, int n, cudaStream_t stream);
void GpuEleMul(double *dx, std::complex<float> *dy, int n, cudaStream_t stream);
template <typename T>
void app8_del2_gpu(T * __restrict__ a,
                   T *b,
                   const int dimx,
                   const int dimy,
                   const int dimz,
                   const fdparms_o8<T> &c);
cudaStream_t getGpuStream(void);

#endif

#if HIP_ENABLED
#include "hip/hip_runtime.h"
#include <hip/hip_complex.h>
void init_hip_fd(int max_threads, size_t bufsize);
void GpuFill(double *dptr, int n, double fillval);
void GpuNegate(double *dx, int incx, double *dy, int incy, int n);
void GpuEleMul(double *dx, std::complex<double> *dy, int n, hipStream_t stream);
void GpuEleMul(double *dx, std::complex<float> *dy, int n, hipStream_t stream);
template <typename T>
void app8_del2_gpu(T * __restrict__ a,
                   T *b,
                   const int dimx,
                   const int dimy,
                   const int dimz,
                   const fdparms_o8<T> &c);
hipStream_t getGpuStream(void);
#endif

#endif
