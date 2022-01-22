
#ifndef GPU_FUNCS_H
#define GPU_FUNCS_H 1

template <typename T>
struct fdparms_o8 {
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

#if CUDA_ENABLED
#include <cublas_v2.h>
#include <complex>

void GpuFill(double *dptr, int n, double fillval);
void GpuNegate(double *dx, int incx, double *dy, int incy, int n);
void gramsch_update_psi(double *V,
                        double *C,
                        int N,
                        int eig_start,
                        int eig_stop,
                        cublasHandle_t cublasH);
template <typename T>
double app8_del2_gpu(const T * __restrict__ a,
                   T *b,
                   const int dimx,
                   const int dimy,
                   const int dimz,
                   T h2x,
                   T h2y,
                   T h2z,
                   cudaStream_t cstream);
void GpuEleMul(double *dx, std::complex<double> *dy, int n, cudaStream_t stream);
void GpuEleMul(double *dx, std::complex<float> *dy, int n, cudaStream_t stream);

#endif

#if HIP_ENABLED
#include "hip/hip_runtime.h"
#include <hip/hip_complex.h>
void GpuFill(double *dptr, int n, double fillval);
void GpuNegate(double *dx, int incx, double *dy, int incy, int n);
void GpuEleMul(double *dx, std::complex<double> *dy, int n, hipStream_t stream);
void GpuEleMul(double *dx, std::complex<float> *dy, int n, hipStream_t stream);
template <typename T>
double app8_del2_gpu(const T * __restrict__ a,
                   T *b,
                   const int dimx,
                   const int dimy,
                   const int dimz,
                   double h2x,
                   double h2y,
                   double h2z,
                   hipStream_t cstream);
#endif

#endif
