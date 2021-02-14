
#ifndef GPU_FUNCS_H
#define GPU_FUNCS_H 1

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
double app8_del2_gpu(const double * __restrict__ a,
                   double *b,
                   const int dimx,
                   const int dimy,
                   const int dimz,
                   double h2x,
                   double h2y,
                   double h2z,
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
#endif

#endif
