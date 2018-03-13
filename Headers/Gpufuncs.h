
#ifndef GPU_FUNCS_H
#define GPU_FUNCS_H 1

#if GPU_ENABLED
#include <cublas_v2.h>
#include <cublasXt.h>

void GpuFill(double *dptr, int n, double fillval);
void GpuNegate(double *dx, int incx, double *dy, int incy, int n);
void gramsch_update_psi(double *V,
                        double *G,
                        double *C,
                        int N,
                        int eig_start,
                        int eig_stop,
                        cublasHandle_t cublasH);
#endif

#endif
