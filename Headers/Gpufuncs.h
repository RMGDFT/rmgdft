
#ifndef GPU_FUNCS_H
#define GPU_FUNCS_H 1

void GpuFill(double *dptr, int n, double fillval);
void GpuNegate(double *dx, int incx, double *dy, int incy, int n);

#endif
