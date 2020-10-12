#include "negf_prototypes.h"

/************************** SVN Revision Information **************************
 **    $Id$    **
 ******************************************************************************/


#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <complex>


#if CUDA_ENABLED
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cublas_v2.h>
#endif



void setvector_host_device (int n, int elemsize, void *x_host, int ia, void *x_device, int ib)
{

#if CUDA_ENABLED
    cublasSetVector( n, elemsize, x_host, ia, x_device, ib );
#endif 
}

void getvector_device_host (int n, int elemsize, void *x_device, int ia, void *x_host, int ib)
{

#if CUDA_ENABLED
    cublasGetVector( n, elemsize, x_device, ia, x_host, ib );
#endif 
}

double *memory_ptr_host_device(double *ptr_host, double *ptr_device)
{
#if CUDA_ENABLED
    return ptr_device;
#else
    return ptr_host;
#endif 
}

std::complex<double> *memory_ptr_host_device(std::complex<double> *ptr_host, std::complex<double> *ptr_device)
{
#if CUDA_ENABLED
    return ptr_device;
#else
    return ptr_host;
#endif 
}
