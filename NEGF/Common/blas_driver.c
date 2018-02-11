/************************** SVN Revision Information **************************
 **    $Id$    **
 ******************************************************************************/


#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>


#if GPU_ENABLED
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cublas_v2.h>
#endif

#include "main.h"
#include "Scalapack.h"
#include "blas.h"



void zcopy_driver (int n, complex double *A, int ia, complex double *B, int ib) 
{

#if GPU_ENABLED
    cublasZcopy (ct.cublas_handle, n, A, ia, B, ib);
#else
    zcopy (&n, A, &ia, B, &ib);
#endif
}


void zaxpy_driver (int n, complex double alpha, complex double *A, int ia, complex double *B, int ib) 
{

#if GPU_ENABLED
    cublasZaxpy (ct.cublas_handle, n, &alpha, A, ia, B, ib);
#else
    zaxpy (&n, &alpha, A, &ia, B, &ib);
#endif
}

dzasum_driver(int n, complex double *A, int ia, double *sum)
{
#if GPU_ENABLED
    cublasDzasum (ct.cublas_handle, n, A, ia, sum);
#else
    *sum = dzasum(&n, (double *)A, &ia);
#endif
}


