/************************** SVN Revision Information **************************
 **    $Id$    **
 ******************************************************************************/


#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include "make_conf.h"

#if GPU_ENABLED
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cublas_v2.h>
#endif

#include "main.h"
#include "my_scalapack.h"
#include "blas.h"



void zcopy_driver (int n, complex double *A, int ia, complex double *B, int ib) 
{

#if GPU_ENABLED
    cublasZcopy (ct.cublas_handle, n, A, ia, B, ib);
#else
    zcopy (&n, A, &ia, B, &ib);
#endif
}
