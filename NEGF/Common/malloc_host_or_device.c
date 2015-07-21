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




void *malloc_host_or_device (size_t size)

{

    void *ptr;

#if GPU_ENABLED
    cudaError_t cuerr;
    cuerr = cudaMalloc((void **)&ptr , size);
    if(cuerr != cudaSuccess)
    {
        printf ("cuda allocation faild in malloc_host_or_device= %d \n", size);
        fflush (NULL);
        exit (0);
    }
    return ptr;

#else
    if(NULL == (ptr = malloc(size)))
    {
        printf ("allocation faild in malloc_host_or_device= %d \n", size);
        fflush (NULL);
        exit (0);
    }

    return ptr;
#endif

}

void free_host_or_device(void *ptr)
{
#if GPU_ENABLED
    cudaFree(ptr);
#else
    free(ptr);
#endif
}

