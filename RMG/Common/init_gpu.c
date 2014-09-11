/*   init_gpu.c
 * COPYRIGHT
 *   Copyright (C) 1995  Emil Briggs
 *   Copyright (C) 1998  Emil Briggs, Charles Brabec, Mark Wensell, 
 *                       Dan Sullivan, Chris Rapcewicz, Jerzy Bernholc
 *   Copyright (C) 2001  Emil Briggs, Wenchang Lu,
 *                       Marco Buongiorno Nardelli,Charles Brabec, 
 *                       Mark Wensell,Dan Sullivan, Chris Rapcewicz,
 *                       Jerzy Bernholc
 * FUNCTION
 *   void init_gpu(void)
 *   sets up gpu
 * INPUTS
 *   nothing
 * OUTPUT
 *   nothing
 * PARENTS
 *   init.c
 * CHILDREN
 *   nothing
 * SOURCE
 */



#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "main.h"
#include "common_prototypes.h"

#if GPU_ENABLED
    #include <cuda.h>
    #include <cuda_runtime_api.h>
    #include <cublas_v2.h>

    #if MAGMA_LIBS
        #include <magma/magma.h>
    #endif

void rmg_printout_devices( )
{
    int ndevices, idevice;
    cuDeviceGetCount( &ndevices );
    char name[200];
#if CUDA_VERSION > 3010 
    size_t totalMem;
#else
    unsigned int totalMem;
#endif

    int clock;
    CUdevice dev;

    for(idevice = 0; idevice < ndevices; idevice++ ) {

        cuDeviceGet( &dev, idevice );
        cuDeviceGetName( name, sizeof(name), dev );
        cuDeviceTotalMem( &totalMem, dev );
        cuDeviceGetAttribute( &clock,
        CU_DEVICE_ATTRIBUTE_CLOCK_RATE, dev );
        printf( "device %d: %s, %.1f MHz clock, %.1f MB memory\n",
        idevice, name, clock/1000.f, totalMem/1024.f/1024.f );
    }
#if MAGMA_LIBS
    magma_init();          
#endif
}



//cudaStream_t rmg_default_cuda_stream;
void init_gpu (void)
{

  cublasStatus_t custat;
  int alloc;


  rmg_printout_devices( );
  cudaSetDeviceFlags(cudaDeviceScheduleSpin);
  // Create a stream for the main process
  cudaStreamCreate(&ct.cuda_stream);
//  cublasSetStream(ct.cublas_handle, ct.cuda_stream); 
//  magmablasSetKernelStream(ct.cuda_stream);

}

void finalize_gpu (void)
{
 
    cublasDestroy(ct.cublas_handle);
    cudaFreeHost(ct.gpu_host_temp4);
    cudaFreeHost(ct.gpu_host_temp3);
    cudaFreeHost(ct.gpu_host_temp2);
    cudaFreeHost(ct.gpu_host_temp1);

//    cuCtxDetach( ct.cu_context ); 
 //   cublasShutdown();

}

#endif
