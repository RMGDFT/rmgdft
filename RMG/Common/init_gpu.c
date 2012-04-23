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

#if GPU_ENABLED

void printout_devices( )
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
}




void init_gpu (void)
{

  cublasStatus_t custat;

  cudaDeviceReset();
  if( CUDA_SUCCESS != cuInit( 0 ) ) {
      fprintf(stderr, "CUDA: Not initialized\n" ); exit(-1);
  }
  if( CUDA_SUCCESS != cuDeviceGet( &ct.cu_dev, 0 ) ) {
      fprintf(stderr, "CUDA: Cannot get the device\n"); exit(-1);
  }
//  if( CUDA_SUCCESS != cuCtxCreate( &ct.cu_context, 0, ct.cu_dev ) ) {
//      fprintf(stderr, "CUDA: Cannot create the context\n"); exit(-1);
//  }
//  if( CUBLAS_STATUS_SUCCESS != cublasInit( ) ) {
//      fprintf(stderr, "CUBLAS: Not initialized\n"); exit(-1);
//  }
  printout_devices( );

  if( cudaSuccess != cudaMalloc((void **)&ct.gpu_global_matrix , ct.num_states * ct.num_states * sizeof(REAL) )){
      fprintf (stderr, "!!!! cublasAlloc failed for: gpu_global_matrix\n");
      exit(-1);
  }
  if( cudaSuccess != cudaMalloc((void **)&ct.gpu_states , ct.num_states * P0_BASIS * sizeof(REAL) )){
      fprintf (stderr, "Error: cudaMalloc failed for: gpu_states\n");
      exit(-1);
  }
  if( cudaSuccess != cudaMalloc((void **)&ct.gpu_temp , ct.num_states * P0_BASIS * sizeof(REAL) )){
      fprintf (stderr, "Error: cudaMalloc failed for: gpu_temp\n");
      exit(-1);
  }
  if( cudaSuccess != cudaMallocHost((void **)&ct.gpu_host_temp1, THREADS_PER_NODE * (PX0_GRID + 4) * (PY0_GRID + 4) * (PZ0_GRID + 4) * sizeof(REAL) )){
      fprintf (stderr, "Error: cudaMallocHost failed for: ct.gpu_host_temp\n");
      exit(-1);
  }
  if( cudaSuccess != cudaMallocHost((void **)&ct.gpu_host_temp2, 4 * THREADS_PER_NODE * (PX0_GRID + 4) * (PY0_GRID + 4) * (PZ0_GRID + 4) * sizeof(REAL) )){
      fprintf (stderr, "Error: cudaMallocHost failed for: ct.gpu_host_temp\n");
      exit(-1);
  }
  if( cudaSuccess != cudaMallocHost((void **)&ct.gpu_host_temp3, THREADS_PER_NODE * (PX0_GRID + 4) * (PY0_GRID + 4) * (PZ0_GRID + 4) * sizeof(REAL) )){
      fprintf (stderr, "Error: cudaMallocHost failed for: ct.gpu_host_temp\n");
      exit(-1);
  }
  if( cudaSuccess != cudaMallocHost((void **)&ct.gpu_host_temp4, THREADS_PER_NODE * (PX0_GRID + 4) * (PY0_GRID + 4) * (PZ0_GRID + 4) * sizeof(REAL) )){
      fprintf (stderr, "Error: cudaMallocHost failed for: ct.gpu_host_temp\n");
      exit(-1);
  }
  if( cudaSuccess != cudaMalloc((void **)&ct.gpu_work1, THREADS_PER_NODE * (PX0_GRID + 4) * (PY0_GRID + 4) * (PZ0_GRID + 4) * sizeof(REAL) )){
      fprintf (stderr, "Error: cudaMalloc failed for: ct.gpu_work\n");
      exit(-1);
  }
  if( cudaSuccess != cudaMalloc((void **)&ct.gpu_work2, THREADS_PER_NODE * (PX0_GRID + 4) * (PY0_GRID + 4) * (PZ0_GRID + 4) * sizeof(REAL) )){
      fprintf (stderr, "Error: cudaMalloc failed for: ct.gpu_work\n");
      exit(-1);
  }
  if( cudaSuccess != cudaMalloc((void **)&ct.gpu_work3, THREADS_PER_NODE * (PX0_GRID + 4) * (PY0_GRID + 4) * (PZ0_GRID + 4) * sizeof(REAL) )){
      fprintf (stderr, "Error: cudaMalloc failed for: ct.gpu_work\n");
      exit(-1);
  }
  if( cudaSuccess != cudaMalloc((void **)&ct.gpu_work4, THREADS_PER_NODE * (PX0_GRID + 4) * (PY0_GRID + 4) * (PZ0_GRID + 4) * sizeof(REAL) )){
      fprintf (stderr, "Error: cudaMalloc failed for: ct.gpu_work\n");
      exit(-1);
  }

  custat = cublasCreate(&ct.cublas_handle);
  if( custat != CUBLAS_STATUS_SUCCESS ) {
      fprintf (stderr, "Error cublasCreate failed for: ct.cublas_handle\n");
      exit(-1);
  }

}

void finalize_gpu (void)
{
 
    cublasDestroy(ct.cublas_handle);
    cudaFree(ct.gpu_global_matrix);
    cudaFree(ct.gpu_temp);
    cudaFree(ct.gpu_states);
    cudaFreeHost(ct.gpu_host_temp4);
    cudaFreeHost(ct.gpu_host_temp3);
    cudaFreeHost(ct.gpu_host_temp2);
    cudaFreeHost(ct.gpu_host_temp1);

//    cuCtxDetach( ct.cu_context ); 
 //   cublasShutdown();

}

#endif
