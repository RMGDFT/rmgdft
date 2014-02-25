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
#include "init_var.h"
#include "LCR.h"
#include "pmo.h"

#if GPU_ENABLED

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
}



cudaStream_t rmg_default_cuda_stream;
void init_gpu (void)
{

	cublasStatus_t custat;
	int alloc;
	int i, ntot_row, maxrow;

	rmg_printout_devices( );

	if(alloc < 1024 * 1024) alloc = 1024 * 1024;
	if(alloc < pct.num_local_orbit * get_P0_BASIS()) alloc = pct.num_local_orbit * get_P0_BASIS();

	if( cudaSuccess != cudaMallocHost((void **)&ct.gpu_states, alloc * sizeof(rmg_double_t) )){
		fprintf (stderr, "Error: cudaMallocHost failed for: ct.gpu_states\n");
		exit(-1);
	}

	alloc = pmo.ntot * sizeof(complex double);
	if( cudaSuccess != cudaMalloc((void **)&ct.gpu_Htri , alloc )){
		fprintf (stderr, "!!!! cublasAlloc failed for: gpu_global_matrix\n");
		exit(-1);
	}
	if( cudaSuccess != cudaMalloc((void **)&ct.gpu_Gtri , alloc )){
		fprintf (stderr, "!!!! cublasAlloc failed for: gpu_global_matrix\n");
		exit(-1);
	}

	ntot_row = 0;
	maxrow = 0;
	for(i = 0; i < ct.num_blocks; i++)
	{
		ntot_row += ct.block_dim[i];
		maxrow = max(maxrow, ct.block_dim[i]);
	}

	alloc = ntot_row * maxrow * sizeof(complex double);
	if( cudaSuccess != cudaMalloc((void **)&ct.gpu_Gtem , alloc ) ){
		fprintf (stderr, "Error: cudaMalloc failed for: gpu_Gtem\n");
		exit(-1);
	}

	alloc = maxrow * maxrow * sizeof(complex double);
	if( cudaSuccess != cudaMalloc((void **)&ct.gpu_Hii, alloc )){
		fprintf (stderr, "Error: cudaMalloc failed for: ct.gpu_Hii\n");
		exit(-1);
	}
	if( cudaSuccess != cudaMalloc((void **)&ct.gpu_Gii, alloc )){
		fprintf (stderr, "Error: cudaMalloc failed for: ct.gpu_Gii\n");
		exit(-1);
	}
	if( cudaSuccess != cudaMalloc((void **)&ct.gpu_temp, alloc )){
		fprintf (stderr, "Error: cudaMalloc failed for: ct.gpu_temp\n");
		exit(-1);
	}
	if( cudaSuccess != cudaMalloc((void **)&ct.gpu_Imatrix, alloc )){
		fprintf (stderr, "Error: cudaMalloc failed for: ct.gpu_Imatrix\n");
		exit(-1);
	}

	alloc = maxrow * sizeof(int);
	if( cudaSuccess != cudaMalloc((void **)&ct.gpu_ipiv, alloc )){
		fprintf (stderr, "Error: cudaMalloc failed for: ct.gpu_ipiv\n");
		exit(-1);
	}

	// Make sure enough is allocated for hartree solver on fine grid
	alloc = 1024 ;
	if(alloc < pct.num_local_orbit * get_P0_BASIS()) alloc = pct.num_local_orbit * get_P0_BASIS();

	if( cudaSuccess != cudaMallocHost((void **)&ct.gpu_host_temp1, alloc * sizeof(rmg_double_t) )){
		fprintf (stderr, "Error: cudaMallocHost failed for: ct.gpu_host_temp1\n");
		exit(-1);
	}

	alloc = 1024 ;
	if(alloc < get_P0_BASIS()) alloc = get_P0_BASIS();
	if(alloc < ct.num_states * ct.num_states) alloc = ct.num_states * ct.num_states;

	if( cudaSuccess != cudaMallocHost((void **)&ct.gpu_host_temp2, alloc * sizeof(rmg_double_t) )){
		fprintf (stderr, "Error: cudaMallocHost failed for: ct.gpu_host_temp2\n");
		exit(-1);
	}

	// fdbuf needs to be big enough for hartee potential on fine grid
	alloc = ct.THREADS_PER_NODE * (get_PX0_GRID() + 4) * (get_PY0_GRID() + 4) * (get_PZ0_GRID() + 4);
	if(alloc < ((get_FPX0_GRID() + 2)*(get_FPY0_GRID() + 2)*(get_FPZ0_GRID() + 2))) alloc = (get_FPX0_GRID() + 2)*(get_FPY0_GRID() + 2)*(get_FPZ0_GRID() + 2);
	if( cudaSuccess != cudaMallocHost((void **)&ct.gpu_host_fdbuf1, alloc * sizeof(rmg_double_t) )){
		fprintf (stderr, "Error: cudaMallocHost failed for: ct.gpu_host_fdbuf1\n");
		exit(-1);
	}
	if( cudaSuccess != cudaMallocHost((void **)&ct.gpu_host_fdbuf2, alloc * sizeof(rmg_double_t) )){
		fprintf (stderr, "Error: cudaMallocHost failed for: ct.gpu_host_fdbuf2\n");
		exit(-1);
	}
	if( cudaSuccess != cudaMallocHost((void **)&ct.gpu_host_work, (3 * ct.num_states*ct.num_states + 8*ct.num_states) * sizeof(rmg_double_t) )){
		fprintf (stderr, "Error: cudaMallocHost failed for: ct.gpu_host_temp\n");
		exit(-1);
	}



	custat = cublasCreate(&ct.cublas_handle);
	if( custat != CUBLAS_STATUS_SUCCESS ) {
		fprintf (stderr, "Error cublasCreate failed for: ct.cublas_handle\n");
		exit(-1);
	}

	cudaSetDeviceFlags(cudaDeviceScheduleSpin);
	// Create a stream for the main process
	cudaStreamCreate(&ct.cuda_stream);
	//  cublasSetStream(ct.cublas_handle, ct.cuda_stream); 
	//  magmablasSetKernelStream(ct.cuda_stream);

}

void finalize_gpu (void)
{

	cublasDestroy(ct.cublas_handle);
	cudaFree(ct.gpu_global_matrix);
	cudaFree(ct.gpu_Gtem);
	cudaFree(ct.gpu_work1);
	cudaFree(ct.gpu_work2);
	cudaFree(ct.gpu_work3);
	cudaFree(ct.gpu_work4);
	cudaFreeHost(ct.gpu_host_temp4);
	cudaFreeHost(ct.gpu_host_temp3);
	cudaFreeHost(ct.gpu_host_temp2);
	cudaFreeHost(ct.gpu_host_temp1);

	//    cuCtxDetach( ct.cu_context ); 
	//   cublasShutdown();

}

#endif
