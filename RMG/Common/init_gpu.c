/*
 *
 * Copyright 2014 The RMG Project Developers. See the COPYRIGHT file 
 * at the top-level directory of this distribution or in the current
 * directory.
 * 
 * This file is part of RMG. 
 * RMG is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * any later version.
 *
 * RMG is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
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
        #include <magma.h>
    #endif

void rmg_printout_devices( )
{
    int ndevices, idevice;
    cuDeviceGetCount( &ndevices );
    char name[200];
    size_t totalMem;

    int clock;
    CUdevice dev;

    for(idevice = 0; idevice < ndevices; idevice++ ) {

        cuDeviceGet( &dev, idevice );
        cuDeviceGetName( name, sizeof(name), dev );
        cuDeviceTotalMem( &totalMem, dev );
        ct.gpu_mem[idevice] = totalMem;
        cuDeviceGetAttribute( &clock, CU_DEVICE_ATTRIBUTE_CLOCK_RATE, dev );
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

//  cublasStatus_t custat;

  rmg_printout_devices( );
  cudaSetDeviceFlags(cudaDeviceScheduleSpin);
  // Create a stream for the main process
  //cudaStreamCreate(&ct.cuda_stream);
  //cublasSetStream(ct.cublas_handle, ct.cuda_stream); 

}

void finalize_gpu (void)
{
    cublasDestroy(ct.cublas_handle);
    cublasXtDestroy(ct.cublasXt_handle);

//    cuCtxDetach( ct.cu_context ); 
 //   cublasShutdown();

}

#endif
