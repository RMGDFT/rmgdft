#include <iostream>     
#include <algorithm>    
#include <cfloat>
#include <math.h>       
#include <mpi.h>       
#include "RmgException.h"
#include "transition.h"
#include "GlobalSums.h"
#include "Atomic.h"

#if HIP_ENABLED
#include <hip/hip_runtime.h>
#endif

template void  TiledM_to_glob<double>(double *matrix_glob, double *tiledM, int numst, MPI_Comm tiled_comm);
template void  TiledM_to_glob<std::complex<double>>(std::complex<double> *matrix_glob, std::complex<double> *tiledM, int numst, MPI_Comm tiled_comm);

template <typename T> void  TiledM_to_glob(T *matrix_glob, T *tiledM, int numst, MPI_Comm tiled_comm)
{
    int nprocs, my_rank;
    MPI_Comm_size(tiled_comm, &nprocs);
    MPI_Comm_rank(tiled_comm, &my_rank);
    if(numst%nprocs != 0) 
    {
        rmg_printf("\n numst nprocs %d %d", numst, nprocs);
        fflush(NULL);
        rmg_error_handler(__FILE__, __LINE__, "numst must be divisible by nprocs ");
    }

    bool matrix_glob_dev = false;
    bool tiledM_dev = false;
#if CUDA_ENABLED
    cudaPointerAttributes attr;
    cudaError_t cudaerr;
    cudaerr = cudaPointerGetAttributes(&attr, matrix_glob);
    if(cudaerr == cudaSuccess && attr.type == cudaMemoryTypeDevice) matrix_glob_dev = true;
    cudaerr = cudaPointerGetAttributes(&attr, tiledM);
    if(cudaerr == cudaSuccess && attr.type == cudaMemoryTypeDevice) tiledM_dev = true;
#elif HIP_ENABLED
    hipPointerAttribute_t attr;
    hipError_t hiperr;
    hiperr = hipPointerGetAttributes(&attr, matrix_glob);
    if(hiperr == hipSuccess && attr.type == hipMemoryTypeDevice) matrix_glob_dev = true;
    hiperr = hipPointerGetAttributes(&attr, tiledM);
    if(hiperr == hipSuccess && attr.type == hipMemoryTypeDevice) tiledM_dev = true;
#endif

    if(tiledM_dev != matrix_glob_dev)
    {
        std::cout<< tiledM_dev << matrix_glob_dev << std::endl;
        rmg_error_handler(__FILE__, __LINE__, "both matrix must be in device or both in host  ");
    }

    size_t sendcount = numst * numst/nprocs * sizeof(T)/sizeof(double);
    if(!tiledM_dev)
    {
        MPI_Allgather(tiledM, sendcount, MPI_DOUBLE, matrix_glob, sendcount, MPI_DOUBLE, tiled_comm);
    }
    else
    {
#if (CUDA_ENABLED || HIP_ENABLED) && USE_NCCL
        ncclAllGather(tiledM, matrix_glob, sendcount, ncclDouble, ct.nccl_local_comm, 0);
//keep those if we have trouble in the future
//#elif (CUDA_ENABLED) && USE_NCCL
//   use allreduce instead of allgather
//    size_t matrix_size = numst * numst * sizeof(T);
//    cudaMemset(matrix_glob, 0, matrix_size);
//    matrix_size = numst * numst/nprocs * sizeof(T);
//    RmgMemcpy(&matrix_glob[my_rank * numst * numst/nprocs], tiledM, matrix_size );
//    sendcount = numst * numst * sizeof(T)/sizeof(double);
//    ncclAllReduce(matrix_glob, matrix_glob, sendcount, ncclDouble, ncclSum, ct.nccl_local_comm, 0);

  // use cpu MPI
  //std::complex<double> *matrix_cpu = new std::complex<double>[numst*numst];
  //std::complex<double> *tm_cpu = new std::complex<double>[numst*numst/nprocs];

  //size_t matrix_size = numst * numst/nprocs * sizeof(std::complex<double>);
  //RmgMemcpy(tm_cpu, tiledM, matrix_size );
  //MPI_Allgather(tm_cpu, sendcount, MPI_DOUBLE, matrix_cpu, sendcount, MPI_DOUBLE, tiled_comm);
  //matrix_size = numst * numst * sizeof(std::complex<double>);
  //RmgMemcpy(matrix_glob, matrix_cpu, matrix_size);
  //delete [] matrix_cpu;
  //delete [] tm_cpu;

#else
	rmg_error_handler(__FILE__, __LINE__, "compile with -DUSE_NCCL=1  ");
#endif
    }

}
