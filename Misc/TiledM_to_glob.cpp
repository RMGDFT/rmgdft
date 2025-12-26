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
#else
        rmg_error_handler(__FILE__, __LINE__, "compile with -DUSE_NCCL=1  ");
#endif
    }

}
