#include <iostream>     
#include <algorithm>    
#include <cfloat>
#include <math.h>       
#include <mpi.h>       
#include "RmgException.h"
#include "transition.h"
#include "GlobalSums.h"
#include "Atomic.h"

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

    size_t sendcount = numst * numst/nprocs * sizeof(T)/sizeof(double);
#if (CUDA_ENABLED || HIP_ENABLED) && USE_NCCL
    ncclAllGather(tiledM, matrix_glob, sendcount, ncclDouble, ct.nccl_local_comm, 0);
#else
    MPI_Allgather(tiledM, sendcount, MPI_DOUBLE, matrix_glob, sendcount, MPI_DOUBLE, tiled_comm);
#endif

}
