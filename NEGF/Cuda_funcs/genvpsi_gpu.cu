#include "arch.h"

#if GPU_ENABLED

__global__ void genvpsi_kernel(const double *psi, const double *veff, double *psiv, const int TILE_DIM,
        const int num_local_orbital, const int p0basis)
{


    int x = blockIdx.x * TILE_DIM + threadIdx.x;
    int y = blockIdx.y * TILE_DIM + threadIdx.y;

    if( x < num_local_orbital &&  y < p0basis)
        psiv[ x * p0basis + y] = psi [x * p0basis +y] * veff[y];
}



// C wrapper functions that call the cuda kernels above
extern "C" void genvpsi_gpu(const double *psi, const double *veff, double *psiv, 
        const int num_local_orbital, const int p0basis)
{



    int TILE_DIM = 32, BLOCK_ROWS = 8;
    dim3 dimGrid;
    dim3 dimBlock;
    dimGrid.x = (num_local_orbital + TILE_DIM -1)/TILE_DIM;
    dimGrid.y = (p0basis + TILE_DIM -1)/TILE_DIM;
    dimBlock.x = TILE_DIM;
    dimBlock.y = TILE_DIM;


    genvpsi_kernel<<<dimGrid, dimBlock>>>(psi, veff, psiv, TILE_DIM, num_local_orbital, p0basis);



}

#endif
