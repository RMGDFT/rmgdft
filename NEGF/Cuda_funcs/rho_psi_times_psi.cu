#include "arch.h"

#if GPU_ENABLED

#define THREAD_PER_BLOCK 32
#define THREAD_PER_BLOCKX 32

__device__ inline void atomicAdd_double(double* address, double val)
{
    unsigned long long int* address_as_ull =
                              (unsigned long long int*)address;
    unsigned long long int old = *address_as_ull, assumed;
    do {
        assumed = old;
        old = atomicCAS(address_as_ull, assumed,
                        __double_as_longlong(val +
                               __longlong_as_double(assumed)));
    } while (assumed != old);
}

__global__ void rho_psi_times_psi_kernel(const double *psi1, const double *psi2, double *rho, 
        const int num_local_orbital, const int p0basis)
{


    int x = blockIdx.x * THREAD_PER_BLOCKX +threadIdx.x ;
    int y = blockIdx.y * THREAD_PER_BLOCK + threadIdx.y;

    __shared__ double temp[THREAD_PER_BLOCK * THREAD_PER_BLOCKX];

    temp[threadIdx.y * THREAD_PER_BLOCKX + threadIdx.x] = 0.0;
    if( y < num_local_orbital &&  x < p0basis)
        temp[threadIdx.y * THREAD_PER_BLOCKX + threadIdx.x] = psi1[ y * p0basis + x] * psi2 [y * p0basis +x];
    __syncthreads();

    if( threadIdx.y == 0)
    {
        double sum = 0.0;
        for(int i = 0; i < THREAD_PER_BLOCK;  i++) sum += temp[i * THREAD_PER_BLOCKX + threadIdx.x];
        if(x < p0basis) atomicAdd_double(&rho[x], sum);
    }

}



// C wrapper functions that call the cuda kernels above
extern "C" void rho_psi_times_psi(const double *psi1, const double *psi2, double *rho, 
        const int num_local_orbital, const int p0basis)
{



    dim3 dimGrid;
    dim3 dimBlock;
    dimGrid.x = (p0basis + THREAD_PER_BLOCKX -1)/ THREAD_PER_BLOCKX;
    dimGrid.y = (num_local_orbital + THREAD_PER_BLOCK-1)/THREAD_PER_BLOCK;
    dimBlock.x = THREAD_PER_BLOCKX;
    dimBlock.y = THREAD_PER_BLOCK;


    rho_psi_times_psi_kernel<<<dimGrid, dimBlock>>>(psi1, psi2, rho, num_local_orbital, p0basis);



}


#endif


