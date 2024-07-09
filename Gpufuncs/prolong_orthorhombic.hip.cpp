/*
 *
 * Copyright 2024 The RMG Project Developers. See the COPYRIGHT file 
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


#if HIP_ENABLED

#include <hip/hip_runtime.h>
#include <hip/hip_runtime.h>
#include <hip/hip_runtime_api.h>
#include <hip/hip_ext.h>
#include <stdio.h>
#include <unistd.h>
#include "Gpufuncs.h"
#include <iostream>
#include <vector>
//#include "Prolong.h"
#define MAX_PROLONG_RATIO 4
#define MAX_PROLONG_ORDER 12


//extern "C" void *BeginRmgTimer(const char *what);
//extern "C" void EndRmgTimer(void *ptr);

typedef struct {
  double a[MAX_PROLONG_RATIO][MAX_PROLONG_ORDER];
} pcoeff; 

template <typename T, int images>
void prolong_ortho_gpu(T *full,
                   T *half,
                   const int dimx,
                   const int dimy,
                   const int dimz,
                   double (&a)[MAX_PROLONG_RATIO][MAX_PROLONG_ORDER]);


// Version with shared memory slice
template <typename T, int images>
__global__ void prolong_ortho_kernel(T * full, 
                                     const T *half, 
                                     const int dimx,
                                     const int dimy,
                                     const int dimz,
                                     pcoeff a)

{
    extern __shared__ __align__(sizeof(T)) unsigned char sbuf[];
    T *fulla0 = reinterpret_cast<T *>(sbuf);
    T *fulla1 = fulla0 + ((dimy + 2*images) * (dimz + 2*images));
    T *fullb0 = fulla1 + ((dimy + 2*images) * (dimz + 2*images));
    T *fullb1 = fullb0 + ((2*dimy) * (dimz + 2*images));


    // iy and iz map to the x and y coordinates of the thread
    // within a block of the coarse grid array
    //const int iy1 = blockIdx.x * blockDim.x + threadIdx.x;
//    const int iy1 = blockIdx.y;
    //const int iz1 = blockIdx.y * blockDim.y + threadIdx.y;
    const int iz1 = threadIdx.y;

    const int incy = dimz + 2*images;
    const int incx = (dimy + 2*images) * incy;
    const int incy2 = 2*dimz;
    const int incx2 = 2*dimy * incy2;

    const int ic = images - 1;

    // lambda to clean up code
    auto stencil = [&](const T *ptr, const int stride) {
        T sum(0.0);
        for(int k = 0;k < 2*images;k++) sum = sum + a.a[1][k] * ptr[k*stride];
        return sum;
    };

    for (int ix = 0; ix < dimx; ix++)
    {
        // This first block is the only time we access global memory
        // The chunks are (dimy +2*images)*(dimz + 2*images) and
        // the stencil is applied along the z-direction as well.
        __syncthreads();
        for (int iy1 = 0; iy1 < dimy + 2*images; iy1++)
        {
            const T *halfptr = &half[(ix + 1) * incx + iy1 * incy];
            fulla0[iy1 * incy + iz1] = a.a[0][ic] * halfptr[ic*incx + iz1];
            fulla1[iy1 * incy + iz1] = stencil(halfptr + iz1, incx);
        }

        // This next block applies the stencil along the y-direction
        // y only has to go from (0,dimy) but we use the iy1 index
        // which goes to (dimy + 2*images) to avoid thread divergence
        //__syncthreads();
        for (int iy1 = 0; iy1 < dimy; iy1++)
        {
            T *full_tmp = &fulla0[(iy1 + 1) * incy + iz1];
            fullb0[(2 * iy1 + 0) * incy + iz1] = a.a[0][ic] * full_tmp[ic*incy];
            fullb0[(2 * iy1 + 1) * incy + iz1] = stencil(full_tmp, incy);

            full_tmp = &fulla1[(iy1 + 1) * incy + iz1];
            fullb1[(2 * iy1 + 0) * incy + iz1] = a.a[0][ic] * full_tmp[ic*incy];
            fullb1[(2 * iy1 + 1) * incy + iz1] = stencil(full_tmp, incy);
        }
        __syncthreads();

        // This last block writes back to global memory
        for (int iy1 = 0; iy1 < 2*dimy; iy1++)
        {
            if(iz1 < dimz)
            {
                T *full_tmp = &fullb0[iy1 * incy + iz1 + 1];
                full[2*ix * incx2 + iy1 * incy2 + 2 * iz1 + 0] += a.a[0][ic] * full_tmp[ic];
                full[2*ix * incx2 + iy1 * incy2 + 2 * iz1 + 1] += stencil(full_tmp, 1);

                full_tmp = &fullb1[iy1 * incy + iz1 + 1];
                full[(2*ix + 1) * incx2 + iy1 * incy2 + 2 * iz1 + 0] += a.a[0][ic] * full_tmp[ic];
                full[(2*ix + 1) * incx2 + iy1 * incy2 + 2 * iz1 + 1] += stencil(full_tmp, 1);
            }
        }
    }
    __syncthreads();

}
#endif

std::vector<void *> abufs;
std::vector<void *> hbufs;
std::vector<double *> rbufs;
std::vector<double> host_accumulator;

template <typename T>
void init_orthorhombic_gpu_prolong(int dimx, int dimy, int dimz)
{
    int order = MAX_PROLONG_ORDER;
    size_t rbufsize = 8*dimx*dimy*dimz*sizeof(double);
    size_t bufsize = (dimx + order)*(dimy + order)*(dimz + order)*sizeof(T);

    // Check if just clearing the accumulators
    if(rbufs.size() > 0)
    {
        GpuFill(rbufs[getThreadId()], 8*dimx*dimy*dimz, 0.0);
        std::fill(host_accumulator.begin(), host_accumulator.end(), 0.0);
        return;
    }

    int max_threads = getThreadNum();
    abufs.resize(max_threads);
    rbufs.resize(max_threads);
    hbufs.resize(max_threads);
    host_accumulator.resize(8*dimx*dimy*dimz);
    for(int i=0;i < max_threads;i++)
    {
        hipMalloc((void **)&abufs[i], bufsize);
        hipMallocHost((void **)&hbufs[i], rbufsize);
        hipMalloc((void **)&rbufs[i], rbufsize);
    }
    hipDeviceSynchronize();
    GpuFill(rbufs[getThreadId()], 8*dimx*dimy*dimz, 0.0);
    std::fill(host_accumulator.begin(), host_accumulator.end(), 0.0);
}

void fetch_gpu_density(double *rho)
{
    hipStream_t stream = getGpuStream();
    int numthreads = getThreadNum();
    int fbasis = host_accumulator.size();
    for(int tid=0;tid < numthreads;tid++)
    {
        hipMemcpyAsync(hbufs[tid], rbufs[tid], fbasis*sizeof(double), hipMemcpyDeviceToHost, stream);
        hipStreamSynchronize(stream);
        double *hptr = (double *)hbufs[tid];
        for(int idx=0;idx<fbasis;idx++) host_accumulator[idx] += hptr[idx];
    }
    for(int idx=0;idx<fbasis;idx++) rho[idx] = host_accumulator[idx];
}

template void init_orthorhombic_gpu_prolong<double>(int, int, int);
template void init_orthorhombic_gpu_prolong<std::complex<double>>(int, int, int);
template void init_orthorhombic_gpu_prolong<float>(int, int, int);
template void init_orthorhombic_gpu_prolong<std::complex<float>>(int, int, int);


template void prolong_ortho_gpu<double,4>(double * , double *, int, int, int, double (&)[MAX_PROLONG_RATIO][MAX_PROLONG_ORDER]);
template void prolong_ortho_gpu<std::complex<double>,4>(std::complex<double> * , std::complex<double> *, int, int, int, double (&)[MAX_PROLONG_RATIO][MAX_PROLONG_ORDER]);

template void prolong_ortho_gpu<double,5>(double * , double *, int, int, int, double (&)[MAX_PROLONG_RATIO][MAX_PROLONG_ORDER]);
template void prolong_ortho_gpu<std::complex<double>,5>(std::complex<double> * , std::complex<double> *, int, int, int, double (&)[MAX_PROLONG_RATIO][MAX_PROLONG_ORDER]);

template void prolong_ortho_gpu<double,6>(double * , double *, int, int, int, double (&)[MAX_PROLONG_RATIO][MAX_PROLONG_ORDER]);
template void prolong_ortho_gpu<std::complex<double>,6>(std::complex<double> * , std::complex<double> *, int, int, int, double (&)[MAX_PROLONG_RATIO][MAX_PROLONG_ORDER]);


void prolong_test_gpu(double *full, double *half, int dimx, int dimy, int dimz, double (&a)[MAX_PROLONG_RATIO][MAX_PROLONG_ORDER])
{
  prolong_ortho_gpu<double,6>(full, half, dimx, dimy, dimz, a);
}
#include "RmgTimer.h"
template <typename T, int images>
void prolong_ortho_gpu(T *full, 
                   T *half, 
                   const int dimx,
                   const int dimy,
                   const int dimz,
                   double (&a)[MAX_PROLONG_RATIO][MAX_PROLONG_ORDER])
{

    pcoeff agpu;
    for(int i=0;i<MAX_PROLONG_RATIO;i++)
    {
        for(int j=0;j<MAX_PROLONG_ORDER;j++)
        {
            agpu.a[i][j] = a[i][j];
        }
    }


    dim3 Grid, Block;
    hipStream_t stream = getGpuStream();
    int tid = getThreadId();
    int fbasis = 8*dimx*dimy*dimz;
    int sbasis = (dimx+2*images)*(dimy+2*images)*(dimz+2*images);


    Grid.x = 1;
    Grid.y = 1;
    Block.x = 1;
    Block.y = dimz + 2*images;
    int smem_siz = 2*(dimy + 2*images) * (dimz + 2*images) +
                   2*(2*dimy) * (dimz + 2*images);
    smem_siz *= sizeof(T);

    T *hptr = (T *)hbufs[tid];
    for(int idx=0;idx<sbasis;idx++) hptr[idx] = half[idx];
    hipStreamSynchronize(stream);
    hipMemcpyAsync(abufs[tid], hbufs[tid], sbasis*sizeof(T), hipMemcpyHostToDevice, stream);

    hipLaunchKernelGGL(HIP_KERNEL_NAME(prolong_ortho_kernel<T, images>), Grid, Block, smem_siz, stream,
               (T *)rbufs[tid], (T *)abufs[tid], dimx, dimy, dimz, agpu);

    hipStreamSynchronize(stream);
RmgTimer RT("Copy density");
    hipMemcpyAsync(hbufs[tid], rbufs[tid], fbasis*sizeof(T), hipMemcpyDeviceToHost, stream);
    hipStreamSynchronize(stream);
    for(int idx=0;idx<fbasis;idx++) full[idx] = hptr[idx];

}
