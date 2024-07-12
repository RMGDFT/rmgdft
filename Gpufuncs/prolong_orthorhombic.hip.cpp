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
#include "Prolong.h"
#include "const.h"
#include "params.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "rmg_control.h"


typedef struct {
  float a[MAX_PROLONG_RATIO][MAX_PROLONG_ORDER];
} pcoeff; 

template <typename T, int ord>
void prolong_ortho_gpu(T *full,
                   T *half,
                   const int dimx,
                   const int dimy,
                   const int dimz);

template <typename T, int images>
void prolong_ortho_gpu_internal(T *full,
                   T *half,
                   const int dimx,
                   const int dimy,
                   const int dimz,
                   float (&a)[MAX_PROLONG_RATIO][MAX_PROLONG_ORDER]);


// Version with shared memory slice
template <typename T, int images>
__global__ void prolong_ortho_kernel(T * full, 
                                     const T *half, 
                                     const int zstart,
                                     const int zlen,
                                     const int dimx,
                                     const int dimy,
                                     const int dimz,
                                     int tid,
                                     pcoeff a)

{
    extern __shared__ __align__(sizeof(T)) unsigned char sbuf[];
    T *fulla0 = reinterpret_cast<T *>(sbuf);
    T *fulla1 = fulla0 + ((dimy + 2*images) * (zlen + 2*images));
    T *fullb0 = fulla1 + ((dimy + 2*images) * (zlen + 2*images));
    T *fullb1 = fullb0 + ((2*dimy) * (zlen + 2*images));

    const int iz1 = threadIdx.y;

    const int incy = dimz + 2*images;
    const int incx = (dimy + 2*images) * incy;
    const int sincy = zlen + 2*images;
    const int sincx = (dimy + 2*images) * sincy;

    const int incy2 = 2*dimz;
    const int incx2 = 2*dimy * incy2;

    const int ic = images - 1;

    // lambda to clean up code
    auto stencil = [&](const T *ptr, const int stride) {
        T sum(0.0);
        for(int k = 0;k < 2*images;k++) sum = sum + a.a[1][k] * ptr[k*stride];
        return sum;
    };

    for(int ix = blockIdx.x * blockDim.x + threadIdx.x;ix < dimx;ix += gridDim.x * blockDim.x)
    {
        // This first block is the only time we access global memory
        // The chunks are (dimy +2*images)*(dimz + 2*images) and
        // the stencil is applied along the z-direction as well.
        for (int iy1 = 0; iy1 < dimy + 2*images; iy1++)
        {
            __syncthreads();
            for(int iz1 = blockIdx.y*blockDim.y + threadIdx.y;iz1 < zlen+2*images;iz1 += gridDim.y*blockDim.y)
            {
                const T *halfptr = &half[(ix + 1) * incx + iy1 * incy];
                fulla0[iy1 * sincy + iz1] = a.a[0][ic] * halfptr[ic*incx + zstart + iz1];
            }
        }
        __syncthreads();
        for (int iy1 = 0; iy1 < dimy + 2*images; iy1++)
        {
            __syncthreads();
            for(int iz1 = blockIdx.y*blockDim.y + threadIdx.y;iz1 < zlen+2*images;iz1 += gridDim.y*blockDim.y)
            {
                const T *halfptr = &half[(ix + 1) * incx + iy1 * incy];
                fulla1[iy1 * sincy + iz1] = stencil(halfptr + zstart + iz1, incx);
            }
        }
        __syncthreads();

        for (int iy1 = 0; iy1 < dimy; iy1++)
        {
            T *full_tmp = &fulla0[(iy1 + 1) * sincy + iz1];
            fullb0[(2 * iy1 + 0) * sincy + iz1] = a.a[0][ic] * full_tmp[ic*sincy];
            fullb0[(2 * iy1 + 1) * sincy + iz1] = stencil(full_tmp, sincy);
        }
        __syncthreads();
        for (int iy1 = 0; iy1 < dimy; iy1++)
        {
            T *full_tmp = &fulla1[(iy1 + 1) * sincy + iz1];
            fullb1[(2 * iy1 + 0) * sincy + iz1] = a.a[0][ic] * full_tmp[ic*sincy];
            fullb1[(2 * iy1 + 1) * sincy + iz1] = stencil(full_tmp, sincy);
        }
        __syncthreads();
        // This last block writes back to global memory
        for (int iy1 = 0; iy1 < 2*dimy; iy1++)
        {
            for(int iz1 = blockIdx.y*blockDim.y + threadIdx.y;iz1 < zlen;iz1 += gridDim.y*blockDim.y)
            {
                T *full_tmp = &fullb0[iy1 * sincy + iz1 + 1];
                full[2*ix * incx2 + iy1 * incy2 + 2 * (zstart+iz1) + 0] = a.a[0][ic] * full_tmp[ic];
                full[2*ix * incx2 + iy1 * incy2 + 2 * (zstart+iz1) + 1] = stencil(full_tmp, 1);
            }
        }
        __syncthreads();
        for (int iy1 = 0; iy1 < 2*dimy; iy1++)
        {
            for(int iz1 = blockIdx.y*blockDim.y + threadIdx.y;iz1 < zlen;iz1 += gridDim.y*blockDim.y)
            {
                T *full_tmp = &fullb1[iy1 * sincy + iz1 + 1];
                full[(2*ix + 1) * incx2 + iy1 * incy2 + 2 * (zstart+iz1) + 0] = a.a[0][ic] * full_tmp[ic];
                full[(2*ix + 1) * incx2 + iy1 * incy2 + 2 * (zstart+iz1) + 1] = stencil(full_tmp, 1);
            }
        }
    }
}
#endif

std::vector<void *> abufs;
std::vector<void *> hbufs;
std::vector<float *> rbufs;

template <typename T>
void init_orthorhombic_gpu_prolong(int dimx, int dimy, int dimz)
{
    int order = MAX_PROLONG_ORDER;
    size_t rbufsize = 8*dimx*dimy*dimz*sizeof(double);
    size_t bufsize = (dimx + order)*(dimy + order)*(dimz + order)*sizeof(T);

    // Check if just clearing the accumulators
    if(rbufs.size() > 0)
    {
        //GpuFill(rbufs[getThreadId()], 8*dimx*dimy*dimz, 0.0);
        return;
    }

    int max_threads = getThreadNum();
    abufs.resize(max_threads);
    rbufs.resize(max_threads);
    hbufs.resize(max_threads);
    for(int i=0;i < max_threads;i++)
    {
        hipMalloc((void **)&abufs[i], bufsize);
        hipMallocHost((void **)&hbufs[i], rbufsize);
        hipMalloc((void **)&rbufs[i], rbufsize);
    }
    hipDeviceSynchronize();
}


template void init_orthorhombic_gpu_prolong<double>(int, int, int);
template void init_orthorhombic_gpu_prolong<std::complex<double>>(int, int, int);
template void init_orthorhombic_gpu_prolong<float>(int, int, int);
template void init_orthorhombic_gpu_prolong<std::complex<float>>(int, int, int);

template void Prolong::prolong_ortho_gpu<float,6>(float * , float *, int, int, int);
template void Prolong::prolong_ortho_gpu<std::complex<float>,6>(std::complex<float> * , std::complex<float> *, int, int, int);

template void Prolong::prolong_ortho_gpu<float,8>(float * , float *, int, int, int);
template void Prolong::prolong_ortho_gpu<std::complex<float>,8>(std::complex<float> * , std::complex<float> *, int, int, int);

template void Prolong::prolong_ortho_gpu<float,10>(float * , float *, int, int, int);
template void Prolong::prolong_ortho_gpu<std::complex<float>,10>(std::complex<float> * , std::complex<float> *, int, int, int);

template void Prolong::prolong_ortho_gpu<float,12>(float * , float *, int, int, int);
template void Prolong::prolong_ortho_gpu<std::complex<float>,12>(std::complex<float> * , std::complex<float> *, int, int, int);

template void prolong_ortho_gpu_internal<float,3>(float * , float *, int, int, int, float (&)[MAX_PROLONG_RATIO][MAX_PROLONG_ORDER]);
template void prolong_ortho_gpu_internal<std::complex<float>,3>(std::complex<float> * , std::complex<float> *, int, int, int, float (&)[MAX_PROLONG_RATIO][MAX_PROLONG_ORDER]);

template void prolong_ortho_gpu_internal<float,4>(float * , float *, int, int, int, float (&)[MAX_PROLONG_RATIO][MAX_PROLONG_ORDER]);
template void prolong_ortho_gpu_internal<std::complex<float>,4>(std::complex<float> * , std::complex<float> *, int, int, int, float (&)[MAX_PROLONG_RATIO][MAX_PROLONG_ORDER]);

template void prolong_ortho_gpu_internal<float,5>(float * , float *, int, int, int, float (&)[MAX_PROLONG_RATIO][MAX_PROLONG_ORDER]);
template void prolong_ortho_gpu_internal<std::complex<float>,5>(std::complex<float> * , std::complex<float> *, int, int, int, float (&)[MAX_PROLONG_RATIO][MAX_PROLONG_ORDER]);

template void prolong_ortho_gpu_internal<float,6>(float * , float *, int, int, int, float (&)[MAX_PROLONG_RATIO][MAX_PROLONG_ORDER]);
template void prolong_ortho_gpu_internal<std::complex<float>,6>(std::complex<float> * , std::complex<float> *, int, int, int, float (&)[MAX_PROLONG_RATIO][MAX_PROLONG_ORDER]);


#include "RmgTimer.h"

template <typename T, int ord>
void Prolong::prolong_ortho_gpu(T *full, 
                   T *half, 
                   const int dimx,
                   const int dimy,
                   const int dimz)
{
    if constexpr (std::is_same_v<T, float> && ord==6)
        prolong_ortho_gpu_internal<float, 3>(full, half, dimx, dimy, dimz, af);
    if constexpr (std::is_same_v<T, float> && ord==8)
        prolong_ortho_gpu_internal<float, 4>(full, half, dimx, dimy, dimz, af);
    if constexpr (std::is_same_v<T, float> && ord==10)
        prolong_ortho_gpu_internal<float, 5>(full, half, dimx, dimy, dimz, af);
    if constexpr (std::is_same_v<T, float> && ord==12)
        prolong_ortho_gpu_internal<float, 6>(full, half, dimx, dimy, dimz, af);

}

template <typename T, int images>
void prolong_ortho_gpu_internal(T *full, 
                   T *half, 
                   const int dimx,
                   const int dimy,
                   const int dimz,
                   float (&a)[MAX_PROLONG_RATIO][MAX_PROLONG_ORDER])
{

    pcoeff agpu;
    for(int i=0;i<MAX_PROLONG_RATIO;i++)
    {
        for(int j=0;j<MAX_PROLONG_ORDER;j++)
        {
            agpu.a[i][j] = a[i][j];
        }
    }

    dim3 Grid, Block, Grid1, Block1;
    hipStream_t stream = getGpuStream();
    std::vector<int> zstart, zlen, smem_sizes;
    int smem_limit = ct.smemSize[ct.hip_dev] - 4092;

    auto smem_needed = [&](const int dimy, int dimz) {
        int val = 2*(dimy + 2*images) * (dimz + 2*images) +
                   2*(2*dimy) * (dimz + 2*images);
        return val * sizeof(T);
    };

    int chunksize = dimz;
    // Get max possible chunksize
    while(smem_needed(dimy, chunksize) > smem_limit) chunksize--;
    int chunks = dimz / chunksize;
    if(dimz % chunksize) chunks++;
    for(int i=0;i < chunks;i++)
    {
        zstart.emplace_back(i*chunksize);
        zlen.emplace_back(chunksize);
        if((zstart[i] + zlen[i]) > dimz) zlen[i] = dimz - zstart[i];
        smem_sizes.emplace_back(smem_needed(dimy, chunksize));
        //printf("DDDD  %d  %d  %d  %d  %d\n",i, dimz, zstart[i], smem_sizes[i], zlen[i]);fflush(NULL);
    }

    int tid = getThreadId();
    int fbasis = 8*dimx*dimy*dimz;
    int sbasis = (dimx+2*images)*(dimy+2*images)*(dimz+2*images);

    Grid.x = dimx;
    Grid.y = 1;

    T *hptr = (T *)hbufs[tid];
    RmgTimer *RT = new RmgTimer("Prolong copy to");
    std::copy(half, half+sbasis, hptr);
    hipMemcpyAsync(abufs[tid], hbufs[tid], sbasis*sizeof(T), hipMemcpyHostToDevice, stream);
    hipStreamSynchronize(stream);
    delete RT;

    RT = new RmgTimer("Prolong kernel");

    for(int i=0;i < chunks;i++)
    {
        Block.x = 1;
        Block.y = (zlen[i] + 2*images);
        hipLaunchKernelGGL(HIP_KERNEL_NAME(prolong_ortho_kernel<T, images>), 
                   Grid,
                   Block,
                   smem_sizes[i],
                   stream,
                   (T *)rbufs[tid],
                   (T *)abufs[tid],
                   zstart[i],
                   zlen[i],
                   dimx, dimy, dimz, tid, agpu);
    }

    hipStreamSynchronize(stream);
    delete RT;

    RT = new RmgTimer("Prolong copy back");
    hipMemcpyAsync(hbufs[tid], rbufs[tid], fbasis*sizeof(T), hipMemcpyDeviceToHost, stream);
    hipStreamSynchronize(stream);
    std::copy(hptr, hptr+fbasis, full);
    delete RT;
}
