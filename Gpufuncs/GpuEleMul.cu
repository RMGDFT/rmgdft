/*
 *
 * Copyright 2018 The RMG Project Developers. See the COPYRIGHT file 
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



#if CUDA_ENABLED
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/transform.h>
#include <thrust/copy.h>
#include <thrust/fill.h>
#include <thrust/sequence.h>
#include <thrust/complex.h>
#include <thrust/complex.h>
#include <cuComplex.h>


__global__ void MulVec(double *dx, cuDoubleComplex *dy, int n)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    for (int i = idx; i < n; i += gridDim.x * blockDim.x) dy[i] = make_cuDoubleComplex(cuCreal(dy[i]) * dx[i], cuCimag(dy[i]) * dx[i]);
}

__global__ void MulVec1(double *dx, cuFloatComplex *dy, int n)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    for (int i = idx; i < n; i += gridDim.x * blockDim.x) dy[i] = make_cuFloatComplex(cuCrealf(dy[i]) * dx[i], cuCimagf(dy[i]) * dx[i]);
}

void GpuEleMul(double *dx, std::complex<double> *dy, int n, cudaStream_t stream)
{
    cudaStreamSynchronize(stream);
    int blockSize = 128;
    int numBlocks = (n + blockSize - 1) / n;
    MulVec<<<numBlocks, blockSize, 0, stream>>>(dx, (cuDoubleComplex *)dy, n);
    cudaStreamSynchronize(stream);
}

void GpuEleMul(double *dx, std::complex<float> *dy, int n, cudaStream_t stream)
{
    cudaStreamSynchronize(stream);
    int blockSize = 128;
    int numBlocks = (n + blockSize - 1) / n;
    MulVec1<<<numBlocks, blockSize, 0, stream>>>(dx, (cuFloatComplex *)dy, n);
    cudaStreamSynchronize(stream);
}
  

#if 0
// Fills a device pointer with a floating point value
void GpuFill(double *dptr, int n, double fillval)
{
    cudaDeviceSynchronize();
    int blockSize = 256;
    int numBlocks = (n + blockSize - 1) / n;
    FillMatrix<<<numBlocks, blockSize>>>(dptr, n, fillval);
    cudaDeviceSynchronize();
}

struct ElementWiseOp : public thrust::binary_function<double, thrust::complex<double> , thrust::complex<double> >
{
    __device__
    thrust::complex<double> operator()(double & v1, thrust::complex<double> & v2)
    {
        thrust::complex<double> res;
        res = v1 * v2;
        return res;
    }
};

void GpuEleMul(double *dx, std::complex<double> *dy, int n, cudaStream_t stream)
{
    //thrust::device_ptr<double> dx_dev = thrust::device_pointer_cast(dx);
    thrust::device_ptr<double> dx_dev(dx);
    //thrust::device_ptr<thrust::complex<double>> dy_dev = thrust::device_pointer_cast((thrust::complex<double> *)dy);
    thrust::device_ptr<thrust::complex<double>> dy_dev((thrust::complex<double> *)dy);

//    thrust::complex<double> * dz = (thrust::complex<double> *)dy;
    thrust::transform(thrust::cuda::par.on(stream), dx_dev, dx_dev + n, dy_dev, dy_dev, ElementWiseOp());
}

#endif
#endif
