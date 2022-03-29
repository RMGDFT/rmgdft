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


#include "GpuAlloc.h"
#include "rmg_error.h"
#include "transition.h"
#include "ErrorFuncs.h"
#include "Gpufuncs.h"
#include "RmgMatrix.h"
#include "blas.h"

#if CUDA_ENABLED
#include <cuda_runtime.h>
#include <cusolverDn.h>
#include <cusolverMg.h>

void DsygvdMgDriver(double *A, double *B, double *eigs, int n)
{
    RmgTimer RT("6-Diag-DsygvdMg");
    static int count = 0;
    cusolverMgHandle_t cusolverMg_handle = NULL;
    int current_device = 0;
    Cuda_error(cudaGetDevice(&current_device));

    int block_size = 128;
    int nGpus = ct.num_gpu_devices;
    int *device_list = ct.gpu_device_ids;
    const cusolverEigType_t itype = CUSOLVER_EIG_TYPE_1;
    const cusolverEigMode_t jobz = CUSOLVER_EIG_MODE_VECTOR; // compute eigenvectors.
    const cublasFillMode_t  uplo = CUBLAS_FILL_MODE_LOWER;

    cudaLibMgMatrixDesc_t desc_gpu;
    cudaLibMgGrid_t grid_gpu;

    cusolverMgGridMapping_t mapping = CUDALIBMG_GRID_MAPPING_COL_MAJOR;

    Cusolver_status(cusolverMgCreate(&cusolverMg_handle) );
    Cusolver_status(cusolverMgDeviceSelect(cusolverMg_handle, nGpus,
                ct.gpu_device_ids) );

    // enable peer access
    if(count == 0) {
        for(int igpu = 0; igpu < nGpus; igpu++)
        {
            Cuda_error(cudaSetDevice(igpu));
            for(int jgpu = 0; jgpu < nGpus; jgpu++)
            {
                if(igpu != jgpu)
                {
                    int canAccess = 0;
                    Cuda_error(cudaDeviceCanAccessPeer(&canAccess, igpu, jgpu));
                    if(canAccess)
                    {
                        Cuda_error(cudaDeviceEnablePeerAccess(jgpu, 0));
                    }
                }
            }
        }
        count++;
    }
    // setuo desc for nxn matrix 
    RmgTimer *RT1 = new RmgTimer("6-Diag-DsygvdMg: init");
    Cusolver_status(cusolverMgCreateDeviceGrid(&grid_gpu, 1, nGpus, device_list, mapping));
    Cusolver_status(cusolverMgCreateMatrixDesc(&desc_gpu, n, n, n, block_size, CUDA_R_64F, grid_gpu));

    //allocate mem for distributed matrix A B on devices.
    int num_blks = (n+block_size-1)/block_size;
    int max_nblks_per_gpu =(num_blks +nGpus -1)/nGpus;
    size_t size = n * block_size * max_nblks_per_gpu * sizeof(double);
    size_t len_lastblk = n - (num_blks-1) * block_size; 

    std::vector<double *> A_dist(nGpus, NULL);
    std::vector<double *> B_dist(nGpus, NULL);
    for(int igpu = 0; igpu < nGpus; igpu++)
    {
        Cuda_error(cudaSetDevice(device_list[igpu]));
        Cuda_error(cudaMalloc(&A_dist[igpu], size)); 
        Cuda_error(cudaMalloc(&B_dist[igpu], size)); 
        Cuda_error(cudaMemset(A_dist[igpu], 0, size));
        Cuda_error(cudaMemset(B_dist[igpu], 0, size));
    }
    Cuda_error(cudaDeviceSynchronize());

    delete RT1;

    // copy A and B from host to Device, the matrix is nxn on host and is distributed on devices

    RT1 = new RmgTimer("6-Diag-DsygvdMg: memcpy");
    for(int igpu = 0; igpu < nGpus; igpu++)
    {
        Cuda_error(cudaSetDevice(device_list[igpu]));
        for(int ib = igpu; ib < num_blks; ib+=nGpus)
        {

            size_t num_items = block_size * n * sizeof(double);
            if(ib == num_blks -1) num_items = len_lastblk * n * sizeof(double);
            double *h_A = &A[ib*block_size * n]; 
            double *d_A = A_dist[igpu] + (size_t)(ib/nGpus) *block_size * n; 
            double *h_B = &B[ib*block_size * n]; 
            double *d_B = &B_dist[igpu][(ib/nGpus) *block_size * n]; 
            Cuda_error(cudaMemcpy(d_A, h_A, num_items, cudaMemcpyHostToDevice));
            Cuda_error(cudaMemcpy(d_B, h_B, num_items, cudaMemcpyHostToDevice));
        }
    }
    Cuda_error(cudaDeviceSynchronize());
    delete RT1;

    int IA = 1, JA = 1;
    int64_t lwork_potrf, lwork_syevd, lwork;
    Cusolver_status(
            cusolverMgPotrf_bufferSize(cusolverMg_handle, uplo, n,
                reinterpret_cast<void **>(B_dist.data()), IA, JA,                                           
                desc_gpu, CUDA_R_64F, &lwork_potrf));
    Cusolver_status(cusolverMgSyevd_bufferSize(
                cusolverMg_handle, jobz, uplo,
                n, reinterpret_cast<void **>(A_dist.data()), IA, JA,
                desc_gpu, reinterpret_cast<void *>(eigs), CUDA_R_64F, CUDA_R_64F, &lwork_syevd));
    lwork_syevd = lwork_syevd;
    lwork = std::max(lwork_potrf, lwork_syevd);
    std::vector<double *> work_gpu(nGpus, NULL);

    if(pct.gridpe == 0) std::cout << "lwork in DsygvdMg" << lwork *sizeof(double)/1024.0/1024.0 << "MB" <<std::endl;
    for(int igpu = 0; igpu < nGpus; igpu++)
    {
        Cuda_error(cudaSetDevice(device_list[igpu]));
        int idevice;
        Cuda_error(cudaGetDevice(&idevice));
        Cuda_error(cudaMalloc(&work_gpu[igpu], lwork * sizeof(double))); 
        Cuda_error(cudaMemset(work_gpu[igpu], 0, lwork * sizeof(double))); 
    }

    Cuda_error(cudaDeviceSynchronize());
    RT1 = new RmgTimer("6-Diag-DsygvdMg: MgPotrf");
    int info;
    Cusolver_status(cusolverMgPotrf(
                cusolverMg_handle, uplo, n, reinterpret_cast<void **>(B_dist.data()), IA, JA,
                desc_gpu, CUDA_R_64F, reinterpret_cast<void **>(work_gpu.data()),
                lwork_potrf, &info /* host */
                ));
    Cuda_error(cudaDeviceSynchronize());
    delete RT1;

    RT1 = new RmgTimer("6-Diag-DsygvdMg: memcpy");
    // copy B=L * L^G from device to host, the matrix is nxn on host and is distributed on devices
    for(int igpu = 0; igpu < nGpus; igpu++)
    {
        for(int ib = igpu; ib < num_blks; ib+=nGpus)
        {

            size_t num_items = block_size * n * sizeof(double);
            if(ib == num_blks -1) num_items = len_lastblk * n * sizeof(double);
            double *h_B = &B[ib*block_size * n]; 
            double *d_B =  B_dist[igpu]+ib/nGpus *block_size * n; 
            Cuda_error(cudaMemcpy(h_B, d_B, num_items, cudaMemcpyDeviceToHost));
        }
    }

    Cuda_error(cudaDeviceSynchronize());
    delete RT1;
    //    dpotrf("L", &n, B, &n, &info);
    RT1 = new RmgTimer("6-Diag-DsygvdMg: dsygst");
    int itype_ds = 1;
    dsygst(&itype_ds, "L", &n, A, &n, B, &n, &info );
    delete RT1;

    //    dsyevd ("V", "L", &n, A, &n, eigs, work, &lwork_a, iwork, &liwork, &info);

    Cuda_error(cudaDeviceSynchronize());
    RT1 = new RmgTimer("6-Diag-DsygvdMg: memcpy");
    for(int igpu = 0; igpu < nGpus; igpu++)
    {
        for(int ib = igpu; ib < num_blks; ib+=nGpus)
        {

            size_t num_items = block_size * n * sizeof(double);
            if(ib == num_blks -1) num_items = len_lastblk * n * sizeof(double);
            double *h_A = &A[ib*block_size * n]; 
            double *d_A =  A_dist[igpu] + ib/nGpus *block_size * n; 
            Cuda_error(cudaMemcpy(d_A, h_A, num_items, cudaMemcpyHostToDevice));
        }
    }
    Cuda_error(cudaDeviceSynchronize());
    delete RT1;
    //Cuda_error(cudaMemcpy(A_dist[0], A, sizeof(double)*n*n, cudaMemcpyHostToDevice));

    RT1 = new RmgTimer("6-Diag-DsygvdMg: MgSyevd");
    Cuda_error(cudaSetDevice(current_device));
    Cusolver_status(cusolverMgSyevd(
                cusolverMg_handle, (cusolverEigMode_t)jobz, CUBLAS_FILL_MODE_LOWER,
                n, reinterpret_cast<void **>(A_dist.data()),  
                IA, JA, desc_gpu, reinterpret_cast<void **>(eigs), CUDA_R_64F, CUDA_R_64F, 
                reinterpret_cast<void **>(work_gpu.data()), lwork_syevd, &info));

    Cuda_error(cudaDeviceSynchronize());
    delete RT1;
    RT1 = new RmgTimer("6-Diag-DsygvdMg: memcpy");
    for(int igpu = 0; igpu < nGpus; igpu++)
    {
        for(int ib = igpu; ib < num_blks; ib+=nGpus)
        {

            size_t num_items = block_size * n * sizeof(double);
            if(ib == num_blks -1) num_items = len_lastblk * n * sizeof(double);
            double *h_A = &A[ib*block_size * n]; 
            double *d_A = &A_dist[igpu][ib/nGpus *block_size * n]; 
            Cuda_error(cudaMemcpy(h_A, d_A, num_items, cudaMemcpyDeviceToHost));
        }
    }
    Cuda_error(cudaDeviceSynchronize());
    delete RT1;

    RT1 = new RmgTimer("6-Diag-DsygvdMg: dtrsm");
    double one = 1.0;
    dtrsm( "Left", "L", "T", "N", &n, &n, &one, B, &n , A, &n );
    delete RT1;
    for(int igpu = 0; igpu < nGpus; igpu++)
    {
        Cuda_error(cudaSetDevice(device_list[igpu]));
        Cuda_error(cudaFree(work_gpu[igpu]));
        Cuda_error(cudaFree(B_dist[igpu]));
        Cuda_error(cudaFree(A_dist[igpu]));
    }
    Cuda_error(cudaSetDevice(current_device));

}

#else

void DsygvdMgDriver(double *A, double *B, double *eigs, int n)
{
    rmg_error_handler (__FILE__, __LINE__, " cusolverDsygvdMg not programmed without CUDA.");
}
#endif
