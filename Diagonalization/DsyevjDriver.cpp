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

void DsyevjDriver(double *A, double *eigs, double *work, int worksize, int n, int ld)
{

    cusolverStatus_t cu_status;
    int lwork, *devInfo;
    const cusolverEigMode_t jobz = CUSOLVER_EIG_MODE_VECTOR; // compute eigenvectors.
    const cublasFillMode_t  uplo = CUBLAS_FILL_MODE_LOWER;
    syevjInfo_t dsyevj_params = NULL;

    cu_status = cusolverDnCreateSyevjInfo(&dsyevj_params);
    if(cu_status != CUSOLVER_STATUS_SUCCESS) rmg_error_handler (__FILE__, __LINE__, " cusolverDnCreateDsyevjInfo failed.");

    cu_status = cusolverDnDsyevj_bufferSize(ct.cusolver_handle, jobz, uplo, n, A, n, work, &lwork, dsyevj_params);
    if(cu_status != CUSOLVER_STATUS_SUCCESS) rmg_error_handler (__FILE__, __LINE__, " cusolverDnDsyevj_bufferSize failed.");
    if(lwork > worksize) rmg_error_handler (__FILE__, __LINE__, " DsyevjDriver: provided workspace too small.");

    RmgGpuError(__FILE__, __LINE__, gpuMalloc((void **)&devInfo, sizeof(int) ), "Problem with gpuMalloc");

    cu_status = cusolverDnDsyevj(ct.cusolver_handle, jobz, uplo, n, A, n, eigs, work, lwork, devInfo, dsyevj_params);
    int info;
    gpuMemcpy(&info, devInfo, sizeof(int), gpuMemcpyDeviceToHost);
    if(cu_status != CUSOLVER_STATUS_SUCCESS) rmg_error_handler (__FILE__, __LINE__, " cusolverDnDsyevj failed.");

    gpuFree(devInfo);
    if (dsyevj_params) cusolverDnDestroySyevjInfo(dsyevj_params);
    DeviceSynchronize();
}

#elif HIP_ENABLED

#include  <rocsolver/rocsolver.h>

void DsyevjDriver(double *A, double *eigs, double *work, int worksize, int n, int ld)
{

    const rocblas_esort sortdir = rocblas_esort_ascending;
    const rocblas_evect jobz = rocblas_evect_original;
    const rocblas_fill uplo = rocblas_fill_upper;
    const double abstol = 1.0e-9;
    double *dev_residual = NULL;
    double *residual;
    double *devResidual;
    const rocblas_int max_sweeps = 100;
    int n_sweeps, *dev_n_sweeps = NULL;
    rocblas_int *devInfo;
    int info;
    rocblas_status status;

    RmgGpuError(__FILE__, __LINE__, gpuMalloc((void **)&devInfo, sizeof(int) ), "Problem with gpuMalloc");
    RmgGpuError(__FILE__, __LINE__, gpuMalloc((void **)&devResidual, sizeof(double) ), "Problem with gpuMalloc");
    RmgGpuError(__FILE__, __LINE__, gpuMalloc((void **)&dev_n_sweeps, sizeof(int) ), "Problem with gpuMalloc");

    status = rocsolver_dsyevj(ct.roc_handle,
                              sortdir,
                              jobz,
                              uplo,
                              n,
                              A,
                              ld,
                              abstol,
                              devResidual,
                              max_sweeps,
                              dev_n_sweeps,
                              eigs,
                              devInfo);

    gpuMemcpy(&info, devInfo, sizeof(int), gpuMemcpyDeviceToHost);
    gpuMemcpy(&residual, devResidual, sizeof(double), gpuMemcpyDeviceToHost);
    gpuMemcpy(&n_sweeps, dev_n_sweeps, sizeof(int), gpuMemcpyDeviceToHost);
    gpuFree(dev_n_sweeps);
    gpuFree(devResidual);
    gpuFree(devInfo);
    if(status != 0) rmg_error_handler (__FILE__, __LINE__, " rocsolver_dsygvj failed.");
    //printf("RRRRR  %d  %d  %e  %d\n",status, n_sweeps, residual, info);


}

#else

void DsyevjDriver(double *A, double *eigs, double *work, int worksize, int n, int ld)
{
    // Redirect to Dsyevd since the Jacobi driver is not standard in CPU libraries
    DsyevdDriver(A, eigs, work, worksize, n, ld);

}

#endif
