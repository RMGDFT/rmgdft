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

void DsygvjDriver(double *A, double *B, double *eigs, double *work, int worksize, int n, int ld)
{

    cusolverStatus_t cu_status;
    int lwork, *devInfo;
    const cusolverEigType_t itype = CUSOLVER_EIG_TYPE_1;
    const cusolverEigMode_t jobz = CUSOLVER_EIG_MODE_VECTOR; // compute eigenvectors.
    const cublasFillMode_t  uplo = CUBLAS_FILL_MODE_LOWER;
    syevjInfo_t dsygvj_params = NULL;

    cu_status = cusolverDnCreateSyevjInfo(&dsygvj_params);
    if(cu_status != CUSOLVER_STATUS_SUCCESS) rmg_error_handler (__FILE__, __LINE__, " cusolverDnCreateDsygvjInfo failed.");

    cu_status = cusolverDnDsygvj_bufferSize(ct.cusolver_handle, itype, jobz, uplo, n, A, n, B, n, eigs, &lwork, dsygvj_params);
    if(cu_status != CUSOLVER_STATUS_SUCCESS) rmg_error_handler (__FILE__, __LINE__, " cusolverDnDsyevd_bufferSize failed.");
    if(lwork > worksize) rmg_error_handler (__FILE__, __LINE__, " DsygvjDriver: provided workspace too small.");

    RmgGpuError(__FILE__, __LINE__, gpuMalloc((void **)&devInfo, sizeof(int) ), "Problem with gpuMalloc");
    double abstol = 1.0e-5;
    abstol = std::min(abstol, ct.scf_accuracy);
    cusolverDnXsyevjSetTolerance( dsygvj_params, abstol);

    cu_status = cusolverDnDsygvj(ct.cusolver_handle, itype, jobz, uplo, n, A, n, B, n, eigs, work, worksize, devInfo, dsygvj_params);
    int info;
    gpuMemcpy(&info, devInfo, sizeof(int), gpuMemcpyDeviceToHost);
    if(cu_status != CUSOLVER_STATUS_SUCCESS) rmg_error_handler (__FILE__, __LINE__, " cusolverDnDsygvj failed.");

    gpuFree(devInfo);
    if (dsygvj_params) cusolverDnDestroySyevjInfo(dsygvj_params);

}

#elif HIP_ENABLED

void DsygvjDriver(double *A, double *B, double *eigs, double *work, int worksize, int n, int ld)
{
    const rocblas_eform itype = rocblas_eform_ax;
    const rocblas_esort sortdir = rocblas_esort_ascending;
    const rocblas_evect jobz = rocblas_evect_original;
    const rocblas_fill uplo = rocblas_fill_lower;
    double abstol = 1.0e-5;
    abstol = std::min(abstol, ct.scf_accuracy);
    double *devResidual = NULL;
    rocblas_int max_sweeps = 15;
    int n_sweeps;
    int info;
    double residual;
    int *dev_n_sweeps;
    rocblas_int *devInfo;
    rocblas_status status;
    RmgGpuError(__FILE__, __LINE__, gpuMalloc((void **)&devInfo, sizeof(int) ), "Problem with gpuMalloc");
    RmgGpuError(__FILE__, __LINE__, gpuMalloc((void **)&devResidual, sizeof(double) ), "Problem with gpuMalloc");
    RmgGpuError(__FILE__, __LINE__, gpuMalloc((void **)&dev_n_sweeps, sizeof(int) ), "Problem with gpuMalloc");

    double tstart = my_crtc();
    status = rocsolver_dsygvj(ct.roc_handle,
                             itype,
                             jobz,
                             uplo,
                             n,
                             A,
                             ld,
                             B,
                             ld,
                             abstol,
                             devResidual,
                             max_sweeps,
                             dev_n_sweeps,
                             eigs,
                             devInfo);
    if(ct.verbose && pct.gridpe==0)
        printf("\nrocsolver_dsygvj time = %14.6f\n", my_crtc() - tstart);
    gpuMemcpy(&info, devInfo, sizeof(int), gpuMemcpyDeviceToHost);
    gpuMemcpy(&residual, devResidual, sizeof(double), gpuMemcpyDeviceToHost);
    gpuMemcpy(&n_sweeps, dev_n_sweeps, sizeof(int), gpuMemcpyDeviceToHost);
    gpuFree(dev_n_sweeps);
    gpuFree(devResidual);
    gpuFree(devInfo);
    if(status != 0) rmg_error_handler (__FILE__, __LINE__, " rocsolver_dsygvj failed.");
    if(ct.verbose && pct.gridpe==0)
        printf("rocsolver_dsygvj  %d  %d  %e  %d\n",status, n_sweeps, residual, info);
}


#else

void DsygvjDriver(double *A, double *B, double *eigs, double *work, int worksize, int n, int ld)
{
    // Redirect to Dsygvd since the Jacobi driver is not standard in CPU libraries
    DsygvdDriver(A, B, eigs, work, worksize, n, ld);

}
#endif
