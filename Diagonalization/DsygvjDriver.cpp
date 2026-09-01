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

#include "Gpufuncs.h"
#include "RmgMatrix.h"
#include "blas.h"

#if CUDA_ENABLED
#include <cuda_runtime.h>
#include <cusolverDn.h>

void DsygvjDriver(double *A, double *B, double *eigs, double *work, int worksize, int n, int ld)
{

    int lwork, *devInfo;
    const cusolverEigType_t itype = CUSOLVER_EIG_TYPE_1;
    const cusolverEigMode_t jobz = CUSOLVER_EIG_MODE_VECTOR; // compute eigenvectors.
    const cublasFillMode_t  uplo = CUBLAS_FILL_MODE_LOWER;
    syevjInfo_t dsygvj_params = NULL;

    rmg::error(cusolverDnCreateSyevjInfo(&dsygvj_params));

    rmg::error(cusolverDnDsygvj_bufferSize(ct.cusolver_handle, itype, jobz, uplo, n, A, n, B, n, eigs, &lwork, dsygvj_params));
    if(lwork > worksize) rmg::error(" DsygvjDriver: provided workspace too small.");

    gpuMalloc((void **)&devInfo, sizeof(int));
    double abstol = 1.0e-5;
    abstol = std::min(abstol, ct.scf_accuracy)/100.0;
    rmg::error(cusolverDnXsyevjSetTolerance( dsygvj_params, abstol));

    rmg::error(cusolverDnDsygvj(ct.cusolver_handle, itype, jobz, uplo, n, A, n, B, n, eigs, work, worksize, devInfo, dsygvj_params));
    int info;
    gpuMemcpy(&info, devInfo, sizeof(int), gpuMemcpyDeviceToHost);

    gpuFree(devInfo);
    if (dsygvj_params) rmg::error(cusolverDnDestroySyevjInfo(dsygvj_params));

}

#elif HIP_ENABLED

void DsygvjDriver(double *A, double *B, double *eigs, double *work, int worksize, int n, int ld)
{
    const rocblas_eform itype = rocblas_eform_ax;
    const rocblas_evect jobz = rocblas_evect_original;
    const rocblas_fill uplo = rocblas_fill_lower;
    double abstol = 1.0e-20;
    double *devResidual = NULL;
    rocblas_int max_sweeps = 50;
    int n_sweeps;
    int info;
    double residual;
    int *dev_n_sweeps;
    rocblas_int *devInfo;
    rocblas_status status;
    gpuMalloc((void **)&devInfo, sizeof(int));
    gpuMalloc((void **)&devResidual, sizeof(double));
    gpuMalloc((void **)&dev_n_sweeps, sizeof(int));

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
    if(status != 0) rmg::error(" rocsolver_dsygvj failed.");
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
