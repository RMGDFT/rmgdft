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

    cu_status = cusolverDnDsygvj(ct.cusolver_handle, itype, jobz, uplo, n, A, n, B, n, eigs, work, worksize, devInfo, dsygvj_params);
    int info;
    gpuMemcpy(&info, devInfo, sizeof(int), gpuMemcpyDeviceToHost);
    if(cu_status != CUSOLVER_STATUS_SUCCESS) rmg_error_handler (__FILE__, __LINE__, " cusolverDnDsygvj failed.");

    gpuFree(devInfo);
    if (dsygvj_params) cusolverDnDestroySyevjInfo(dsygvj_params);

}

#else

void DsygvjDriver(double *A, double *B, double *eigs, double *work, int worksize, int n, int ld)
{
    // Redirect to Dsygvd since the Jacobi driver is not standard in CPU libraries
    DsygvdDriver(A, B, eigs, work, worksize, n, ld);

}
#endif
