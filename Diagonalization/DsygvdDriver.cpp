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

void DsygvdDriver_lapack(double *A, double *B, double *eigs, double *work, int worksize, int n, int ld);

#if CUDA_ENABLED
#include <cuda_runtime.h>
#include <cusolverDn.h>

void DsygvdDriver(double *A, double *B, double *eigs, double *work, int worksize, int n, int ld)
{

    cusolverStatus_t cu_status;
    int lwork, *devInfo;
    const cusolverEigType_t itype = CUSOLVER_EIG_TYPE_1;
    const cusolverEigMode_t jobz = CUSOLVER_EIG_MODE_VECTOR; // compute eigenvectors.
    const cublasFillMode_t  uplo = CUBLAS_FILL_MODE_LOWER;


    cu_status = cusolverDnDsygvd_bufferSize(ct.cusolver_handle, itype, jobz, uplo, n, A, n, B, n, eigs, &lwork);
    if(cu_status != CUSOLVER_STATUS_SUCCESS) rmg_error_handler (__FILE__, __LINE__, " cusolverDnDsyevd_bufferSize failed.");
    if(lwork > worksize) 
    {
        cudaFree(work);
        Cuda_error(cudaMalloc((void **)&work, lwork * sizeof(double)));
    }
    RmgGpuError(__FILE__, __LINE__, gpuMalloc((void **)&devInfo, sizeof(int) ), "Problem with gpuMalloc");

    cu_status = cusolverDnDsygvd(ct.cusolver_handle, itype, jobz, uplo, n, A, n, B, n, eigs, work, lwork, devInfo);
    int info;
    gpuMemcpy(&info, devInfo, sizeof(int), gpuMemcpyDeviceToHost);
    if(cu_status != CUSOLVER_STATUS_SUCCESS) rmg_error_handler (__FILE__, __LINE__, " cusolverDnDsygvd failed.");

    gpuFree(devInfo);
}

#elif HIP_ENABLED
#if MAGMA_LIBS
#include <magma_v2.h>

void DsygvdDriver(double *A, double *B, double *eigs, double *work, int worksize, int n, int ld)
{
    int itype = 1, info;
    int liwork = 6*n;
    int *iwork = new int[liwork];

    magma_dsygvd(itype, MagmaVec, MagmaLower, n, (double *)A, ld, (double *)B, ld, eigs, work, worksize, iwork, liwork, &info);

//    double vl = 0.0, vu = 0.0;
//    int eigs_found;
//    magma_dsygvdx_2stage(itype, MagmaVec, MagmaRangeAll, MagmaLower, n,
//             (double *)A, ld, (double *)B, ld, vl, vu, 1, n, &eigs_found,
//    eigs, work, worksize, iwork, liwork, &info);

    delete [] iwork;
}
#else
#include <rocsolver.h>

void DsygvdDriver(double *A, double *B, double *eigs, double *work, int worksize, int n, int ld)
{
    if(ct.subdiag_driver == SUBDIAG_LAPACK)
    {
        DsygvdDriver_lapack(A, B, eigs, work, worksize, n, ld);
        return;
    }

    rocblas_status status;
    rocblas_int *devInfo;
    const rocblas_evect jobz = rocblas_evect_original; // compute eigenvectors.
    const rocblas_fill uplo = rocblas_fill_lower;
    const rocblas_eform itype = rocblas_eform_ax;

    gpuSetDevice(ct.hip_dev);
    RmgGpuError(__FILE__, __LINE__, gpuMalloc((void **)&devInfo, sizeof(int) ), "Problem with gpuMalloc");
    status = rocsolver_dsygvd(ct.roc_handle, itype, jobz, uplo, n, A, ld,
                             B, ld, eigs, work, devInfo);
    int info;
    gpuMemcpy(&info, devInfo, sizeof(int), gpuMemcpyDeviceToHost);
    if(status != rocblas_status_success) rmg_error_handler (__FILE__, __LINE__, " rocsolver_dsygv failed.");

    gpuFree(devInfo);
}
#endif

#else

void DsygvdDriver(double *A, double *B, double *eigs, double *work, int worksize, int n, int ld)
{
    DsygvdDriver_lapack(A, B, eigs, work, worksize, n, ld);
}
#endif

void DsygvdDriver_lapack(double *A, double *B, double *eigs, double *work, int worksize, int n, int ld)
{
    char *cuplo = "l", *jobz="V";
    int lwork, info=0, *iwork, liwork, ione=1;

    liwork = 6*n;
    iwork = new int[liwork];

    lwork = worksize;

    dsygvd(&ione, jobz, cuplo, &n, A, &n, B, &n, eigs, work, &lwork, iwork, &liwork, &info);

    if(info)
        rmg_error_handler (__FILE__, __LINE__, " dsyevd failed.");

    delete [] iwork;
}
