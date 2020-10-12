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

void ZhegvdDriver(std::complex<double> *A, std::complex<double> *B, double *eigs, double *work, int worksize, int n, int ld)
{

    cusolverStatus_t cu_status;
    int lwork, *devInfo;
    const cusolverEigType_t itype = CUSOLVER_EIG_TYPE_1;
    const cusolverEigMode_t jobz = CUSOLVER_EIG_MODE_VECTOR; // compute eigenvectors.
    const cublasFillMode_t  uplo = CUBLAS_FILL_MODE_LOWER;
    cuDoubleComplex *zwork;

    cu_status = cusolverDnZhegvd_bufferSize(ct.cusolver_handle, itype, jobz, uplo, n, (cuDoubleComplex *)A, n, (cuDoubleComplex *)B, n, eigs, &lwork);
    if(cu_status != CUSOLVER_STATUS_SUCCESS) rmg_error_handler (__FILE__, __LINE__, " cusolverDnZheevd_bufferSize failed.");

    if(work == NULL)
    {
        zwork = (cuDoubleComplex *)GpuMallocManaged(lwork * sizeof(std::complex<double>));
    }
    else
    {
        if(lwork > worksize) rmg_error_handler (__FILE__, __LINE__, " ZhegvdDriver: provided workspace too small.");
        zwork = (cuDoubleComplex *)work;
    }

    RmgCudaError(__FILE__, __LINE__, cudaMalloc((void **)&devInfo, sizeof(int) ), "Problem with cudaMalloc");

    cu_status = cusolverDnZhegvd(ct.cusolver_handle, itype, jobz, uplo, n, (cuDoubleComplex *)A, n, (cuDoubleComplex *)B, n, eigs, (cuDoubleComplex *)zwork, lwork, devInfo);
    int info;
    cudaMemcpy(&info, devInfo, sizeof(int), cudaMemcpyDeviceToHost);
    if(cu_status != CUSOLVER_STATUS_SUCCESS) rmg_error_handler (__FILE__, __LINE__, " cusolverDnZhegvd failed.");

    cudaFree(devInfo);
    if(work == NULL) GpuFreeManaged(zwork);
}

#else

void ZhegvdDriver(std::complex<double> *A, std::complex<double> *B, double *eigs, double *work, int worksize, int n, int ld)
{
    char *cuplo = "l", *jobz="V";
    int lwork, info=0, *iwork, liwork, itype=1;

    liwork = 6*n;
    iwork = new int[liwork];
    lwork = worksize;
    int lrwork = 2*n*n + 6*n; 
    double *rwork = new double[lrwork];

    zhegvd(&itype, jobz, cuplo, &n, (double *)A, &n, (double *)B, &n, eigs, work, &lwork, rwork, &lrwork, iwork, &liwork, &info);

    if(info)
        rmg_error_handler (__FILE__, __LINE__, " zhegvd failed.");

    delete [] rwork;
    delete [] iwork;
}
#endif
