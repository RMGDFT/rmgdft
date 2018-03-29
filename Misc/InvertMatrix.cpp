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

template void InvertMatrix<double>(double *, double *, int);
template void InvertMatrix<std::complex<double>>(std::complex<double> *, std::complex<double> *, int);

// Inverts square matrix A and returns result in B. For the GPU case A and B must be
// located in either device or managed memory.

#if GPU_ENABLED

template <typename DataType> void InvertMatrix(DataType *A, DataType *B, int n)
{

    cusolverStatus_t cu_status;
    int Lwork;
    int *devIpiv, *devInfo;
    cublasOperation_t trans = CUBLAS_OP_N;
    DataType *Workspace;

    RmgCudaError(__FILE__, __LINE__, cudaMalloc((void **)&devIpiv, sizeof(int) *n), "Problem with cudaMalloc");
    RmgCudaError(__FILE__, __LINE__, cudaMalloc((void **)&devInfo, sizeof(int) ), "Problem with cudaMalloc");

    if(typeid(DataType) == typeid(double))
    {
        GpuFill((double *)B, n*n, 0.0);
        cu_status = cusolverDnDgetrf_bufferSize(ct.cusolver_handle, n, n, (double *)A, n, &Lwork);
        if(cu_status != CUSOLVER_STATUS_SUCCESS) rmg_error_handler (__FILE__, __LINE__, " cusolverDnDgetrf_bufferSize failed.");
        RmgCudaError(__FILE__, __LINE__, cudaMalloc((void **) &Workspace, sizeof(double) * std::max(Lwork, n)), "Problem with cudaMalloc");

        // Create unitary matrix
        GpuFill((double *)Workspace, n, 1.0);
        cublasDcopy(ct.cublas_handle, n, (double *)Workspace, 1, (double *)B, n+1);

        cu_status = cusolverDnDgetrf(ct.cusolver_handle, n, n, (double *)A, n, (double *)Workspace, devIpiv, devInfo );
        if(cu_status != CUSOLVER_STATUS_SUCCESS) rmg_error_handler (__FILE__, __LINE__, " cusolverDnDgetrf failed.");

        cu_status = cusolverDnDgetrs(ct.cusolver_handle, trans, n, n, (double *)A, n, devIpiv, (double *)B, n, devInfo );
        if(cu_status != CUSOLVER_STATUS_SUCCESS) rmg_error_handler (__FILE__, __LINE__, " cusolverDnDgetrs failed.");
    }
    else if(typeid(DataType) == typeid(std::complex<double>))
    {
        GpuFill((double *)B, 2*n*n, 0.0);
        cu_status = cusolverDnZgetrf_bufferSize(ct.cusolver_handle, n, n, (cuDoubleComplex *)A, n, &Lwork);
        if(cu_status != CUSOLVER_STATUS_SUCCESS) rmg_error_handler (__FILE__, __LINE__, " cusolverZnDgetrf_bufferSize failed.");
        RmgCudaError(__FILE__, __LINE__, cudaMalloc((void **) &Workspace, 2*sizeof(double) * std::max(Lwork, n)), "Problem with cudaMalloc");

        // Create unitary matrix
        GpuFill((double *)Workspace, n, 1.0);
        cublasDcopy(ct.cublas_handle, n, (double *)Workspace, 1, (double *)B, 2*(n+1));

        cu_status = cusolverDnZgetrf(ct.cusolver_handle, n, n, (cuDoubleComplex *)A, n, (cuDoubleComplex *)Workspace, devIpiv, devInfo );
        if(cu_status != CUSOLVER_STATUS_SUCCESS) rmg_error_handler (__FILE__, __LINE__, " cusolverZnDgetrf failed.");

        cu_status = cusolverDnZgetrs(ct.cusolver_handle, trans, n, n, (cuDoubleComplex *)A, n, devIpiv, (cuDoubleComplex *)B, n, devInfo );
        if(cu_status != CUSOLVER_STATUS_SUCCESS) rmg_error_handler (__FILE__, __LINE__, " cusolverDnDgetrs failed.");
    }

    cudaFree(Workspace);
    cudaFree(devInfo);
    cudaFree(devIpiv);

}

#else

template <typename DataType> void InvertMatrix(DataType *A, DataType *B, int n)
{

    DataType ONE_t(1.0);
    DataType ZERO_t(0.0);
    int *ipiv = new int[2*n]();
    int info = 0;

    // Create unitary matrix
    for (int idx = 0; idx < n*n; idx++) B[idx] = ZERO_t;
    for (int idx = 0; idx < n; idx++) B[idx * n + idx] = ONE_t;

    if(typeid(DataType) == typeid(double))
    {
        dgesv (&n, &n, (double *)A, &n, ipiv, (double *)B, &n, &info);
    }
    else if(typeid(DataType) == typeid(std::complex<double>))
    {
        zgesv (&n, &n, (double *)A, &n, ipiv, (double *)B, &n, &info);
    }
    if (info) {
        rmg_printf ("\n PE %d: p{d,z}gesv failed, info is %d", pct.gridpe, info);
        rmg_error_handler (__FILE__, __LINE__, " p{d,z}gesv failed");
    }

    delete [] ipiv;
}

#endif
