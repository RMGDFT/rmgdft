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

void DgetrftrsDriver(int n, int m, double *A, double *B)
{

    int info = 0;

    int *d_Ipiv = nullptr; /* pivoting sequence */
    int *d_info = nullptr; /* error info */

    int lwork = 0;            /* size of workspace */
    double *d_work = nullptr; /* device workspace for getrf */

    Cuda_error(cudaMalloc(reinterpret_cast<void **>(&d_Ipiv), sizeof(int) * n));
    Cuda_error(cudaMalloc(reinterpret_cast<void **>(&d_info), sizeof(int)));

    /* query working space of getrf */
    Cusolver_status(cusolverDnDgetrf_bufferSize(ct.cusolver_handle, n, n, A, n, &lwork));

    Cuda_error(cudaMalloc(reinterpret_cast<void **>(&d_work), sizeof(double) * lwork));

    /*  LU factorization */
    Cusolver_status(cusolverDnDgetrf(ct.cusolver_handle, n, n, A, n, d_work, d_Ipiv, d_info));


    /*
     *  solve A*X = B
     */
    Cusolver_status(cusolverDnDgetrs(ct.cusolver_handle, CUBLAS_OP_N, n, m, /* nrhs */
                                        A, n, d_Ipiv, B, n, d_info));

    /* free resources */
    Cuda_error(cudaFree(d_Ipiv));
    Cuda_error(cudaFree(d_info));
    Cuda_error(cudaFree(d_work));

}

#elif HIP_ENABLED
#if MAGMA_LIBS
#include <magma_v2.h>

void DgetrftrsDriver(int n, int m, double *A, double *B)
{
    rmg_error_handler (__FILE__, __LINE__, " dgestrs  not programmed.");
}
#else
#include <rocsolver.h>

void DgetrftrsDriver(int n, int m, double *A, double *B)
{
    rmg_error_handler (__FILE__, __LINE__, " dgetrs not programmed.");
}
#endif

#else

void DgetrftrsDriver(int n, int m, double *A, double *B)
{
    rmg_error_handler (__FILE__, __LINE__, " dgetrs not programmed.");

}
#endif
