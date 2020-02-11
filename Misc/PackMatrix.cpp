/*
 *
 * Copyright 2020 The RMG Project Developers. See the COPYRIGHT file 
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

#include <complex>
#include "blas.h"
#include "main.h"

#if GPU_ENABLED
    #include <cuda.h>
    #include <cuda_runtime_api.h>
    #include <cublas_v2.h>
#endif

// These are used to pack and unpack symmetric and Hermetian matrices from a full NxN
// format to a packed format. Lower triangular format is assumed for the NxN matrix.
void PackSqToTr(char *uplo, int N, double *Sq, int lda, double *Tr)
{
#if GPU_ENABLED
    cublasFillMode_t cu_uplo = CUBLAS_FILL_MODE_LOWER;
    if(!strcmp(uplo, "l")) cu_uplo = CUBLAS_FILL_MODE_LOWER;
    if(!strcmp(uplo, "L")) cu_uplo = CUBLAS_FILL_MODE_LOWER;
    if(!strcmp(uplo, "u")) cu_uplo = CUBLAS_FILL_MODE_UPPER;
    if(!strcmp(uplo, "U")) cu_uplo = CUBLAS_FILL_MODE_UPPER;
    cublasDtrttp ( ct.cublas_handle, cu_uplo, N, Sq, lda, Tr);
    cudaDeviceSynchronize();
#else
    int info; 
    dtrttp(uplo, &N, Sq, &lda, Tr, &info);
#endif

}

void PackSqToTr(char *uplo, int N, std::complex<double> *Sq, int lda, std::complex<double> *Tr)
{
#if GPU_ENABLED
    cublasFillMode_t cu_uplo = CUBLAS_FILL_MODE_LOWER;
    if(!strcmp(uplo, "l")) cu_uplo = CUBLAS_FILL_MODE_LOWER;
    if(!strcmp(uplo, "L")) cu_uplo = CUBLAS_FILL_MODE_LOWER;
    if(!strcmp(uplo, "u")) cu_uplo = CUBLAS_FILL_MODE_UPPER;
    if(!strcmp(uplo, "U")) cu_uplo = CUBLAS_FILL_MODE_UPPER;
    cublasZtrttp ( ct.cublas_handle, cu_uplo, N, (cuDoubleComplex*)Sq, lda, (cuDoubleComplex*)Tr);
    cudaDeviceSynchronize();
#else
    int info; 
    ztrttp(uplo, &N, Sq, &lda, Tr, &info);
#endif
}

void UnPackSqToTr(char *uplo, int N, double *Sq, int lda, double *Tr)
{
#if GPU_ENABLED
    cublasFillMode_t cu_uplo = CUBLAS_FILL_MODE_LOWER;
    if(!strcmp(uplo, "l")) cu_uplo = CUBLAS_FILL_MODE_LOWER;
    if(!strcmp(uplo, "L")) cu_uplo = CUBLAS_FILL_MODE_LOWER;
    if(!strcmp(uplo, "u")) cu_uplo = CUBLAS_FILL_MODE_UPPER;
    if(!strcmp(uplo, "U")) cu_uplo = CUBLAS_FILL_MODE_UPPER;
    cublasDtpttr ( ct.cublas_handle, cu_uplo, N, Tr, Sq, lda);
    cudaDeviceSynchronize();
#else
    int info;
    dtpttr(uplo, &N, Tr, Sq, &lda, &info);
#endif
}

void UnPackSqToTr(char *uplo, int N, std::complex<double> *Sq, int lda, std::complex<double> *Tr)
{
#if GPU_ENABLED
    cublasFillMode_t cu_uplo = CUBLAS_FILL_MODE_LOWER;
    if(!strcmp(uplo, "l")) cu_uplo = CUBLAS_FILL_MODE_LOWER;
    if(!strcmp(uplo, "L")) cu_uplo = CUBLAS_FILL_MODE_LOWER;
    if(!strcmp(uplo, "u")) cu_uplo = CUBLAS_FILL_MODE_UPPER;
    if(!strcmp(uplo, "U")) cu_uplo = CUBLAS_FILL_MODE_UPPER;
    cublasZtpttr ( ct.cublas_handle, cu_uplo, N, (cuDoubleComplex*)Tr, (cuDoubleComplex*)Sq, lda);
    cudaDeviceSynchronize();
#else
    int info;
    ztpttr(uplo, &N, Tr, Sq, &lda, &info);
#endif
}
