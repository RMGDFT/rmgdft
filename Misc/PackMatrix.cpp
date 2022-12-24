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
#include "transition.h"


#if CUDA_ENABLED
    #include <cuda.h>
    #include <cuda_runtime_api.h>
    #include <cublas_v2.h>
#endif

// These are used to pack and unpack symmetric and Hermetian matrices from a full NxN
// format to a packed format. Lower triangular format is assumed for the NxN matrix.
void PackSqToTr(char *uplo, int N, double *Sq, int lda, double *Tr)
{
#if CUDA_ENABLED
    if(!ct.use_cublasxt)
    {
        cublasFillMode_t cu_uplo = CUBLAS_FILL_MODE_LOWER;
        if(!strcmp(uplo, "l")) cu_uplo = CUBLAS_FILL_MODE_LOWER;
        if(!strcmp(uplo, "L")) cu_uplo = CUBLAS_FILL_MODE_LOWER;
        if(!strcmp(uplo, "u")) cu_uplo = CUBLAS_FILL_MODE_UPPER;
        if(!strcmp(uplo, "U")) cu_uplo = CUBLAS_FILL_MODE_UPPER;
        cublasDtrttp ( ct.cublas_handle, cu_uplo, N, Sq, lda, Tr);
        DeviceSynchronize();
        return;
    }
#endif
    int info; 
    dtrttp(uplo, &N, Sq, &lda, Tr, &info);
}

void PackSqToTr(char *uplo, int N, double *Sq, int lda, float *Tr)
{
    float tbuf[4096];
    float *fSq = (float *)Sq;
    double *tSq = Sq;
    size_t blocks = ((size_t) N * (size_t)lda) / 4096;
    size_t rem = ((size_t) N * (size_t)lda) % 4096;
    for(size_t i=0;i < blocks;i++)
    {
        for(int j=0;j < 4096;j++) tbuf[j] = (float)tSq[j];
        std::copy(tbuf, tbuf + 4096,fSq);
        tSq += 4096;
        fSq += 4096;
    }
    for(size_t j=0;j < rem;j++) tbuf[j] = (float)tSq[j];
    for(size_t j=0;j < rem;j++) fSq[j] = tbuf[j];
    
    fSq = (float *)Sq;

#if CUDA_ENABLED
    if(!ct.use_cublasxt)
    {
        cublasFillMode_t cu_uplo = CUBLAS_FILL_MODE_LOWER;
        if(!strcmp(uplo, "l")) cu_uplo = CUBLAS_FILL_MODE_LOWER;
        if(!strcmp(uplo, "L")) cu_uplo = CUBLAS_FILL_MODE_LOWER;
        if(!strcmp(uplo, "u")) cu_uplo = CUBLAS_FILL_MODE_UPPER;
        if(!strcmp(uplo, "U")) cu_uplo = CUBLAS_FILL_MODE_UPPER;
        cublasStrttp ( ct.cublas_handle, cu_uplo, N, fSq, lda, Tr);
        DeviceSynchronize();
        return;
    }
#endif
    int info; 
    strttp(uplo, &N, fSq, &lda, Tr, &info);
}


void PackSqToTr(char *uplo, int N, std::complex<double> *Sq, int lda, std::complex<double> *Tr)
{
#if CUDA_ENABLED
    if(!ct.use_cublasxt)
    {
        cublasFillMode_t cu_uplo = CUBLAS_FILL_MODE_LOWER;
        if(!strcmp(uplo, "l")) cu_uplo = CUBLAS_FILL_MODE_LOWER;
        if(!strcmp(uplo, "L")) cu_uplo = CUBLAS_FILL_MODE_LOWER;
        if(!strcmp(uplo, "u")) cu_uplo = CUBLAS_FILL_MODE_UPPER;
        if(!strcmp(uplo, "U")) cu_uplo = CUBLAS_FILL_MODE_UPPER;
        cublasZtrttp ( ct.cublas_handle, cu_uplo, N, (cuDoubleComplex*)Sq, lda, (cuDoubleComplex*)Tr);
        DeviceSynchronize();
        return;
    }
#endif
    int info; 
    ztrttp(uplo, &N, Sq, &lda, Tr, &info);
}

void PackSqToTr(char *uplo, int N, std::complex<double> *Sq, int lda, std::complex<float> *Tr)
{
    float tbuf[1024];
    float *fSq = (float *)Sq;
    double *tSq = (double *)Sq;
    size_t blocks = (2*(size_t) N * (size_t)lda) / 1024;
    size_t rem = (2*(size_t) N * (size_t)lda) % 1024;
    for(size_t i=0;i < blocks;i++)
    {
        for(int j=0;j < 1024;j++) tbuf[j] = (float)tSq[j];
        for(int j=0;j < 1024;j++) fSq[j] = tbuf[j];
        tSq += 1024;
        fSq += 1024;
    }
    for(size_t j=0;j < rem;j++) tbuf[j] = (float)tSq[j];
    for(size_t j=0;j < rem;j++) fSq[j] = tbuf[j];
    
    fSq = (float *)Sq;
#if CUDA_ENABLED
    if(!ct.use_cublasxt)
    {
        cublasFillMode_t cu_uplo = CUBLAS_FILL_MODE_LOWER;
        if(!strcmp(uplo, "l")) cu_uplo = CUBLAS_FILL_MODE_LOWER;
        if(!strcmp(uplo, "L")) cu_uplo = CUBLAS_FILL_MODE_LOWER;
        if(!strcmp(uplo, "u")) cu_uplo = CUBLAS_FILL_MODE_UPPER;
        if(!strcmp(uplo, "U")) cu_uplo = CUBLAS_FILL_MODE_UPPER;
        cublasCtrttp ( ct.cublas_handle, cu_uplo, N, (cuComplex*)fSq, lda, (cuComplex*)Tr);
        DeviceSynchronize();
        return;
    }
#endif
    int info; 
    ctrttp(uplo, &N, (std::complex<float> *)fSq, &lda, Tr, &info);
}


void UnPackSqToTr(char *uplo, int N, double *Sq, int lda, double *Tr)
{
#if CUDA_ENABLED
    if(!ct.use_cublasxt)
    {
        cublasFillMode_t cu_uplo = CUBLAS_FILL_MODE_LOWER;
        if(!strcmp(uplo, "l")) cu_uplo = CUBLAS_FILL_MODE_LOWER;
        if(!strcmp(uplo, "L")) cu_uplo = CUBLAS_FILL_MODE_LOWER;
        if(!strcmp(uplo, "u")) cu_uplo = CUBLAS_FILL_MODE_UPPER;
        if(!strcmp(uplo, "U")) cu_uplo = CUBLAS_FILL_MODE_UPPER;
        cublasDtpttr ( ct.cublas_handle, cu_uplo, N, Tr, Sq, lda);
        DeviceSynchronize();
        return;
    }
#endif
    int info;
    dtpttr(uplo, &N, Tr, Sq, &lda, &info);
}

void UnPackSqToTr(char *uplo, int N, double *Sq, int lda, float *Tr)
{
    float *fSq = (float *)Sq;
#if CUDA_ENABLED
    if(!ct.use_cublasxt)
    {
        cublasFillMode_t cu_uplo = CUBLAS_FILL_MODE_LOWER;
        if(!strcmp(uplo, "l")) cu_uplo = CUBLAS_FILL_MODE_LOWER;
        if(!strcmp(uplo, "L")) cu_uplo = CUBLAS_FILL_MODE_LOWER;
        if(!strcmp(uplo, "u")) cu_uplo = CUBLAS_FILL_MODE_UPPER;
        if(!strcmp(uplo, "U")) cu_uplo = CUBLAS_FILL_MODE_UPPER;
        cublasStpttr ( ct.cublas_handle, cu_uplo, N, Tr, fSq, lda);
        DeviceSynchronize();
        return;
    }
#endif
    int info;
    stpttr(uplo, &N, Tr, fSq, &lda, &info);

    double *dTr = (double *)Tr;
    size_t stop = (size_t)N * (size_t)lda;
    for(size_t i=0;i < stop;i++) dTr[i] = (double)fSq[i];
    //for(size_t i=0;i < stop;i++) Sq[i] = dTr[i];
    std::copy(dTr, dTr + stop, Sq);
}

void UnPackSqToTr(char *uplo, int N, std::complex<double> *Sq, int lda, std::complex<double> *Tr)
{
#if CUDA_ENABLED
    if(!ct.use_cublasxt)
    {
        cublasFillMode_t cu_uplo = CUBLAS_FILL_MODE_LOWER;
        if(!strcmp(uplo, "l")) cu_uplo = CUBLAS_FILL_MODE_LOWER;
        if(!strcmp(uplo, "L")) cu_uplo = CUBLAS_FILL_MODE_LOWER;
        if(!strcmp(uplo, "u")) cu_uplo = CUBLAS_FILL_MODE_UPPER;
        if(!strcmp(uplo, "U")) cu_uplo = CUBLAS_FILL_MODE_UPPER;
        cublasZtpttr ( ct.cublas_handle, cu_uplo, N, (cuDoubleComplex*)Tr, (cuDoubleComplex*)Sq, lda);
        DeviceSynchronize();
    }
#endif
    int info;
    ztpttr(uplo, &N, Tr, Sq, &lda, &info);
}

void UnPackSqToTr(char *uplo, int N, std::complex<double> *Sq, int lda, std::complex<float> *Tr)
{
    float *fSq = (float *)Sq;
#if CUDA_ENABLED
    if(!ct.use_cublasxt)
    {
        cublasFillMode_t cu_uplo = CUBLAS_FILL_MODE_LOWER;
        if(!strcmp(uplo, "l")) cu_uplo = CUBLAS_FILL_MODE_LOWER;
        if(!strcmp(uplo, "L")) cu_uplo = CUBLAS_FILL_MODE_LOWER;
        if(!strcmp(uplo, "u")) cu_uplo = CUBLAS_FILL_MODE_UPPER;
        if(!strcmp(uplo, "U")) cu_uplo = CUBLAS_FILL_MODE_UPPER;
        cublasCtpttr ( ct.cublas_handle, cu_uplo, N, (cuComplex*)Tr, (cuComplex*)fSq, lda);
        DeviceSynchronize();
    }
#endif
    int info;
    ctpttr(uplo, &N, Tr, (std::complex<float> *)fSq, &lda, &info);

    double *dTr = (double *)Tr;
    double *dSq = (double *)Sq;
    size_t stop = 2 * (size_t)N * (size_t)lda;
    for(size_t i=0;i < stop;i++) dTr[i] = (double)fSq[i];
    for(size_t i=0;i < stop;i++) dSq[i] = dTr[i];
}
