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
#include "RmgMatrix.h"


#if CUDA_ENABLED
    #include <cuda.h>
    #include <cuda_runtime_api.h>
    #include <cublas_v2.h>
#endif


// Pack and convert versions
void dtrttp(char *uplo, int n_in, double *a, float *ap)
{
    size_t n = (size_t)n_in;
    if(!strcmp(uplo, "l") || !strcmp(uplo, "L"))
    {
        size_t k = 0;
        for(size_t j=0;j < n;j++)
        {
           for(size_t i=j;i < n;i++)
           {
              ap[k] = (float)a[j*n + i];
              k++;
           }
        }
    }
    else
    {
        size_t k = 0;
        for(size_t j=0;j < n;j++)
        {
            for(size_t i=0;i < j;i++)
            {
                ap[k] = (float)a[j*n + i];
                k++;
            }
        }
    }
}

void ztrttp(char *uplo, int n_in, std::complex<double> *a, std::complex<float> *ap)
{
    size_t n = (size_t)n_in;
    if(!strcmp(uplo, "l") || !strcmp(uplo, "L"))
    {
        size_t k = 0;
        for(size_t j=0;j < n;j++)
        {
            for(size_t i=j;i < n;i++)
            {
                ap[k] = a[j*n + i];
                k++;
            }
        }
    }
    else
    {
        size_t k = 0;
        for(size_t j=0;j < n;j++)
        {
            for(size_t i=0;i < j;i++)
            {
                ap[k] = a[j*n + i];
                k++;
            }
        }
    }
}

void dtpttr(char *uplo, int n_in, float *ap, double *a)
{
    size_t n = (size_t)n_in;

    if(!strcmp(uplo, "l") || !strcmp(uplo, "L"))
    {
        size_t k = 0;
        for(size_t j=0;j < n;j++)
        {
           for(size_t i=j;i < n;i++)
           {
              a[j*n + i] = (double)ap[k];
              k++;
           }
        }
    }
    else
    {
        size_t k = 0;
        for(size_t j=0;j < n;j++)
        {
            for(size_t i=0;i < j;i++)
            {
                a[j*n + i] = (double)ap[k];
                k++;
            }
        }
    }
}

void ztpttr(char *uplo, int n_in, std::complex<float> *ap, std::complex<double> *a)
{
    size_t n = (size_t)n_in;

    if(!strcmp(uplo, "l") || !strcmp(uplo, "L"))
    {
        size_t k = 0;
        for(size_t j=0;j < n;j++)
        {
           for(size_t i=j;i < n;i++)
           {
              a[j*n + i] = ap[k];
              k++;
           }
        }
    }
    else
    {
        size_t k = 0;
        for(size_t j=0;j < n;j++)
        {
            for(size_t i=0;i < j;i++)
            {
                a[j*n + i] = ap[k];
                k++;
            }
        }
    }
}


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
    dtrttp(uplo, N, Sq, Tr);
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
    ztrttp(uplo, N, Sq, Tr);
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
    dtpttr(uplo, N, Tr, Sq);
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
    ztpttr(uplo, N, Tr, Sq);
}

#if 0
template void TransposeMatrix(double *, int, int);
template void TransposeMatrix(std::complex<double> *, int, int);

// Inplace matrix transpose with conjugation for complex data
template <typename T> void TransposeMatrix(T *A, int n, int m)
{
    double alpha = 1.0;
    char *trans = "t";
    char * ordering = "R";
    if(typeid(T) == typeid(double))
    {
        dimatcopy(ordering, trans, &m, &n, &alpha, (double *)A, &n, &m);
    }
    else if(typeid(T) == typeid(std::complex<double>))
    {
        zimatcopy(ordering, trans, &m, &n, &alpha, (double *)A, &n, &m);
    }
    else
    {
        rmg_error_handler(__FILE__, __LINE__, " Data type not implemented.\n");
    }
}
#endif
