/*
 *
 * Copyright 2014 The RMG Project Developers. See the COPYRIGHT file 
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


#include <stdio.h>
#include <stdlib.h>
#include <complex>
#include <typeinfo>
#include <string.h>

#include "blas.h"



void MyZgemm(char *transa, char *transb, int m, int n, int k, std::complex<double> *alpha, 
        std::complex<double> *A, int lda, std::complex<double> *B, int ldb, 
        std::complex<double> *beta, std::complex<double> *C, int ldc)

{

    int alloc_a, alloc_b;
    double *a_ptr, *b_ptr, a_herm, b_herm, scale;
    char *transa_d, *transb_d;
    int ione = 1, itwo = 2;
    a_herm = 1.0;
    b_herm = 1.0;

    if( (!strcmp("n", transa)) || (!strcmp("N", transa)) ) 
    {
        alloc_a = lda * k;
        transa_d = "N";
    }
    else
    {
        alloc_a = lda * m;
        transa_d = "T";
    }

    if( (!strcmp("n", transb)) || (!strcmp("N", transb)) ) 
    {
        alloc_b = ldb * n;
        transb_d = "N";
    }
    else
    {
        alloc_b = ldb * k;
        transb_d = "T";
    }

    // if hermitian, imaginary part will be negative of transpose
    if( (!strcmp("c", transa)) || (!strcmp("C", transa)) )  a_herm = -1.0;
    if( (!strcmp("c", transb)) || (!strcmp("C", transb)) )  b_herm = -1.0;

    double *Ar = new double[2*alloc_a];
    double *Ai = &Ar[alloc_a];
    double *Br = new double[2*alloc_b];
    double *Bi = &Br[alloc_b];

    double *ABr = new double[2 * m * n];
    double *ABi = &ABr[m*n];

    double zero(0.0), one(1.0);
    a_ptr = (double *)A;
    b_ptr = (double *)B;

    //  separte the complex matrix into two real matrix 
    dcopy(&alloc_a, a_ptr, &itwo, Ar, &ione);
    dcopy(&alloc_a, &a_ptr[1], &itwo, Ai, &ione);
    dcopy(&alloc_b, b_ptr, &itwo, Br, &ione);
    dcopy(&alloc_b, &b_ptr[1], &itwo, Bi, &ione);


    scale = - a_herm * b_herm; 
    dgemm(transa_d, transb_d, &m, &n, &k, &one, (double *)Ar, &lda, Br, &ldb, &zero, ABr, &m);
    dgemm(transa_d, transb_d, &m, &n, &k, &scale, (double *)Ai, &lda, Bi, &ldb, &one, ABr, &m);

    dgemm(transa_d, transb_d, &m, &n, &k, &a_herm, (double *)Ai, &lda, Br, &ldb, &zero, ABi, &m);
    dgemm(transa_d, transb_d, &m, &n, &k, &b_herm, (double *)Ar, &lda, Bi, &ldb, &one, ABi, &m);


    for(int j = 0; j < n; j++)
        for(int i = 0; i < m; i++)
        {
            double tem_r = std::real(*alpha) * ABr[j * m + i] - std::imag(*alpha) *ABi[j * m + i];
            double tem_i = std::real(*alpha) * ABi[j * m + i] + std::imag(*alpha) *ABr[j * m + i];

            C[j * ldc + i] = *beta * C[j * ldc + i] + std::complex<double>(tem_r, tem_i);
        }

    delete []Ar;
    delete []Br;
    delete []ABr;

}
