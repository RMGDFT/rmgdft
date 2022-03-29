/*
 *
 * Copyright (c) 2018, Emil Briggs
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the <organization> nor the
 *       names of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.

 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * 
*/

#ifndef RMG_RmgMatrix_H
#define RMG_RmgMatrix_H 1

#if __cplusplus

#include <complex>
template <typename DataType> void InvertMatrix(DataType *A, DataType *B, int n);
void DsyevjDriver(double *A, double *eigs, double *work, int worksize, int n, int ld);
void DsyevdDriver(double *A, double *eigs, double *work, int worksize, int n, int ld);
void DsygvdDriver(double *A, double *B, double *eigs, double *work, int worksize, int n, int ld);
void DsygvjDriver(double *A, double *B, double *eigs, double *work, int worksize, int n, int ld);
void ZhegvdDriver(std::complex<double> *A, std::complex<double> *B, double *eigs, double *work, int worksize, int n, int ld);
void DgetrftrsDriver(int n, int m, double *A, double *B);
void DsygvdMgDriver(double *A, double *B, double *eigs, int n);

void PackSqToTr(char *uplo, int N, double *Sq, int lda, double *Tr);
void PackSqToTr(char *uplo, int N, std::complex<double> *Sq, int lda, std::complex<double> *Tr);
void UnPackSqToTr(char *uplo, int N, double *Sq, int lda, double *Tr);
void UnPackSqToTr(char *uplo, int N, std::complex<double> *Sq, int lda, std::complex<double> *Tr);

#endif
#endif


