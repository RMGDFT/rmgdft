/*
 *
 * Copyright (c) 2014, Emil Briggs
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

#ifndef RMG_RmgGemm_H
#define RMG_RmgGemm_H 1

#if __cplusplus

#include <complex>
template <typename DataType> void RmgGemm(char *transa, char *transb, int m, int n, int k,
                             DataType alpha, DataType *A, int lda, DataType *B, int ldb, DataType beta,
                             DataType *C, int ldc);

template <typename DataType> void RmgSymm(char *side, char *uplo, int m, int n,
                             DataType alpha, DataType *A, int lda, DataType *B, int ldb, DataType beta,
                             DataType *C, int ldc);

template <typename DataType> void RmgSyrkx(char *uplo, char *trans, int n, int k,
                             DataType alpha, DataType *A, int lda, DataType *B, int ldb, DataType beta,
                             DataType *C, int ldc);
void MyZgemm(char *transa, char *transb, int m, int n, int k, std::complex<double> *alpha, 
        std::complex<double> *A, int lda, std::complex<double> *B, int ldb,
        std::complex<double> *beta, std::complex<double> *C, int ldc);

template <typename DataType> void RmgGemmStridedBatched(char *transa, char *transb, int m, int n, int k,
                             DataType alpha, DataType *A, int lda, size_t strdedA, DataType *B, int ldb, size_t strideB, DataType beta,
                             DataType *C, int ldc, size_t strideC, int batchCount);
#endif
#endif





