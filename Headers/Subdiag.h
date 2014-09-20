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

#ifndef RMG_Subdiag_H
#define RMG_Subdiag_H 1

#if __cplusplus
#include "BaseGrid.h"
#include "Lattice.h"
#include "TradeImages.h"
#include "const.h"
#include "RmgTimer.h"
#include "rmgtypedefs.h"
#include "params.h"
#include "typedefs.h"
#include "rmg_error.h"
#include "Kpoint.h"

template <typename KpointType> void Subdiag(Kpoint<KpointType> *kptr, 
                                            double *vh, 
                                            double *vnuc,
                                            double *vxc,
                                            int subdiag_driver);
template <typename KpointType>
void ApplyOperators (Kpoint<KpointType> *kptr, int istate, KpointType *a_psi, KpointType *b_psi, double *vtot);
template <typename KpointType>
void Subdiag_Lapack (Kpoint<KpointType> *kptr, KpointType *Aij, KpointType *Bij, KpointType *Sij, double *eigs, KpointType *eigvectors);
template <typename KpointType>
void Subdiag_Magma (Kpoint<KpointType> *kptr, KpointType *Aij, KpointType *Bij, KpointType *Sij, double *eigs, KpointType *eigvectors, KpointType *gpu_eigvectors);
template <typename KpointType>
void Subdiag_Scalapack (Kpoint<KpointType> *kptr, KpointType *Aij, KpointType *Bij, KpointType *Sij, double *eigs, KpointType *eigvectors);
template <typename KpointType>
int FoldedSpectrumCpu(Kpoint<KpointType> *kptr, int n, KpointType *A, int lda, KpointType *B, int ldb,
                double *eigs, double *work, int lwork, int *iwork, int liwork, KpointType *C);
void FoldedSpectrumIterator(double *A, int n, double *eigs, int k, double *X, double alpha, int iterations);
void FoldedSpectrumIterator(std::complex<double> *A, int n, double *eigs, int k, std::complex<double> *X, double alpha, int iterations, int istart);

template <typename DataType>
void FoldedSpectrumMinverse(DataType *A, DataType *Ainv, int n, int iterations);

template <typename DataType>
void FoldedSpectrumGSE(DataType *A, DataType *B, DataType *Z, int n, int istart, int istop, int *fs_eigcounts, int *fs_eigstart, int iterations);

int Rmg_dsygvd_gpu(int n, double *a, int lda, double *b, int ldb,
                double *w, double *work, int lwork, int *iwork, int liwork, double *wa);
int Rmg_zhegvd_gpu(int n, std::complex<double> *a, int lda, std::complex<double> *b, int ldb,
                double *eigs, double *work, int lwork, double *rwork, int lrwork, int *iwork, int liwork, double *wa);
void FoldedSpectrumSetup(int n, int NPES, int THISPE, 
                         int *eigstart, int *eigstop, int *eigstep,
                         int *n_start, int *n_win,
                         int *fs_eigstart, int *fs_eigstop, int *fs_eigcounts);
template <typename KpointType>
void FoldedSpectrumOrtho(int n, int eig_start, int eig_stop, int *fs_eigcounts, int *fs_eigstart, KpointType *V, KpointType *B);



#endif
#endif


