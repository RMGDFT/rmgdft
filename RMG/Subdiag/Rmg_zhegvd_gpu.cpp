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

#include <complex>
#include "FiniteDiff.h"
#include "const.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "rmg_error.h"
#include "RmgTimer.h"
#include "GlobalSums.h"
#include "Kpoint.h"
#include "Subdiag.h"
#include "RmgGemm.h"
#include "GpuAlloc.h"
#include "ErrorFuncs.h"
#include "blas.h"

#include "prototypes.h"
#include "common_prototypes.h"
#include "common_prototypes1.h"
#include "transition.h"

#if (GPU_ENABLED && MAGMA_LIBS)
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cublas_v2.h>
#include <magma.h>


// GPU specific versions with itype=1,jobz=v,uplo=l
// assumes that a and b are already in gpu memory.
// Magma does not provide a routine that works as
// required so we put one together using magma routines
// and the cublas version of dtrsm.
int Rmg_zhegvd_gpu(int n, std::complex<double> *a, int lda, std::complex<double> *b, int ldb,
                double *eigs, double *work, int lwork, double *rwork, int lrwork, int *iwork, int liwork, double *wa)
{
        int itype=1, info=0;
        cublasOperation_t cu_transT = CUBLAS_OP_T, cu_transN = CUBLAS_OP_N, cu_transC = CUBLAS_OP_C;
        cublasSideMode_t side=CUBLAS_SIDE_LEFT;
        cublasFillMode_t cuplo=CUBLAS_FILL_MODE_LOWER;
        cublasDiagType_t diag=CUBLAS_DIAG_NON_UNIT;

        //  Form a Cholesky factorization of B.
        //  This routine is buggy
        magma_zpotrf_gpu(MagmaLower, n, (magmaDoubleComplex *)b, ldb, &info);
        //cublasGetVector(n*n, sizeof( double ), b, ione, wa, ione );
        //dpotrf_("L", &n, wa, &n, &info);
        //cublasSetVector(n*n, sizeof( double ), wa, ione, b, ione );
        if( info != 0 ) {
                rmg_error_handler(__FILE__, __LINE__, "dpotrf failure");
        }

        //  Transform problem to standard eigenvalue problem and solve.
        //   dsygst_( &itype, uplo, &n, a, &lda, b, &ldb, &info );
        //   dsyevd_( jobz, uplo, &n, a, &lda, w, work, &lwork, iwork, &liwork, &info );

        magma_zhegst_gpu(itype, MagmaLower, n, (magmaDoubleComplex *)a, lda, (magmaDoubleComplex *)b, ldb, &info);
        if( info != 0 ) {
            rmg_error_handler(__FILE__, __LINE__, "dsygst failure");
        }

        magma_zheevd_gpu(MagmaVec, MagmaLower, n, (magmaDoubleComplex *)a, lda, eigs,
                        (magmaDoubleComplex *)wa,  n,
                        (magmaDoubleComplex *)work, lwork,
                        rwork, lrwork,
                        iwork, liwork,
                        &info);

        if( info != 0 ) {
                rmg_error_handler(__FILE__, __LINE__, "dsyevd failure");
        }

        //   For A*x=(lambda)*B*x and A*B*x=(lambda)*x;
        //        backtransform eigenvectors: x = inv(L)**T*y or inv(U)*y
        //   dtrsm_( "Leftx", uplo, trans, "Non-unit", &n, &n, &rone, b, &ldb, a, &lda );
        //

        std::complex<double> zone(1.0, 0.0);
        cublasZtrsm_v2 (ct.cublas_handle,side, cuplo, cu_transC, diag, n, n, (const cuDoubleComplex *)&zone, (cuDoubleComplex *)b, ldb, (cuDoubleComplex *)a, lda);


        return 0;
}

#endif
