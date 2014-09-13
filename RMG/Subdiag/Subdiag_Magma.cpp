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

#if GPU_ENABLED
    #include <cuda.h>
    #include <cuda_runtime_api.h>
    #include <cublas_v2.h>

    #if MAGMA_LIBS
        #include <magma/magma.h>
    #endif
#endif

int rmg_dsygvd_gpu(int n, double *a, int lda, double *b, int ldb,
                double *w, double *work, int lwork, int *iwork, int liwork, double *wa);
int rmg_zhegvd_gpu(int n, std::complex<double> *a, int lda, std::complex<double> *b, int ldb,
                double *eigs, double *work, int lwork, double *rwork, int lrwork, int *iwork, int liwork, double *wa);



template void Subdiag_Magma<double> (Kpoint<double> *kptr, double *Aij, double *Bij, double *Sij, double *eigs, double *eigvectors, double *gpu_eigvectors);
template void Subdiag_Magma<std::complex<double> > (Kpoint<std::complex<double>> *kptr, std::complex<double> *Aij, std::complex<double> *Bij, std::complex<double> *Sij, double *eigs, std::complex<double> *eigvectors, std::complex<double> *gpu_eigvectors);

// eigvectors holds Bij on input and the eigenvectors of the matrix on output
template <typename KpointType>
void Subdiag_Magma (Kpoint<KpointType> *kptr, KpointType *Aij, KpointType *Bij, KpointType *Sij, double *eigs, KpointType *eigvectors, KpointType *gpu_eigvectors)
{

#if !MAGMA_LIBS
    rmg_printf("This version of RMG was not built with MAGMA support. Redirecting to LAPACK.");
    Subdiag_Lapack(kptr, Aij, Bij, Sij, eigs, eigvectors);
    return;
#endif

#if MAGMA_LIBS
#if GPU_ENABLED
    KpointType ONE_t(1.0);
    int num_states = kptr->nstates;
    int ione = 1;
    int factor = 1;
    if(!ct.is_gamma) factor=2;

    cublasStatus_t custat;


    KpointType *gpuAij = (KpointType *)GpuMalloc(num_states * num_states * sizeof(KpointType));
    KpointType *gpuBij = (KpointType *)GpuMalloc(num_states * num_states * sizeof(KpointType));
    KpointType *gpuCij = gpu_eigvectors;
    KpointType *gpuSij = (KpointType *)GpuMalloc(num_states * num_states * sizeof(KpointType));
    KpointType *Cij = new KpointType[num_states * num_states]();


    if(!ct.norm_conserving_pp || (ct.norm_conserving_pp && ct.discretization == MEHRSTELLEN_DISCRETIZATION)) {


        // Create unitary matrix
        for (int idx = 0; idx < num_states; idx++) {
            Cij[idx * num_states + idx] = ONE_t;
        }

        // Transfer it to the GPU
        custat = cublasSetVector(num_states * num_states , sizeof(KpointType), Cij, ione, gpuCij, ione );

        // Transfer eigvectors which holds Bij to the gpuBij
        custat = cublasSetVector(num_states * num_states , sizeof(KpointType), eigvectors, ione, gpuBij, ione );

        // Transfer Aij to gpuAij
        custat = cublasSetVector(num_states * num_states , sizeof(KpointType), Aij, ione, gpuAij, ione );

        // Transfer Sij to gpuSij
        custat = cublasSetVector(num_states * num_states , sizeof(KpointType), Sij, ione, gpuSij, ione );

        // Invert Bij
        int *ipiv = new int[2*num_states];
        for(int idx=0;idx < num_states;idx++) ipiv[idx] = 0;
        int info = 0;
        if(ct.is_gamma) {

            // Inverse of B should be in Cij
            RmgTimer *RT1 = new RmgTimer("Diagonalization: Invert Bij");
            magma_dgesv_gpu (num_states, num_states, (double *)gpuBij, num_states, ipiv, (double *)gpuCij, num_states, &info);
            delete(RT1);


        }
        else {

            // Inverse of B should be in Cij
            RmgTimer *RT1 = new RmgTimer("Diagonalization: Invert Bij");
            magma_zgesv_gpu (num_states, num_states, (magmaDoubleComplex *)gpuBij, num_states, ipiv, (magmaDoubleComplex *)gpuCij, num_states, &info);
            delete(RT1);


        }
        if (info) {
            rmg_printf ("\n PE %d: p{d,z}gesv failed, info is %d", pct.gridpe, info);
            rmg_error_handler (__FILE__, __LINE__, " p{d,z}gesv failed");
        }
        delete [] ipiv;


        /*Multiply inverse of B and and A */
        /*B^-1*A */
        KpointType alpha(1.0);
        KpointType beta(0.0);;

        RmgTimer *RT1 = new RmgTimer("Diagonalization: matrix setup");
        RmgGemm ("n", "n", num_states, num_states, num_states, alpha,
                        Cij, num_states, Aij, num_states, beta, Bij,
                        num_states, gpuCij, gpuAij, gpuBij, false, false, false, false);

        /*Multiply the result with Sij, result is in Cij */
        RmgGemm ("n", "n", num_states, num_states, num_states, alpha,
                        Sij, num_states, Bij, num_states, beta, gpuCij,
                        num_states, gpuSij, gpuBij, gpuCij, false, false, false, false);
        delete(RT1);

    }
    else {

        // For norm conserving S=B so no need to invert and S*(B-1)*A=A so just copy A into gpuCij
        // to pass to eigensolver. Also need Sij on GPU
        custat = cublasSetVector(num_states * num_states , sizeof(KpointType), Aij, ione, gpuCij, ione );

        // Transfer Sij to gpuSij
        custat = cublasSetVector(num_states * num_states , sizeof(KpointType), Sij, ione, gpuSij, ione );

    }

    RmgTimer *RT1 = new RmgTimer("Diagonalization: magma");

    int *ifail = new int[num_states];
    int lwork = 3 * num_states * num_states + 8 * num_states;
    int liwork = 6 * num_states + 4;
    int eigs_found;
    double *work = new double[2*lwork];
    int *iwork = new int[2*liwork];
    double vx = 0.0;
    double tol = 1e-15;

    if(ct.is_gamma) {

        int info = rmg_dsygvd_gpu(num_states, (double *)gpuCij, num_states, (double *)gpuSij, num_states,
                              eigs, work, lwork, iwork, liwork, (double *)Cij);
        // We have to transfer this back in order to rotate the betaxpsi
        custat = cublasGetVector(num_states * num_states, sizeof( KpointType ), gpuCij, 1, eigvectors, 1 );
        RmgCudaError(__FILE__, __LINE__, custat, "Problem transferring eigenvector matrix from GPU to system memory.");

    }
    else {

        int lrwork = 2 * num_states * num_states + 6 * num_states;
        double *rwork = new double[2 * lrwork];
        int info = rmg_zhegvd_gpu(num_states, (std::complex<double> *)gpuCij, num_states, (std::complex<double> *)gpuSij, num_states,
                              eigs, work, lwork, rwork, lrwork, iwork, liwork, (double *)Cij);
        // We have to transfer this back in order to rotate the betaxpsi
        custat = cublasGetVector(num_states * num_states, sizeof( KpointType ), gpuCij, 1, eigvectors, 1 );
        RmgCudaError(__FILE__, __LINE__, custat, "Problem transferring eigenvector matrix from GPU to system memory.");
        delete [] rwork;

    }

    delete [] iwork;
    delete [] work;
    delete [] ifail;


    delete(RT1);
    delete(Cij);

    GpuFree(gpuSij);
    GpuFree(gpuBij);
    GpuFree(gpuAij);

#endif
#endif
     
}

#if GPU_ENABLED
#if MAGMA_LIBS
// GPU specific versions with itype=1,jobz=v,uplo=l
// assumes that a and b are already in gpu memory.
// Magma does not provide a routine that works as
// required so we put one together using magma routines
// and the cublas version of dtrsm.
int rmg_dsygvd_gpu(int n, double *a, int lda, double *b, int ldb,
                double *eigs, double *work, int lwork, int *iwork, int liwork, double *wa)
{
        int itype=1, info=0;
        double rone = 1.0;
        cublasOperation_t cu_transT = CUBLAS_OP_T, cu_transN = CUBLAS_OP_N;
        cublasSideMode_t side=CUBLAS_SIDE_LEFT;
        cublasFillMode_t cuplo=CUBLAS_FILL_MODE_LOWER;
        cublasDiagType_t diag=CUBLAS_DIAG_NON_UNIT;

        //  Form a Cholesky factorization of B.
        //  This routine is buggy
        magma_dpotrf_gpu(MagmaLower, n, b, ldb, &info);
        //cublasGetVector(n*n, sizeof( double ), b, ione, wa, ione );
        //dpotrf_("L", &n, wa, &n, &info);
        //cublasSetVector(n*n, sizeof( double ), wa, ione, b, ione );
        if( info != 0 ) {
                rmg_error_handler(__FILE__, __LINE__, "dpotrf failure");
        }

        //  Transform problem to standard eigenvalue problem and solve.
        //   dsygst_( &itype, uplo, &n, a, &lda, b, &ldb, &info );
        //   dsyevd_( jobz, uplo, &n, a, &lda, w, work, &lwork, iwork, &liwork, &info );

        magma_dsygst_gpu(itype, MagmaLower, n, a, lda, b, ldb, &info);
        if( info != 0 ) {
            rmg_error_handler(__FILE__, __LINE__, "dsygst failure");
        }

        magma_dsyevd_gpu(MagmaVec, MagmaLower, n, a, lda, eigs,
                        wa,  n,
                        work, lwork,
                        iwork, liwork,
                        &info);

        if( info != 0 ) {
                rmg_error_handler(__FILE__, __LINE__, "dsyevd failure");
        }

        //   For A*x=(lambda)*B*x and A*B*x=(lambda)*x;
        //        backtransform eigenvectors: x = inv(L)**T*y or inv(U)*y
        //   dtrsm_( "Leftx", uplo, trans, "Non-unit", &n, &n, &rone, b, &ldb, a, &lda );
        //

        cublasDtrsm_v2 (ct.cublas_handle,side, cuplo, cu_transT, diag, n, n, &rone, b, ldb, a, lda);


        return 0;
}
int rmg_zhegvd_gpu(int n, std::complex<double> *a, int lda, std::complex<double> *b, int ldb,
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
#endif
