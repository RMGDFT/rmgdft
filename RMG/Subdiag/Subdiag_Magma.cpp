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

#if GPU_ENABLED
    #if MAGMA_LIBS
        #include <magma.h>
        #include <magmablas.h>
    #endif

    #include <cuda.h>
    #include <cuda_runtime_api.h>
    #include <cublas_v2.h>

#endif

int rmg_dsygvd_gpu(int n, double *a, int lda, double *b, int ldb,
                double *w, double *work, int lwork, int *iwork, int liwork, double *wa);
int rmg_zhegvd_gpu(int n, std::complex<double> *a, int lda, std::complex<double> *b, int ldb,
                double *eigs, double *work, int lwork, double *rwork, int lrwork, int *iwork, int liwork, double *wa);



template char * Subdiag_Magma<double> (Kpoint<double> *kptr, double *Aij, double *Bij, double *Sij, double *eigs, double *eigvectors);
template char * Subdiag_Magma<std::complex<double> > (Kpoint<std::complex<double>> *kptr, std::complex<double> *Aij, std::complex<double> *Bij, std::complex<double> *Sij, double *eigs, std::complex<double> *eigvectors);

// eigvectors holds Bij on input and the eigenvectors of the matrix on output
template <typename KpointType>
char * Subdiag_Magma (Kpoint<KpointType> *kptr, KpointType *Aij, KpointType *Bij, KpointType *Sij, double *eigs, KpointType *eigvectors)
{

#if !MAGMA_LIBS
    rmg_printf("This version of RMG was not built with MAGMA support. Redirecting to LAPACK.");
    return Subdiag_Lapack(kptr, Aij, Bij, Sij, eigs, eigvectors);
#endif

#if MAGMA_LIBS
#if GPU_ENABLED
    KpointType ONE_t(1.0);
    KpointType ZERO_t(1.0);

    static char *trans_t = "t";
    static char *trans_n = "n";

    int num_states = kptr->nstates;
    int ione = 1;
    bool use_folded = ((ct.use_folded_spectrum && (ct.scf_steps > 6)) || (ct.use_folded_spectrum && (ct.runflag == RESTART)));

    cublasStatus_t custat;

    // Magma is not parallel across MPI procs so only have the local master proc on a node perform
    // the diagonalization. Then broadcast the eigenvalues and vectors to the remaining local nodes.
    if(pct.is_local_master || use_folded) {

        KpointType *gpuAij = (KpointType *)GpuMalloc(num_states * num_states * sizeof(KpointType));
        KpointType *gpuBij = (KpointType *)GpuMalloc(num_states * num_states * sizeof(KpointType));
        KpointType *gpuCij = (KpointType *)GpuMalloc(num_states * num_states * sizeof(KpointType));
        KpointType *gpuSij = (KpointType *)GpuMalloc(num_states * num_states * sizeof(KpointType));
        KpointType *Cij = (KpointType *)GpuMallocHost(num_states * num_states * sizeof(KpointType));


        if(!ct.norm_conserving_pp || (ct.norm_conserving_pp && ct.discretization == MEHRSTELLEN_DISCRETIZATION)) {

            // Transfer eigvectors which holds Bij to the gpuBij
            custat = cublasSetVector(num_states * num_states , sizeof(KpointType), eigvectors, ione, gpuBij, ione );

            // Transfer Aij to gpuAij
            custat = cublasSetVector(num_states * num_states , sizeof(KpointType), Aij, ione, gpuAij, ione );

            // Transfer Sij to gpuSij
            custat = cublasSetVector(num_states * num_states , sizeof(KpointType), Sij, ione, gpuSij, ione );

            // Invert Bij
            int *ipiv = new int[2*num_states]();
            int info = 0;
            if(ct.is_gamma) {

                // Create unitary matrix on the gpu
                magmablas_dlaset(MagmaFull, num_states, num_states, 0.0, 1.0, (double *)gpuCij, num_states);

                // Inverse of B should be in Cij
                RmgTimer *RT1 = new RmgTimer("Diagonalization: Invert Bij");
                magma_dgesv_gpu (num_states, num_states, (double *)gpuBij, num_states, ipiv, (double *)gpuCij, num_states, &info);
                delete(RT1);


            }
            else {

                // Create unitary matrix on the gpu
                magmablas_zlaset(MagmaFull, num_states, num_states, MAGMA_Z_MAKE(0.0,0.0), MAGMA_Z_MAKE(1.0,0.0), (magmaDoubleComplex *)gpuCij, num_states);

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
                            gpuCij, num_states, gpuAij, num_states, beta, gpuBij,
                            num_states, gpuCij, gpuAij, gpuBij, false, false, false, false);

            /*Multiply the result with Sij, result is in Cij */
            RmgGemm ("n", "n", num_states, num_states, num_states, alpha,
                            gpuSij, num_states, gpuBij, num_states, beta, gpuCij,
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
        double *work = (double *)GpuMallocHost(lwork * sizeof(KpointType));
        int *iwork = new int[2*liwork];
        double *work2 = new double[2*lwork];
        double vx = 0.0;
        double tol = 1e-15;

        if(ct.is_gamma) {

            if(use_folded) {

                custat = cublasGetVector(num_states * num_states, sizeof( KpointType ), gpuCij, 1, eigvectors, 1 );
                custat = cublasGetVector(num_states * num_states, sizeof( KpointType ), gpuSij, 1, Sij, 1 );
                custat = cublasGetVector(num_states * num_states, sizeof( KpointType ), gpuAij, 1, Aij, 1 );
                GpuFree(gpuSij);
                GpuFree(gpuBij);
                GpuFree(gpuAij);

                FoldedSpectrum<double> ((Kpoint<double> *)kptr, num_states, (double *)eigvectors, num_states, (double *)Sij, num_states, eigs, work2, lwork, iwork, liwork, (double *)Aij, SUBDIAG_MAGMA);
//                custat = cublasSetVector(num_states * num_states , sizeof(KpointType), eigvectors, ione, gpu_eigvectors, ione );

                delete [] work2;
                delete [] iwork;
                GpuFreeHost(work);
                delete [] ifail;


                delete(RT1);
                GpuFreeHost(Cij);

                return trans_t;

            }
            else {

                int info = Rmg_dsygvd_gpu(num_states, (double *)gpuCij, num_states, (double *)gpuSij, num_states,
                                      eigs, work, lwork, iwork, liwork, (double *)Cij);
                // We have to transfer this back in order to rotate the betaxpsi
                custat = cublasGetVector(num_states * num_states, sizeof( KpointType ), gpuCij, 1, eigvectors, 1 );
                RmgCudaError(__FILE__, __LINE__, custat, "Problem transferring eigenvector matrix from GPU to system memory.");

            }

        }
        else {

            int lrwork = 2 * num_states * num_states + 6 * num_states;
            double *rwork = new double[2 * lrwork];
            int info = Rmg_zhegvd_gpu(num_states, (std::complex<double> *)gpuCij, num_states, (std::complex<double> *)gpuSij, num_states,
                                  eigs, work, lwork, rwork, lrwork, iwork, liwork, (double *)Cij);
            // We have to transfer this back in order to rotate the betaxpsi
            custat = cublasGetVector(num_states * num_states, sizeof( KpointType ), gpuCij, 1, eigvectors, 1 );
            RmgCudaError(__FILE__, __LINE__, custat, "Problem transferring eigenvector matrix from GPU to system memory.");
            delete [] rwork;

        }

        delete [] work2;
        delete [] iwork;
        GpuFreeHost(work);
        delete [] ifail;


        delete(RT1);
        GpuFreeHost(Cij);

        GpuFree(gpuSij);
        GpuFree(gpuCij);
        GpuFree(gpuBij);
        GpuFree(gpuAij);

        
    } // end if is_local_master

    // If only one proc on this host participated broadcast results to the rest
    if((pct.procs_per_host > 1) && !use_folded) {
        int factor = 2;
        if(ct.is_gamma) factor = 1;
        MPI_Bcast(eigvectors, factor * num_states*num_states, MPI_DOUBLE, 0, pct.local_comm);
        MPI_Bcast(eigs, num_states, MPI_DOUBLE, 0, pct.local_comm);

    }

    return trans_n;

#endif
#endif
     
}

