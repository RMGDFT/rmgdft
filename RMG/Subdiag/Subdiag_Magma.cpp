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
        //#include <magmablas.h>
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

#if !GPU_ENABLED
    rmg_printf("This version of RMG was not built with GPU support so MAGMA cannot be used. Redirecting to LAPACK.");
    return Subdiag_Lapack(kptr, Aij, Bij, Sij, eigs, eigvectors);
#endif

#if MAGMA_LIBS
#if GPU_ENABLED
    KpointType ONE_t(1.0);
    KpointType ZERO_t(1.0);
    KpointType *NULLptr = NULL;


    static char *trans_t = "t";
    static char *trans_n = "n";


    int num_states = kptr->nstates;
    int ione = 1;
    bool use_folded = ((ct.use_folded_spectrum && (ct.scf_steps > 6)) || (ct.use_folded_spectrum && (ct.runflag == RESTART)));

    cublasStatus_t custat;

    // Magma is not parallel across MPI procs so only have the local master proc on a node perform
    // the diagonalization. Then broadcast the eigenvalues and vectors to the remaining local nodes.
    // If folded spectrum is selected we only want the local master to participate on each node as
    // long as there are at least 12 nodes.
    int NPES = Rmg_G->get_PE_X() * Rmg_G->get_PE_Y() * Rmg_G->get_PE_Z();
    int nodes = NPES / pct.procs_per_host;

    if(pct.is_local_master || (use_folded && (nodes < 12))) {


        if(!ct.norm_conserving_pp || (ct.norm_conserving_pp && ct.discretization == MEHRSTELLEN_DISCRETIZATION)) {

            // Invert Bij
            int *ipiv = new int[2*num_states]();
            int info = 0;
            if(ct.is_gamma) {

                // Create unitary matrix
                for(int i = 0;i < num_states*num_states;i++) eigvectors[i] = 0.0;
                for(int i = 0;i < num_states;i++) eigvectors[i + i*num_states] = 1.0;

                // Inverse of B should be in eigvectors after this call
                RmgTimer *RT1 = new RmgTimer("4-Diagonalization: Invert Bij");
                magma_dgesv (num_states, num_states, (double *)Bij, num_states, ipiv, (double *)eigvectors, num_states, &info);
                delete(RT1);


            }
            else {

                // Create unitary matrix
                KpointType Zero(0.0);
                KpointType One(0.0);
                for(int i = 0;i < num_states*num_states;i++) eigvectors[i] = Zero;
                for(int i = 0;i < num_states;i++) eigvectors[i + i*num_states] = One;

                // Inverse of B should be in eigvectors after this call
                RmgTimer *RT1 = new RmgTimer("4-Diagonalization: Invert Bij");
                magma_zgesv (num_states, num_states, (magmaDoubleComplex *)Bij, num_states, ipiv, (magmaDoubleComplex *)eigvectors, num_states, &info);
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

            RmgTimer *RT1 = new RmgTimer("4-Diagonalization: matrix setup");
            RmgGemm ("n", "n", num_states, num_states, num_states, alpha,
                            eigvectors, num_states, Aij, num_states, beta, Bij,
                            num_states, NULLptr, NULLptr, NULLptr, false, false, false, false);

            /*Multiply the result with Sij, result is in eigvectors */
            RmgGemm ("n", "n", num_states, num_states, num_states, alpha,
                            Sij, num_states, Bij, num_states, beta, eigvectors,
                            num_states, NULLptr, NULLptr, NULLptr, false, false, false, false);
            delete(RT1);

        }
        else {

            // For norm conserving S=B so no need to invert and S*(B-1)*A=A so just copy A into eigvectors
            for(int i = 0;i < num_states * num_states;i++) eigvectors[i] = Aij[i];

        }

        RmgTimer *RT1 = new RmgTimer("4-Diagonalization: magma");

        int *ifail = new int[num_states];
        int liwork = 6 * num_states + 4;
        int eigs_found;
        int *iwork = new int[2*liwork];
        double vx = 0.0;
        double tol = 1e-15;

        if(ct.is_gamma) {

            if(use_folded) {

                int lwork = num_states * num_states / 3 + num_states;
                double *work = (double *)GpuMallocHost(lwork * sizeof(KpointType));        
                FoldedSpectrum<double> (kptr->G, num_states, (double *)eigvectors, num_states, (double *)Sij, num_states, eigs, work, lwork, iwork, liwork, SUBDIAG_MAGMA);
                GpuFreeHost(work);

            }
            else {

                int info;
                int itype = 1;
                int lwork = 3 * num_states * num_states + 8 * num_states;
                double *work = (double *)GpuMallocHost(lwork * sizeof(KpointType));
                magma_dsygvd(itype, MagmaVec, MagmaLower, num_states, (double *)eigvectors, num_states, (double *)Sij, num_states, eigs, work, lwork, iwork, liwork, &info);
                GpuFreeHost(work);

            }

        }
        else {

            int info;
            int itype = 1;
            int lwork = 3 * num_states * num_states + 8 * num_states;
            double *work = (double *)GpuMallocHost(lwork * sizeof(KpointType));
            int lrwork = 2 * num_states * num_states + 6 * num_states;
            double *rwork = new double[2 * lrwork];
            magma_zhegvd(itype, MagmaVec, MagmaLower, num_states, (cuDoubleComplex *)eigvectors, num_states, (cuDoubleComplex *)Sij, num_states,
                                  eigs, (cuDoubleComplex *)work, lwork, rwork, lrwork, iwork, liwork, &info);

            GpuFreeHost(work);
        }

        delete [] iwork;
        delete [] ifail;


        delete(RT1);

        
    } // end if is_local_master

    // If only one proc on this host participated broadcast results to the rest
    if((pct.procs_per_host > 1) && !(use_folded && (nodes < 12))) {
        int factor = 2;
        if(ct.is_gamma) factor = 1;
        MPI_Bcast(eigvectors, factor * num_states*num_states, MPI_DOUBLE, 0, pct.local_comm);
        MPI_Bcast(eigs, num_states, MPI_DOUBLE, 0, pct.local_comm);

    }

    if(use_folded) return trans_t;
    return trans_n;

#endif
#endif
     
}

