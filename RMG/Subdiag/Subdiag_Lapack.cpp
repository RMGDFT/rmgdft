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
#include <omp.h>
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
#include "blas.h"

#include "prototypes.h"
#include "common_prototypes.h"
#include "common_prototypes1.h"
#include "transition.h"

#if GPU_ENABLED
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cublas_v2.h>
#endif


template char * Subdiag_Lapack<double> (Kpoint<double> *kptr, double *Aij, double *Bij, double *Sij, double *eigs, double *eigvectors);
template char * Subdiag_Lapack<std::complex<double> > (Kpoint<std::complex<double>> *kptr, std::complex<double> *Aij, std::complex<double> *Bij, std::complex<double> *Sij, double *eigs, std::complex<double> *eigvectors);

// eigvectors holds Bij on input and the eigenvectors of the matrix on output
template <typename KpointType>
char * Subdiag_Lapack (Kpoint<KpointType> *kptr, KpointType *Aij, KpointType *Bij, KpointType *Sij, double *eigs, KpointType *eigvectors)
{

    static char *trans_t = "t";
    static char *trans_n = "n";
    int info = 0;
    int num_states = kptr->nstates;
    int ione = 1;
    bool use_folded = ((ct.use_folded_spectrum && (ct.scf_steps > 6)) || (ct.use_folded_spectrum && (ct.runflag == RESTART)));

#if SCALAPACK_LIBS
    // For folded spectrum start with scalapack if available since lapack is so slow on larger problems
    if(ct.use_folded_spectrum && (ct.scf_steps < 6)  && (ct.runflag != RESTART))
        return Subdiag_Scalapack (kptr, Aij, Bij, Sij, eigs, eigvectors);
#endif

    // Lapack is not parallel across MPI procs so only have the local master proc on a node perform
    // the diagonalization. Then broadcast the eigenvalues and vectors to the remaining local nodes.
    // If folded spectrum is selected we only want the local master to participate on each node as
    // long as there are at least 12 nodes.
    int NPES = Rmg_G->get_PE_X() * Rmg_G->get_PE_Y() * Rmg_G->get_PE_Z();
    int nodes = NPES / pct.procs_per_host;

    if(pct.is_local_master || (use_folded && (nodes < 12))) {

        // Increase the resources available to this proc since the others on the local node
        // will be idle
        int nthreads = ct.THREADS_PER_NODE;
        if((pct.procs_per_host > 1) && !use_folded) nthreads = pct.ncpus;
        omp_set_num_threads(nthreads);


        KpointType ONE_t(1.0);
        KpointType *NULLptr = NULL;
    #if GPU_ENABLED
        KpointType ZERO_t(0.0);
        KpointType *Cij = (KpointType *)GpuMallocHost(num_states * num_states * sizeof(KpointType));
        for(int ix = 0;ix < num_states * num_states;ix++) Cij[ix] = ZERO_t;
    #else
        KpointType *Cij = new KpointType[num_states * num_states]();
    #endif

        // Create unitary matrix
        for (int idx = 0; idx < num_states; idx++) {
            Cij[idx * num_states + idx] = ONE_t;
        }

        if(!ct.norm_conserving_pp || (ct.norm_conserving_pp && ct.discretization == MEHRSTELLEN_DISCRETIZATION)) {
            // Invert Bij
            int *ipiv = new int[2*num_states]();
            if(ct.is_gamma) {

                // Inverse of B should be in Cij
                RmgTimer *RT1 = new RmgTimer("4-Diagonalization: Invert Bij");
                dgesv (&num_states, &num_states, (double *)eigvectors, &num_states, ipiv, (double *)Cij, &num_states, &info);
                delete(RT1);

            }
            else {

                // Inverse of B should be in Cij
                RmgTimer *RT1 = new RmgTimer("4-Diagonalization: Invert Bij");
                zgesv (&num_states, &num_states, (double *)eigvectors, &num_states, ipiv, (double *)Cij, &num_states, &info);
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
                            Cij, num_states, Aij, num_states, beta, Bij,
                            num_states, NULLptr, NULLptr, NULLptr,
                            true, true, true, true);

            /*Multiply the result with Sij, result is in Cij */
            RmgGemm ("n", "n", num_states, num_states, num_states, alpha,
                            Sij, num_states, Bij, num_states, beta, Cij,
                            num_states, NULLptr, NULLptr, NULLptr,
                            true, true, true, true);
            delete(RT1);

        }
        else {

            // norm-conserving and central FD
            for(int ix=0;ix < num_states*num_states;ix++) Cij[ix] = Aij[ix];

        }


        RmgTimer *RT1 = new RmgTimer("4-Diagonalization: dsygvx/zhegvx/folded");
        int *ifail = new int[num_states];
        int lwork = 2 * num_states * num_states + 6 * num_states + 2;
        int liwork = 6*num_states;
        int eigs_found;
        double *work2 = new double[2*lwork];
        int *iwork = new int[liwork];
        double vx = 0.0;
        double tol = 1e-15;

        if(ct.is_gamma) {

            if(use_folded) {
                FoldedSpectrum<double> (kptr->G, num_states, (double *)Cij, num_states, (double *)Sij, num_states, eigs, work2, lwork, iwork, liwork, SUBDIAG_LAPACK);
                for(int idx=0;idx< num_states * num_states;idx++)eigvectors[idx] = Cij[idx]; 

            }
            else {

//                dsygvx (&ione, "v", "A", "l", &num_states, (double *)Cij, &num_states, (double *)Sij, &num_states,
//                                &vx, &vx, &ione, &ione,  &tol, &eigs_found, eigs, (double *)eigvectors, &num_states, work2,
//                                &lwork, iwork, ifail, &info);
//
                  dsygvd(&ione, "V", "L", &num_states, (double *)Cij, &num_states, (double *)Sij, &num_states,
                         eigs, work2, &lwork, iwork, &liwork, &info);
                  for(int i=0;i<num_states*num_states;i++)eigvectors[i] = Cij[i];

            }

        }
        else {

            double *rwork = new double[8 * num_states];
            zhegvx (&ione, "v", "A", "l", &num_states, (double *)Cij, &num_states, (double *)Sij, &num_states,
                            &vx, &vx, &ione, &ione,  &tol, &eigs_found, eigs, (double *)eigvectors, &num_states, work2,
                            &lwork, rwork, iwork, ifail, &info);
            delete [] rwork;

        }

        delete [] iwork;
        delete [] work2;
        delete [] ifail;

        delete(RT1);
    #if GPU_ENABLED
        GpuFreeHost(Cij);
    #else
        delete [] Cij;
    #endif


        if (info) {
            rmg_printf ("\n Lapack eigensolver failed, info is %d", info);
            rmg_error_handler (__FILE__, __LINE__, "Lapack eigensolver failed");
        }

        // Reset omp_num_threads
        omp_set_num_threads(ct.THREADS_PER_NODE);

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
}
