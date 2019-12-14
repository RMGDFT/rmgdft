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
#include "RmgMatrix.h"
#include "GpuAlloc.h"
#include "blas.h"

#include "common_prototypes.h"
#include "common_prototypes1.h"
#include "transition.h"


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
    bool use_folded = ((ct.use_folded_spectrum && (ct.scf_steps > 6)) || (ct.use_folded_spectrum && (ct.runflag == RESTART)));

#if GPU_ENABLED
    // If GPU run redirect to GPU routines.
    return Subdiag_Cusolver(kptr, Aij, Bij, Sij, eigs, eigvectors);
#endif

#if SCALAPACK_LIBS
    // For folded spectrum start with scalapack if available since lapack is so slow on larger problems
    if(ct.use_folded_spectrum && (ct.scf_steps < 6)  && (ct.runflag != RESTART))
        return Subdiag_Scalapack (kptr, Aij, Bij, Sij, eigs, eigvectors);
#endif

    // Lapack is not parallel across MPI procs so only have the local master proc on a node perform
    // the diagonalization. Then broadcast the eigenvalues and vectors to the remaining local nodes.
    // If folded spectrum is selected we only want the local master to participate on each node as
    // long as there are at least 12 nodes.
    int nodes = pct.grid_npes / pct.procs_per_host;

    if(pct.is_local_master || (use_folded && (nodes < 12))) {

        // Increase the resources available to this proc since the others on the local node
        // will be idle
        int nthreads = ct.OMP_THREADS_PER_NODE;
        if((pct.procs_per_host > 1) && !use_folded) nthreads = pct.ncpus;
        omp_set_num_threads(nthreads);
#if OPENBLAS_SET_NUM_THREADS
        openblas_set_num_threads(nthreads);
#endif


        // Copy Aij into eigvectors
        for(int ix=0;ix < num_states*num_states;ix++) eigvectors[ix] = Aij[ix];

        RmgTimer *RT1 = new RmgTimer("4-Diagonalization: dsygvx/zhegvx/folded");
        int lwork = 2 * num_states * num_states + 6 * num_states + 2;
        lwork = std::max(lwork, 128000);
        int liwork = 6*num_states;
        double *work2 = new double[2*lwork];
        int *iwork = new int[liwork];

        if(ct.is_gamma) {

            if(use_folded) {
                FoldedSpectrum<double> (kptr->G, num_states, (double *)eigvectors, num_states, (double *)Sij, num_states, (double *)Aij, (double *)Bij, eigs, work2, lwork, iwork, liwork, SUBDIAG_LAPACK);

            }
            else {

                  DsygvdDriver((double *)eigvectors, (double *)Sij, eigs, work2, lwork, num_states);

            }

        }
        else {

            ZhegvdDriver((std::complex<double> *)eigvectors, (std::complex<double> *)Sij, eigs, work2, lwork, num_states);

        }

        delete [] iwork;
        delete [] work2;

        delete(RT1);

        if (info) {
            rmg_printf ("\n Lapack eigensolver failed, info is %d", info);
            rmg_error_handler (__FILE__, __LINE__, "Lapack eigensolver failed");
        }

        // Reset omp_num_threads
        omp_set_num_threads(ct.OMP_THREADS_PER_NODE);
#if OPENBLAS_SET_NUM_THREADS
        openblas_set_num_threads(ct.OMP_THREADS_PER_NODE);
#endif


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
