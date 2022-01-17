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
#include "Gpufuncs.h"
#include "ErrorFuncs.h"
#include "blas.h"

#include "common_prototypes.h"
#include "common_prototypes1.h"
#include "transition.h"
#include "RmgMatrix.h"


#if MAGMA_LIBS
#include <magma_v2.h>


template char * Subdiag_Magma<double> (Kpoint<double> *kptr, double *Aij, double *Bij, double *Sij, double *eigs, double *eigvectors);
template char * Subdiag_Magma<std::complex<double> > (Kpoint<std::complex<double>> *kptr, std::complex<double> *Aij, std::complex<double> *Bij, std::complex<double> *Sij, double *eigs, std::complex<double> *eigvectors);

// eigvectors holds Bij on input and the eigenvectors of the matrix on output
    template <typename KpointType>
char * Subdiag_Magma (Kpoint<KpointType> *kptr, KpointType *Aij, KpointType *Bij, KpointType *Sij, double *eigs, KpointType *eigvectors)
{

#if !MAGMA_LIBS
    rmg_printf("This version of RMG was not built with Magma. Redirecting to LAPACK.");
    return Subdiag_Lapack(kptr, Aij, Bij, Sij, eigs, eigvectors);
#endif

#if HIP_ENABLED

    static char *trans_t = "t";
    static char *trans_n = "n";


    int num_states = kptr->nstates;
    bool use_folded = ((ct.use_folded_spectrum && (ct.scf_steps > 6)) || (ct.use_folded_spectrum && (ct.runflag == RESTART)));

#if SCALAPACK_LIBS
    // For folded spectrum start with scalapack if available since magma is slow on really large problems
    if(ct.use_folded_spectrum && (ct.scf_steps <= 6)  && (ct.runflag != RESTART) && (num_states > 10000))
        return Subdiag_Scalapack (kptr, Aij, Bij, Sij, eigs, eigvectors);
#endif

    RmgTimer *DiagTimer;
    static int call_count, folded_call_count;
    if(use_folded)
    {
        DiagTimer = new RmgTimer("4-Diagonalization: magma folded");
        folded_call_count++;
        rmg_printf("\nDiagonalization using folded magma for step=%d  count=%d\n\n",ct.scf_steps, folded_call_count); 
    }
    else
    {
        DiagTimer = new RmgTimer("4-Diagonalization: magma");
        call_count++;
        rmg_printf("\nDiagonalization using magma for step=%d  count=%d\n\n",ct.scf_steps, call_count); 
fflush(NULL);
    }
    
    // Magma is not parallel across MPI procs so only have the local master proc on a node perform
    // the diagonalization. Then broadcast the eigenvalues and vectors to the remaining local nodes.
    // If folded spectrum is selected we only want the local master to participate on each node as
    // long as there are at least 12 nodes.
    int nodes = pct.grid_npes / pct.procs_per_host;

    //    if(pct.is_local_master || (use_folded && (nodes < 12))) {

    // Copy A into eigvectors
//    memcpy(eigvectors, Aij, (size_t)num_states * (size_t)num_states * sizeof(KpointType));
    KpointType *eigvectors_gpu, *Sij_gpu;
    double *eigs_gpu;
    gpuMallocManaged((void **)&eigvectors_gpu, (size_t)num_states * (size_t)num_states * sizeof(KpointType));
    gpuMallocManaged((void **)&Sij_gpu, (size_t)num_states * (size_t)num_states * sizeof(KpointType));
    gpuMallocManaged((void **)&eigs_gpu, (size_t)num_states * sizeof(KpointType));
    gpuMemcpy(eigvectors_gpu, Aij, (size_t)num_states * (size_t)num_states * sizeof(KpointType), hipMemcpyDefault);
    gpuMemcpy(Sij_gpu, Sij, (size_t)num_states * (size_t)num_states * sizeof(KpointType), hipMemcpyDefault);


    RmgTimer *RT1 = new RmgTimer("4-Diagonalization: dsygvx/zhegvx/folded");

    int *ifail = new int[num_states];
    int liwork = 6 * num_states + 4;
    int *iwork = new int[2*liwork];

    if(ct.is_gamma) {

        if(use_folded) {

            int lwork = num_states * num_states / 3 + num_states;
            lwork = std::max(lwork, 128000);
            double *work, *Aij_gpu, *Bij_gpu;
            hipMalloc((void **)&work, lwork * sizeof(KpointType));
            hipMallocManaged((void **)&Aij_gpu, (size_t)num_states * (size_t)num_states * sizeof(double));
            hipMallocManaged((void **)&Bij_gpu, (size_t)num_states * (size_t)num_states * sizeof(double));
            FoldedSpectrum<double> (kptr->G, num_states, (double *)eigvectors_gpu, num_states, (double *)Sij_gpu, num_states, (double *)Aij_gpu, (double *)Bij_gpu, eigs_gpu, work, lwork, iwork, liwork, SUBDIAG_MAGMA);
            hipFree(Bij_gpu);
            hipFree(Aij_gpu);
            hipFree(work);

        }
        else {

            int lwork = 3 * num_states * num_states + 8 * num_states;
            lwork = std::max(lwork, 128000);
            double *work;
            hipMalloc((void **)&work, lwork * sizeof(KpointType));
            DsygvdDriver((double *)eigvectors_gpu, (double *)Sij_gpu, eigs_gpu, work, lwork, num_states, num_states);
            hipFree(work);

        }

    }
    else {

        int lwork = 3 * num_states * num_states + 8 * num_states;
        lwork = std::max(lwork, 128000);
        ZhegvdDriver((std::complex<double> *)eigvectors_gpu, (std::complex<double> *)Sij_gpu, eigs_gpu, NULL, lwork, num_states, num_states);

    }

    delete [] iwork;
    delete [] ifail;
    delete RT1;

    hipMemcpy(eigvectors, eigvectors_gpu, (size_t)num_states * (size_t)num_states * sizeof(KpointType), hipMemcpyDefault);
    hipMemcpy(eigs, eigs_gpu, (size_t)num_states * sizeof(double), hipMemcpyDefault);
    hipFree(eigs_gpu);
    hipFree(Sij_gpu);
    hipFree(eigvectors_gpu);

    // end if is_local_master

    // If only one proc on this host participated broadcast results to the rest
    //    if((pct.procs_per_host > 1) && !(use_folded && (nodes < 12))) {
    if(1){
        int factor = 2;
        if(ct.is_gamma) factor = 1;
        MPI_Bcast(eigvectors, factor * num_states*num_states, MPI_DOUBLE, 0, pct.local_comm);
        MPI_Bcast(eigs, num_states, MPI_DOUBLE, 0, pct.local_comm);

    }

    delete DiagTimer;
    if(use_folded) return trans_t;
    return trans_n;

#endif

}

#endif
