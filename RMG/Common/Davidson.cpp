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
#include "rmgthreads.h"
#include "RmgTimer.h"
#include "RmgThread.h"
#include "GlobalSums.h"
#include "Kpoint.h"
#include "RmgGemm.h"
#include "Subdiag.h"
#include "Solvers.h"
#include "GpuAlloc.h"
#include "ErrorFuncs.h"

#include "transition.h"

#if GPU_ENABLED
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cublas_v2.h>
#endif


// Davidson diagonalization solver

template void Davidson<double>(Kpoint<double> *, double *);
template void Davidson<std::complex<double> >(Kpoint<std::complex<double>> *, double *);

template <typename OrbitalType>
void Davidson (Kpoint<OrbitalType> *kptr, double *vtot)
{
    RmgTimer RT0("Diagonalization"), *RT1;
    int pbasis = kptr->pbasis;
    char *trans_t = "t";
    char *trans_n = "n";
    char *trans_c = "c";
    char *trans_a;
    if(typeid(OrbitalType) == typeid(std::complex<double>)) {
         trans_a = trans_c;
    }
    else {
        trans_a = trans_t;
    }
    double vel = kptr->L->get_omega() / 
                 ((double)(kptr->G->get_NX_GRID(1) * kptr->G->get_NY_GRID(1) * kptr->G->get_NZ_GRID(1)));

    // For MPI routines
    int factor = 2;
    if(ct.is_gamma) factor = 1;
    OrbitalType *NULLptr = NULL;

    OrbitalType *global_matrix1 = new OrbitalType[ct.max_states * ct.max_states];
    OrbitalType *global_matrix2 = new OrbitalType[ct.max_states * ct.max_states];

#if GPU_ENABLED

    cublasStatus_t custat;
    OrbitalType *h_psi = (OrbitalType *)GpuMallocHost(pbasis * ct.max_states * sizeof(OrbitalType));

#else

    OrbitalType *h_psi = new OrbitalType[pbasis * ct.max_states];

#endif

    // Apply Hamiltonian to current set of eigenvectors. At the current time
    // kptr->ns is also computed in ApplyHamiltonianBlock by AppNls and stored in kptr->ns
    ApplyHamiltonianBlock (kptr, 0, kptr->nstates, h_psi, vtot); 
 
    OrbitalType *s_psi = kptr->ns;

printf("GOT HERE 0\n");fflush(NULL);
    // Generate a new set of expanded eigenvectors from the residuals
    for(int st1 = 0;st1 < kptr->nstates;st1++) {
        for(int idx=0;idx < pbasis;idx++) {
            kptr->Kstates[st1 + kptr->nstates].psi[idx] = h_psi[st1*pbasis] - kptr->Kstates[st1].eig[0] * s_psi[st1*pbasis];
        }
    }    
printf("GOT HERE 1\n");fflush(NULL);

    // Here would be preconditioning

    kptr->nstates = 2*kptr->nstates;
    MgridSubspace(kptr, vtot);
    kptr->nstates /= 2;

#if 0
    // Compute A matrix
    RT1 = new RmgTimer("Davidson: matrix setup/reduce");
    OrbitalType alpha(1.0);
    OrbitalType beta(0.0);
    RmgGemm(trans_a, trans_n, kptr->nstates, kptr->nstates, pbasis, alpha, kptr->orbital_storage, pbasis, h_psi, pbasis, beta, global_matrix1, kptr->nstates, NULLptr, NULLptr, NULLptr, false, true, false, true);

#if HAVE_ASYNC_ALLREDUCE
    // Asynchronously reduce it
    MPI_Request MPI_reqAij;
    MPI_Iallreduce(MPI_IN_PLACE, (double *)global_matrix1, kptr->nstates * kptr->nstates * factor, MPI_DOUBLE, MPI_SUM, pct.grid_comm, &MPI_reqAij);
#else
    MPI_Allreduce(MPI_IN_PLACE, (double *)global_matrix1, kptr->nstates * kptr->nstates * factor, MPI_DOUBLE, MPI_SUM, pct.grid_comm);
#endif

    // Compute S matrix
    OrbitalType alpha1(vel);
    RmgGemm (trans_a, trans_n, kptr->nstates, kptr->nstates, pbasis, alpha1, kptr->orbital_storage, pbasis, kptr->ns, pbasis, beta, global_matrix2, kptr->nstates, NULLptr, NULLptr, NULLptr, false, true, false, true);

#if HAVE_ASYNC_ALLREDUCE
    // Wait for Aij request to finish
    MPI_Wait(&MPI_reqAij, MPI_STATUS_IGNORE);
#endif

#if HAVE_ASYNC_ALLREDUCE
    // Asynchronously reduce Sij request
    MPI_Request MPI_reqSij;
    MPI_Iallreduce(MPI_IN_PLACE, (double *)global_matrix2, kptr->nstates * kptr->nstates * factor, MPI_DOUBLE, MPI_SUM, pct.grid_comm, &MPI_reqSij);
#else
    MPI_Allreduce(MPI_IN_PLACE, (double *)global_matrix2, kptr->nstates * kptr->nstates * factor, MPI_DOUBLE, MPI_SUM, pct.grid_comm);
#endif

#if HAVE_ASYNC_ALLREDUCE
    // Wait for S request to finish and when done store copy in Sij
    MPI_Wait(&MPI_reqSij, MPI_STATUS_IGNORE);
#endif
    delete RT1;
#endif


    // Now we have reduced the original 

#if GPU_ENABLED
    GpuFreeHost(h_psi);
#else
    delete [] h_psi;
    delete [] global_matrix2;
    delete [] global_matrix1;
#endif

   exit(0); 
}


