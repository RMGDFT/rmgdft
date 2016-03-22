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
    RmgTimer RT0("Davidson"), *RT1;

    OrbitalType alpha(1.0);
    OrbitalType beta(0.0);

    int pbasis = kptr->pbasis;
    int nstates = kptr->nstates;
    int notconv = nstates;
    int max_steps = 2;
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
    double *eigs = new double[ct.max_states];
    double *eigsw = new double[ct.max_states];
    bool *converged = new bool[ct.max_states]();
    OrbitalType *global_matrix1 = new OrbitalType[ct.max_states * ct.max_states];
    OrbitalType *global_matrix2 = new OrbitalType[ct.max_states * ct.max_states];
    OrbitalType *vr = new OrbitalType[ct.max_states * ct.max_states]();
    for(int idx = 0;idx < nstates;idx++) vr[idx*ct.max_states + idx] = OrbitalType(1.0);

#if GPU_ENABLED
    cublasStatus_t custat;
    OrbitalType *h_psi = (OrbitalType *)GpuMallocHost(pbasis * ct.max_states * sizeof(OrbitalType));
#else
    OrbitalType *h_psi = new OrbitalType[pbasis * ct.max_states];
#endif

    // short verstion
    OrbitalType *psi = kptr->orbital_storage;

    // Apply Hamiltonian to current set of eigenvectors. At the current time
    // kptr->ns is also computed in ApplyHamiltonianBlock by AppNls and stored in kptr->ns
    double fd_diag = ApplyHamiltonianBlock (kptr, 0, nstates, h_psi, vtot); 
    OrbitalType *s_psi = kptr->ns;

    // Copy current eigs into compact array
    for(int st1 = 0;st1 < nstates;st1++) eigs[st1] = kptr->Kstates[st1].eig[0];

#if 0
    for(int st1=0;st1 < nstates;st1++) {
        for(int idx = 0;idx < pbasis;idx++) {
            kptr->Kstates[st1+nstates].psi[idx] = (h_psi[st1*pbasis + idx] - eigs[st1] * s_psi[st1*pbasis + idx]) / 
                                    (fd_diag + vtot[idx] +kptr->nl_Bweight[idx] - eigs[st1]*kptr->nl_weight[idx]);
        }
        //kptr->Kstates[st1+nstates].normalize(kptr->Kstates[st1+nstates].psi, st1+nstates);
    }
#endif

    for(int steps = 0;steps < max_steps;steps++) {

        // Reorder eigenvectors
        int np = 0;
        for(int st = 0;st < nstates;st++) {

            if(!converged[st]) {

                if(np != st) {
                    for(int idx=0;idx < ct.max_states;idx++) vr[idx + np*ct.max_states] = vr[idx + st*ct.max_states];
                }
                eigsw[nstates + np] = eigs[st];
                np++;                

            }

        }

        int nb1 = nstates;
#if 1
        // expand the basis set with the residuals ( H - e*S )|psi>
        RmgGemm(trans_n, trans_n, pbasis, notconv, nstates, alpha, s_psi, pbasis, vr, ct.max_states, beta, &psi[nb1*pbasis], pbasis, 
                NULLptr, NULLptr, NULLptr, false, true, false, true);

        for(int st1=0;st1 < notconv;st1++) {
            for(int idx=0;idx < pbasis;idx++) psi[(st1 + nb1)*pbasis + idx] = -eigsw[nb1 + st1] * psi[(st1 + nb1)*pbasis + idx];
        }

        RmgGemm(trans_n, trans_n, pbasis, notconv, nstates, alpha, h_psi, pbasis, vr, ct.max_states, alpha, &psi[nb1*pbasis], pbasis, 
                NULLptr, NULLptr, NULLptr, false, true, false, true);



        // Apply preconditioner
        for(int st1=0;st1 < notconv;st1++) {
            for(int idx = 0;idx < pbasis;idx++) {
                psi[(st1 + nb1)*pbasis + idx] /= (fd_diag + vtot[idx] +kptr->nl_Bweight[idx] - eigsw[st1+nb1]*kptr->nl_weight[idx]);
            }
            //kptr->Kstates[st1+nstates].normalize(kptr->Kstates[st1+nstates].psi, st1+nstates);
        }

        // Should we renormalize?
#endif

        // Apply Hamiltonian to the new vectors
kptr->nstates += notconv;
ct.num_states = kptr->nstates;
Betaxpsi (kptr);
kptr->mix_betaxpsi(0);
Subdiag (kptr, vtot, ct.subdiag_driver);
        ApplyHamiltonianBlock (kptr, nb1, notconv, &h_psi[nb1*pbasis], vtot);


        // Update the reduced Hamiltonian
kptr->nstates -= notconv;
ct.num_states = kptr->nstates;




        // Copy updated eigenvalues back
        //for(int st=0;st < nstates;st++) eigs[st] = eigsw[st];

    }

    // Copy eigs from compact array back into state structure
//    for(int st = 0;st < nstates;st++) kptr->Kstates[st].eig[0] = eigs[st];




#if GPU_ENABLED
    GpuFreeHost(h_psi);
#else
    delete [] h_psi;
#endif

    delete [] vr;
    delete [] global_matrix2;
    delete [] global_matrix1;
    delete [] converged;
    delete [] eigsw;
    delete [] eigs;

}


#if 0
    // Compute A matrix
    RT1 = new RmgTimer("Davidson: matrix setup/reduce");
    OrbitalType alpha(1.0);
    OrbitalType beta(0.0);
    RmgGemm(trans_a, trans_n, nstates, nstates, pbasis, alpha, kptr->orbital_storage, pbasis, h_psi, pbasis, beta, global_matrix1, nstates, NULLptr, NULLptr, NULLptr, false, true, false, true);

#if HAVE_ASYNC_ALLREDUCE
    // Asynchronously reduce it
    MPI_Request MPI_reqAij;
    MPI_Iallreduce(MPI_IN_PLACE, (double *)global_matrix1, nstates * nstates * factor, MPI_DOUBLE, MPI_SUM, pct.grid_comm, &MPI_reqAij);
#else
    MPI_Allreduce(MPI_IN_PLACE, (double *)global_matrix1, nstates * nstates * factor, MPI_DOUBLE, MPI_SUM, pct.grid_comm);
#endif

    // Compute S matrix
    OrbitalType alpha1(vel);
    RmgGemm (trans_a, trans_n, nstates, nstates, pbasis, alpha1, kptr->orbital_storage, pbasis, kptr->ns, pbasis, beta, global_matrix2, nstates, NULLptr, NULLptr, NULLptr, false, true, false, true);

#if HAVE_ASYNC_ALLREDUCE
    // Wait for Aij request to finish
    MPI_Wait(&MPI_reqAij, MPI_STATUS_IGNORE);
#endif

#if HAVE_ASYNC_ALLREDUCE
    // Asynchronously reduce Sij request
    MPI_Request MPI_reqSij;
    MPI_Iallreduce(MPI_IN_PLACE, (double *)global_matrix2, nstates * nstates * factor, MPI_DOUBLE, MPI_SUM, pct.grid_comm, &MPI_reqSij);
#else
    MPI_Allreduce(MPI_IN_PLACE, (double *)global_matrix2, nstates * nstates * factor, MPI_DOUBLE, MPI_SUM, pct.grid_comm);
#endif

#if HAVE_ASYNC_ALLREDUCE
    // Wait for S request to finish and when done store copy in Sij
    MPI_Wait(&MPI_reqSij, MPI_STATUS_IGNORE);
#endif
    delete RT1;
#endif
    // Now we have reduced the original 
//    RmgGemm (trans_n, trans_n, pbasis, nstates, nstates, alpha, s_psi, pbasis, global_matrix1, nstates, beta, &kptr->orbital_storage[nstates*pbasis], pbasis, NULLptr, NULLptr, NULLptr, false, true, false, true);

    // Generate initial set of correction vectors. 
