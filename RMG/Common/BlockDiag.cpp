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
#include <cmath>
#include <float.h>
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
#include "Mgrid.h"
#include "RmgException.h"
#include "Subdiag.h"
#include "Solvers.h"
#include "GpuAlloc.h"
#include "ErrorFuncs.h"
#include "RmgParallelFft.h"
#include "TradeImages.h"
#include "packfuncs.h"

#include "transition.h"
#include "blas.h"



template void Kpoint<double>::BlockDiag(double *vtot, double *vxc_psi);
template void Kpoint<std::complex<double>>::BlockDiag(double *vtot, double *vxc_psi);


template void Kpoint<double>::BlockDiagInternal(double *vtot, double *vxc_psi, int first, int N, double *hr, double *sr, double *vr, double *h_psi);
template void Kpoint<std::complex<double>>::BlockDiagInternal(double *vtot, double *vxc_psi, int first, int N, std::complex<double> *hr, std::complex<double> *sr, std::complex<double> *vr, std::complex<double> *h_psi);

template <class KpointType> void Kpoint<KpointType>::BlockDiag(double *vtot, double *vxc_psi)
{
    size_t N = 120;
#if CUDA_ENABLED || HIP_ENABLED || SYCL_ENABLED
    KpointType *h_psi = (KpointType *)RmgMallocHost(pbasis_noncoll * this->nstates * sizeof(KpointType));
    KpointType *hr = (KpointType *)GpuMallocHost(N * N * sizeof(KpointType));
    KpointType *sr = (KpointType *)GpuMallocHost(N * N * sizeof(KpointType));
    KpointType *vr = (KpointType *)GpuMallocHost(N * N * sizeof(KpointType));
#else
    KpointType *h_psi = new KpointType[pbasis_noncoll * this->nstates];
    KpointType *hr = new KpointType[N * N]();
    KpointType *sr = new KpointType[N * N]();
    KpointType *vr = new KpointType[N * N]();
#endif

    this->BlockDiagInternal(vtot, vxc_psi, 0, 60, hr, sr, vr, h_psi);
    this->BlockDiagInternal(vtot, vxc_psi, 60, 119, hr, sr, vr, h_psi);
    DavidsonOrtho(0, this->nstates, pbasis_noncoll, this->orbital_storage);

#if CUDA_ENABLED || HIP_ENABLED || SYCL_ENABLED
    GpuFreeHost(vr);
    GpuFreeHost(sr);
    GpuFreeHost(hr);
    RmgFreeHost(h_psi);
#else
    delete [] vr;
    delete [] sr;
    delete [] hr;
    delete [] h_psi;
#endif

}

template <class KpointType> void Kpoint<KpointType>::BlockDiagInternal(double *vtot, double *vxc_psi, int first, int N, KpointType *hr, KpointType *sr, KpointType *vr, KpointType *h_psi)
{

    RmgTimer RT0("6-BlockDiag"), *RT1;

    KpointType alpha(1.0);
    KpointType beta(0.0);
    KpointType *newsint;

    KpointType *weight = this->nl_weight;
#if HIP_ENABLED || CUDA_ENABLED
    weight = this->nl_weight_gpu;
#endif

    int max_states = this->nstates;

    char *trans_t = "t";
    char *trans_n = "n";
    char *trans_c = "c";
    char *trans_a;
    if(typeid(KpointType) == typeid(std::complex<double>)) {
         trans_a = trans_c;
    }
    else {
        trans_a = trans_t;
    }
    double vel = this->L->get_omega() / 
                 ((double)((size_t)this->G->get_NX_GRID(1) * (size_t)this->G->get_NY_GRID(1) * (size_t)this->G->get_NZ_GRID(1)));
    KpointType alphavel(vel);

    // For MPI routines
    int factor = 2;
    if(ct.is_gamma) factor = 1;

    double *eigs = new double[max_states];

    KpointType *s_psi = &this->ns[first*pbasis_noncoll];
    if(ct.norm_conserving_pp && ct.is_gamma)
        s_psi = &this->orbital_storage[first*pbasis_noncoll];;

    // short version
    KpointType *psi = &this->orbital_storage[first*pbasis_noncoll];

    // Apply Hamiltonian to the new vectors
    RT1 = new RmgTimer("6-BlockDiag: Betaxpsi");
    newsint = this->newsint_local + first * this->BetaProjector->get_num_nonloc_ions() * ct.max_nl * ct.noncoll_factor;
    this->BetaProjector->project(this, newsint, first*ct.noncoll_factor, N*ct.noncoll_factor, weight);
    delete RT1;

    if(ct.ldaU_mode != LDA_PLUS_U_NONE)
    {   
        RmgTimer RTL("6-BlockDiag: ldaUop x psi"); 
        newsint = this->orbitalsint_local + first * this->OrbitalProjector->get_num_nonloc_ions() * 
            this->OrbitalProjector->get_pstride() * ct.noncoll_factor;
        LdaplusUxpsi(this, first, N, newsint);
    }
    RT1 = new RmgTimer("6-BlockDiag: apply hamiltonian");
    ApplyHamiltonianBlock<KpointType> (this, first, N, h_psi, vtot, vxc_psi);
    delete RT1;


    // Update the reduced Hamiltonian and S matrices
    RT1 = new RmgTimer("6-BlockDiag: matrix setup/reduce");
    RmgGemm(trans_a, trans_n, N, N, pbasis_noncoll, alphavel, psi, pbasis_noncoll, &h_psi[first*pbasis_noncoll], pbasis_noncoll, beta, hr, N);

#if HAVE_ASYNC_ALLREDUCE
    // Asynchronously reduce it
    MPI_Request MPI_reqAij;
    if(ct.use_async_allreduce)
        MPI_Iallreduce(MPI_IN_PLACE, (double *)hr, N * N * factor, MPI_DOUBLE, MPI_SUM, pct.grid_comm, &MPI_reqAij);
    else
        BlockAllreduce((double *)hr, (size_t)N*(size_t)N * (size_t)factor, pct.grid_comm);
#else
    BlockAllreduce((double *)hr, (size_t)N*(size_t)N * (size_t)factor, pct.grid_comm);
#endif

    RmgGemm(trans_a, trans_n, N, N, pbasis_noncoll, alphavel, psi, pbasis_noncoll, s_psi, pbasis_noncoll, beta, sr, N);

#if HAVE_ASYNC_ALLREDUCE
    // Wait for Aij request to finish
    if(ct.use_async_allreduce) MPI_Wait(&MPI_reqAij, MPI_STATUS_IGNORE);
#endif

#if HAVE_ASYNC_ALLREDUCE
    // Asynchronously reduce Sij request
    MPI_Request MPI_reqSij;
    if(ct.use_async_allreduce)
        MPI_Iallreduce(MPI_IN_PLACE, (double *)sr, N * N * factor, MPI_DOUBLE, MPI_SUM, pct.grid_comm, &MPI_reqSij);
    else
        BlockAllreduce((double *)sr, (size_t)N*(size_t)N * (size_t)factor, pct.grid_comm);
#else
    BlockAllreduce((double *)sr, (size_t)N*(size_t)N * (size_t)factor, pct.grid_comm);
#endif

#if HAVE_ASYNC_ALLREDUCE
    // Wait for S request to finish
    if(ct.use_async_allreduce) MPI_Wait(&MPI_reqSij, MPI_STATUS_IGNORE);
#endif
    delete RT1;

    std::complex<double> *hr_C = (std::complex<double> *)hr;
    std::complex<double> *sr_C = (std::complex<double> *)sr;

    for(int i=0;i < N;i++) {
        for(int j=i+1;j < N;j++) {

            if(typeid(KpointType) == typeid(std::complex<double>))
            {
                hr_C[j + i*N] = std::conj(hr_C[i + j*N]);
                sr_C[j + i*N] = std::conj(sr_C[i + j*N]);
            }
            else
            {
                hr[j + i*N] = hr[i + j*N];
                sr[j + i*N] = sr[i + j*N];
            }
        }
    }

    RT1 = new RmgTimer("6-BlockDiag: diagonalization");
    int info = GeneralDiag(hr, sr, eigs, vr, N, N, N, ct.subdiag_driver);
    delete RT1;


    // Rotate orbitals
    RT1 = new RmgTimer("6-BlockDiag: rotate orbitals");
#if CUDA_ENABLED || HIP_ENABLED || SYCL_ENABLED
    KpointType *npsi = (KpointType *)RmgMallocHost(N*pbasis_noncoll*sizeof(KpointType));
#else
    KpointType *npsi = new KpointType[N*pbasis_noncoll];
#endif
    RmgGemm(trans_n, trans_n, pbasis_noncoll, N, N, alpha, psi, pbasis_noncoll, vr, N, beta, npsi, pbasis_noncoll);
    for(int idx=0;idx < N*pbasis_noncoll;idx++)psi[idx] = npsi[idx];

#if CUDA_ENABLED || HIP_ENABLED || SYCL_ENABLED
    RmgFreeHost(npsi);
#else
    delete [] npsi;
#endif
    delete RT1;

    // Copy eigs from compact array back into state structure
    for(int st = first;st < N+first;st++) this->Kstates[st].eig[0] = eigs[st-first];


    delete [] eigs;

}

