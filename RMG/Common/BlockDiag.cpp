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
#include <vector>
#include <utility>
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
#include "RmgMatrix.h"

#include "transition.h"
#include "blas.h"



template void Kpoint<double>::BlockDiag(double *vtot, double *vxc_psi);
template void Kpoint<std::complex<double>>::BlockDiag(double *vtot, double *vxc_psi);


template void Kpoint<double>::BlockDiagInternal(double *vtot, double *vxc_psi, int first, int N, double *hr, double *sr, double *vr);
template void Kpoint<std::complex<double>>::BlockDiagInternal(double *vtot, double *vxc_psi, int first, int N, std::complex<double> *hr, std::complex<double> *sr, std::complex<double> *vr);

template <class KpointType> void Kpoint<KpointType>::BlockDiag(double *vtot, double *vxc_psi)
{
    RmgTimer RT0("6-BlockDiag"), *RT1;

    // Distance between blocks of eigenvalues in eV
    double gap_thr = 20.0;

    // Required size of block before we treat it as distinct object
    int count_thr = std::floor(0.05*(double)this->nstates);

    // Identify the start and size of each block
    std::vector<std::pair<int, int>> gaps;
    double last_eig = 1.0e20;
    int count = 0, Nmax = 0, start = 0;
    for(int st=0;st < this->nstates;st++)
    {
        double eig = Ha_eV * this->Kstates[st].eig[0];
        if(((eig - last_eig) > gap_thr) && (count > count_thr))
        {
            gaps.push_back(std::make_pair(start, count));
            start = st;
            Nmax = std::max(Nmax, count);
            count = 0;
        }
        last_eig = eig;
        count++;
    }
    gaps.push_back(std::make_pair(start, this->nstates - start));
    Nmax = std::max(Nmax, this->nstates - start);

    KpointType *hr=NULL, *sr=NULL, *vr=NULL;
    if(ct.subdiag_driver != SUBDIAG_SCALAPACK && ct.subdiag_driver != SUBDIAG_ELPA)
    {
#if CUDA_ENABLED || HIP_ENABLED || SYCL_ENABLED
        KpointType *hr = (KpointType *)GpuMallocHost(Nmax * Nmax * sizeof(KpointType));
        KpointType *sr = (KpointType *)GpuMallocHost(Nmax * Nmax * sizeof(KpointType));
        KpointType *vr = (KpointType *)GpuMallocHost(Nmax * Nmax * sizeof(KpointType));
#else
        KpointType *hr = new KpointType[Nmax * Nmax]();
        KpointType *sr = new KpointType[Nmax * Nmax]();
        KpointType *vr = new KpointType[Nmax * Nmax]();
#endif
    }


    // Loop over blocks.
    for(auto &gap: gaps)
    {
        if(pct.gridpe==0)printf("\nGap start and size  %d  %d\n",gap.first, gap.second);
        this->BlockDiagInternal(vtot, vxc_psi, gap.first, gap.second, hr, sr, vr);
        if(gaps.size() > 1)
        {
            RmgTimer RT2("6-BlockDiag: ortho");
            DavidsonOrtho(gap.first, gap.second, pbasis_noncoll,
                          this->orbital_storage, ct.davidson_2stage_ortho);
        }
    }

    if(ct.subdiag_driver == SUBDIAG_SCALAPACK || ct.subdiag_driver == SUBDIAG_ELPA) return;

#if CUDA_ENABLED || HIP_ENABLED || SYCL_ENABLED
    GpuFreeHost(vr);
    GpuFreeHost(sr);
    GpuFreeHost(hr);
#else
    delete [] vr;
    delete [] sr;
    delete [] hr;
#endif

}

template <class KpointType> void Kpoint<KpointType>::BlockDiagInternal(double *vtot, double *vxc_psi, int first, int N, KpointType *hr, KpointType *sr, KpointType *vr)
{

    RmgTimer *RT1;

    KpointType alpha(1.0);
    KpointType beta(0.0);
    KpointType *newsint;

    // Use the upper part of the wavefunction array as workspace
    KpointType *h_psi = this->orbital_storage + (size_t)this->nstates * (size_t)pbasis_noncoll;

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

    RT1 = new RmgTimer("6-BlockDiag: apply hamiltonian");
    ApplyHamiltonianBlock<KpointType> (this, first, N, h_psi, vtot, vxc_psi);
    delete RT1;

    if(ct.subdiag_driver == SUBDIAG_SCALAPACK || ct.subdiag_driver == SUBDIAG_ELPA)
    {
        Scalapack *SP;
        static std::unordered_map<uint64_t, Scalapack *> SPS;
        uint64_t hash = first << 32 + N;
        if(SPS.contains(hash))
        {
            SP = SPS[hash];
        }
        else
        {
            int last = !ct.use_folded_spectrum;
            SP = new Scalapack(ct.subdiag_groups, pct.thisimg, ct.images_per_node, N,
                    ct.scalapack_block_factor, last, pct.grid_comm);
            SPS.insert({hash, SP});
        }
        Subdiag_Scalapack(this, h_psi, first, N, *SP, true);
        return;
    }

    // Update the reduced Hamiltonian and S matrices
    RT1 = new RmgTimer("6-BlockDiag: matrix setup/reduce");
    RmgGemm(trans_a, trans_n, N, N, pbasis_noncoll, alphavel, psi, pbasis_noncoll, &h_psi[first*pbasis_noncoll], pbasis_noncoll, beta, hr, N);
    PackSqToTr("U", N, hr, N, vr);
    BlockAllreduce((double *)vr, (size_t)(N+2)*N*factor/2, pct.grid_comm);
    UnPackSqToTr("U", N, hr, N, vr);

    RmgGemm(trans_a, trans_n, N, N, pbasis_noncoll, alphavel, psi, pbasis_noncoll, s_psi, pbasis_noncoll, beta, sr, N);
    PackSqToTr("U", N, sr, N, vr);
    BlockAllreduce((double *)vr, (size_t)N*(size_t)N * (size_t)factor, pct.grid_comm);
    UnPackSqToTr("U", N, sr, N, vr);
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


    // Rotate orbitals and use h_psi for scratch space
    RT1 = new RmgTimer("6-BlockDiag: rotate orbitals");
    RmgGemm(trans_n, trans_n, pbasis_noncoll, N, N, alpha, psi, pbasis_noncoll, vr, N, beta, h_psi, pbasis_noncoll);
    for(int idx=0;idx < N*pbasis_noncoll;idx++)psi[idx] = h_psi[idx];

    delete RT1;

    // Copy eigs from compact array back into state structure
    for(int st = first;st < N+first;st++) this->Kstates[st].eig[0] = eigs[st-first];


    delete [] eigs;

}

