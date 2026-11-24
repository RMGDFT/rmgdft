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


#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include "transition.h"
#include "const.h"
#include "State.h"
#include "Kpoint.h"
#include "BaseThread.h"
#include "TradeImages.h"
#include "RmgTimer.h"
#include "RmgThread.h"
#include "rmg_reduce.h"
#include "Subdiag.h"
#include "rmgthreads.h"
#include "packfuncs.h"
#include "typedefs.h"
#include "common_prototypes.h"
#include "common_prototypes1.h"
#include "rmg_error.h"
#include "Kpoint.h"
#include "GatherScatter.h"
#include "Solvers.h"
#include "blas.h"
#include "rmg_ortho.h"


// Solver that uses multigrid preconditioning and subspace rotations
// Part of Kpoint class

template void Kpoint<double>::MgridSubspace(double *, double *vxc_psi);
template void Kpoint<std::complex<double>>::MgridSubspace(double *, double *vxc_psi);
template void Kpoint<double>::MgridSubspace(int, int, int, double *, double *vxc_psi);
template void Kpoint<std::complex<double>>::MgridSubspace(int ,int, int, double *, double *vxc_psi);
template void Kpoint<double>::MgridSubspaceBlocked(double *, double *vxc_psi);
template void Kpoint<std::complex<double>>::MgridSubspaceBlocked(double *, double *vxc_psi);


template <class KpointType> void Kpoint<KpointType>::MgridSubspace (double *vtot_psi, double *vxc_psi)
{
    MgridSubspace (0, this->nstates, ct.non_local_block_size, vtot_psi, vxc_psi);
}

template <class KpointType> void Kpoint<KpointType>::MgridSubspaceBlocked(double *vtot_psi, double *vxc_psi)
{
    int bs = ct.non_local_block_size;

    for(int is=0;is < this->nstates;is++) Kstates[is].dptr = NULL;

    int nb = this->nstates / bs;
    int irem = this->nstates % bs;
    for(int ib=0;ib < nb;ib++)
    {
            MgridSubspaceBlockedGpu (ib*bs, bs, bs, vtot_psi, vxc_psi);
    }
    if(irem)
    {
            MgridSubspaceBlockedGpu(nb*bs, irem, bs, vtot_psi, vxc_psi);
    }


    RmgTimer *RT1 = new RmgTimer("3-MgridSubspace: Diagonalization");
    if(ct.use_block_diag)
        this->BlockDiag(vtot_psi, vxc_psi);
    else
        this->Subdiag (vtot_psi, vxc_psi, ct.subdiag_driver, true);
    delete(RT1);

    // wavefunctions have changed, projectors have to be recalculated */
    RT1 = new RmgTimer("3-MgridSubspace: Beta x psi");
    this->BetaProjector->project(this, this->newsint_local, 0, nstates * ct.noncoll_factor, this->nl_weight);
    delete(RT1);

}

template <class KpointType>
void Kpoint<KpointType>::MgridSubspace (int first, int N, int bs, double *vtot_psi, double *vxc_psi)
{
    RmgTimer RT0("3-MgridSubspace"), *RT1;
    BaseThread *T = BaseThread::getBaseThread(0);
    int cfac = 1;
    if(ct.coalesce_states) cfac = pct.coalesce_factor;
    int my_pe_x, my_pe_y, my_pe_z;
    this->G->pe2xyz(pct.gridpe, &my_pe_x, &my_pe_y, &my_pe_z);
    int my_pe_offset = my_pe_x % cfac;

    KpointType *weight = this->nl_weight;
#if HIP_ENABLED || CUDA_ENABLED
    weight = this->nl_weight_gpu;
#endif
    bool potential_acceleration = (ct.potential_acceleration_constant_step > 0.0);

    int pbasis_noncoll = pbasis * ct.noncoll_factor;

    double *nvtot_psi = vtot_psi;;
    if(cfac > 1)
    {
        nvtot_psi = new double[pbasis * cfac];
        GatherGrid(this->G, pbasis, vtot_psi, nvtot_psi);
    }

    // Set trade images coalesce_factor
    this->T->set_coalesce_factor(cfac);

    int active_threads = rmg::get_active_threads();

    // When grid coalescing is enabled we need nstates to be an integral
    // multiple of (active_threads * pct.coalesce_factor) in order for the
    // coalescing routines to work correctly. So we pad nstates to satisfy
    // that condition where required (stored in mstates).
    int mstates = N / (active_threads * cfac);
    if(N % (active_threads * cfac)) mstates++;
    mstates = mstates * (active_threads * cfac);

    rmg::ortho<KpointType> MGOrtho(this->nstates, pbasis_noncoll);

    // We adjust the block size here for threading and coalescing
    int block_size = bs;
    block_size = block_size / (active_threads * cfac);
    block_size = block_size * (active_threads * cfac);
    int nblocks = mstates / block_size;
    int irem = mstates % block_size;
    if(irem) nblocks++;
    for(int is=first;is < N+first;is++) Kstates[is].eig[1] = Kstates[is].eig[0];

    std::vector<double> deig(20);
    std::fill(deig.begin(), deig.end(), 0.0);
    if(ct.use_rmm_diis)
        for(int is=first;is < N+first;is++) 
            Kstates[is].dptr = new diis<KpointType>(4, cfac*pbasis_noncoll);
std::barrier<> bbarier(active_threads);
    for(int vcycle = 0;vcycle < ct.eig_parm.mucycles;vcycle++)
    {

        // Zero out dvh array if potential acceleration is enabled
        if(potential_acceleration)
        {
           int stop = ct.ndvh * pbasis * cfac;
           for(int i=0;i < stop;i++) this->dvh[i] = 0.0;
           PotentialAccelerationReset(my_pe_offset*active_threads + ct.dvh_skip/cfac);
        }

        // Update betaxpsi        
        RT1 = new RmgTimer("3-MgridSubspace: Beta x psi");
        KpointType *newsint = this->newsint_local + first * this->BetaProjector->get_num_nonloc_ions() *
                    ct.max_nl * ct.noncoll_factor;
        this->BetaProjector->project(this, newsint,
              first*ct.noncoll_factor, mstates * ct.noncoll_factor, weight);
        delete(RT1);

        if(ct.ldaU_mode != LDA_PLUS_U_NONE)
        {
            RmgTimer RTL("3-MgridSubspace: ldaUop x psi");
            LdaplusUxpsi(this, 0, this->nstates, this->orbitalsint_local);
        }

        for(int ib = 0;ib < nblocks;ib++)
        {
            int bofs = ib * block_size;
            RT1 = new RmgTimer("3-MgridSubspace: AppNls");
            AppNls(this, newsint, this->Kstates[first + bofs].psi, this->nv, 
                   &this->ns[(first+bofs) * pbasis_noncoll],
                   bofs, std::min(block_size, mstates - bofs));
            delete(RT1);
            for(int st1=0;st1 < block_size;st1+=active_threads*cfac)
            {
                SCF_THREAD_CONTROL thread_control;

                RT1 = new RmgTimer("3-MgridSubspace: Mg_eig");
                int istart = my_pe_offset*active_threads;
                for(int ist = 0;ist < active_threads;ist++) {
                    int sindex = bofs + st1 + ist + istart;
                    if(sindex >= mstates)
                    {
                        break;
                    }
                    else
                    {
                        thread_control.job = HYBRID_EIG;
                        thread_control.vtot = nvtot_psi;
                        thread_control.vxc_psi = vxc_psi;
                        thread_control.coarse_vtot = NULL;
                        thread_control.vcycle = vcycle;
                        thread_control.sp = &this->Kstates[first + sindex];
                        thread_control.p3 = (void *)this;
// this is a hack since we are backing the nv array up before it's start in order to make
// the code in MgEigState work properly. Fix at some point.
                        if(cfac > 1)
                            thread_control.nv = this->nv - (first+bofs) * pbasis_noncoll;
                        else
                            thread_control.nv = (void *)&this->nv[(st1 + ist + istart) * pbasis_noncoll];
                        thread_control.ns = (void *)&this->ns[(first+sindex) * pbasis_noncoll];  // ns is not blocked!
                        thread_control.basetag = this->Kstates[first + sindex].istate;
                        thread_control.extratag1 = active_threads;
                        thread_control.extratag2 = first + bofs + st1;
                        thread_control.extratag3 = st1 + ist + istart;
                        thread_control.extratag4 = potential_acceleration;
                    }
                    QueueThreadTask(ist, thread_control);
                }

                // Thread tasks are set up so run them
                if(!ct.mpi_queue_mode && active_threads) T->run_thread_tasks(active_threads);
                delete RT1;
            }
            if(!ct.mpi_queue_mode && active_threads) T->run_thread_tasks(active_threads);
        RT1 = new RmgTimer("3-MgridSubspace: Mg_eig");
            if(ct.mpi_queue_mode) T->run_thread_tasks(active_threads, Rmg_Q);
        delete RT1;

        } // end for ib

        RT1 = new RmgTimer("3-MgridSubspace: Mg_eig");
        if(!ct.mpi_queue_mode && active_threads) T->run_thread_tasks(active_threads);
        if(ct.mpi_queue_mode) T->run_thread_tasks(active_threads, Rmg_Q);
        delete RT1;

        for(int is=first;is < N+first;is++)
        {
            Kstates[is].eig[1] = Kstates[is].eig[0];
        }

        RmgTimer RTO("3-MgridSubspace: ortho");
        int stop = std::min(first+block_size, this->nstates-first);
        MGOrtho.orthogonalize(first, stop, this->orbital_storage, true);

    }

    // Set trade images coalesce factor back to 1 for other routines.
    this->T->set_coalesce_factor(1);

    if(cfac > 1)
    {
        delete [] nvtot_psi;
        // Eigenvalues are not copied to all nodes in MgEigState when using coalesced grids.
        GatherEigs(this);
    }

    /* wavefunctions have changed, projectors have to be recalculated
     * but if we are using potential acceleration and not well converged yet
     * it is counterproductive to do so */
    if(!potential_acceleration)
    {
        RT1 = new RmgTimer("3-MgridSubspace: Beta x psi");
        this->BetaProjector->project(this, this->newsint_local, 0, nstates * ct.noncoll_factor, weight);
        delete(RT1);

        if(ct.ldaU_mode != LDA_PLUS_U_NONE)
        {   
            RmgTimer RTL("3-MgridSubspace: ldaUop x psi"); 
            LdaplusUxpsi(this, 0, this->nstates, this->orbitalsint_local);
        }
    }

    if(ct.use_rmm_diis)
        for(int is=first;is < N+first;is++) delete Kstates[is].dptr;

}


#include "gpu_solver_bridge.hpp"
#include "rmg_gpu_taps.hpp"
#include "vloc_gather.hpp"


// Blocked gpu version.
template <class KpointType>
void Kpoint<KpointType>::MgridSubspaceBlockedGpu(int first, int N, int bs, double *vtot_psi, double *vxc_psi)
{
    RmgTimer RT0("3-MgridSubspace"), *RT1;
    FiniteDiff FD(&Rmg_L, ct.alt_laplacian);
//    int my_pe_x, my_pe_y, my_pe_z;
//    this->G->pe2xyz(pct.gridpe, &my_pe_x, &my_pe_y, &my_pe_z);
//    int my_pe_offset = my_pe_x;
    int NX = Rmg_G->get_NX_GRID(1);
    int NY = Rmg_G->get_NY_GRID(1);
    int NZ = Rmg_G->get_NZ_GRID(1);
    double key = Rmg_G->get_hxgrid(1);
    double hx = Rmg_L.get_xside() * Rmg_G->get_hxgrid(1);
    double hy = Rmg_L.get_yside() * Rmg_G->get_hygrid(1);
    double hz = Rmg_L.get_zside() * Rmg_G->get_hzgrid(1);
    double tpiba = 2.0 * PI / L->celldm[0];
    double tpiba2 = tpiba * tpiba;


    GpuSubBox3D my_box = { Rmg_G->get_PX_OFFSET(1), Rmg_G->get_PY_OFFSET(1), Rmg_G->get_PZ_OFFSET(1), Rmg_G->get_PX0_GRID(1), Rmg_G->get_PY0_GRID(1), Rmg_G->get_PZ0_GRID(1) };

//    std::vector< GpuStencilTap > fine_taps_8th;
//    std::vector< std::vector< GpuStencilTap > > coarse_taps_2nd;
// Level 0: 8th-order fine grid
auto fine_taps_8th = build_rmg_taps_from_fd(FD, 8, hx, hy, hz, this->kp.kvec, key);

// Levels 1 .. mg_num_levels - 1: 2nd-order coarse grids (2x coarsening per level)
std::vector< std::vector< GpuStencilTap > > coarse_taps_2nd;
for (int l = 1; l < ct.eig_parm.levels+1; ++l) {
    double scale = static_cast< double >(1 << l); // 2.0, 4.0, ...
    coarse_taps_2nd.push_back(build_rmg_taps_from_fd(FD, 2, scale * hx, scale * hy, scale * hz, this->kp.kvec, key*2.0));
}

    int desired_mpi_batch = 1024;
    int pre_smooth_deg = 3;
    int post_smooth_deg = 3;
    int coarsest_solve_deg = 6;
    KpointType* h_psi_block = this->Kstates[first].psi;
    KpointType* h_vnl_block = this->nv;

    std::vector< double > global_eigvals;
    global_eigvals.resize(ct.num_states);
    for(int is=0;is < ct.num_states;is++) global_eigvals[is] = Kstates[is].eig[0];

    if(!this->mg_bridge)
    {
        this->mg_bridge = new GpuSolverBridge<KpointType>(pct.grid_comm,
                      N,
                      NX,
                      NY, 
                      NZ,
                      my_box,
                      ct.eig_parm.levels+1,
                      fine_taps_8th,
                      coarse_taps_2nd,
                      ct.lambda_min,
                      ct.lambda_max,
                      -0.0001*tpiba2,
                      get_vel(),
                      true,
                      desired_mpi_batch,
                      pre_smooth_deg,
                      post_smooth_deg,
                      coarsest_solve_deg);
    }
    this->mg_bridge->update_stencil_taps(fine_taps_8th, coarse_taps_2nd);

    KpointType *weight = this->nl_weight;
#if HIP_ENABLED || CUDA_ENABLED
    weight = this->nl_weight_gpu;
#endif
    int pbasis_noncoll = pbasis * ct.noncoll_factor;
    double *nvtot_psi = vtot_psi;;

    rmg::ortho<KpointType> MGOrtho(this->nstates, pbasis_noncoll);

    int block_size = bs;

    std::vector<double> deig(20);
    std::fill(deig.begin(), deig.end(), 0.0);

    std::vector< double > h_vloc_padded;
    gather_vloc_padded(pct.grid_comm, vtot_psi, my_box, NX, NY, NZ, h_vloc_padded, /*periodic=*/true);
    if constexpr (std::is_same_v<KpointType, double> || std::is_same_v<KpointType, std::complex<double>>)
    {
         this->mg_bridge->upload_local_potential_64(h_vloc_padded.data());
    }
    else
    {
         this->mg_bridge->upload_local_potential_32(h_vloc_padded.data());
    }

    for(int vcycle = 0;vcycle < ct.eig_parm.mucycles;vcycle++)
    {

        // Update betaxpsi        
        RT1 = new RmgTimer("3-MgridSubspace: Beta x psi");
        KpointType *newsint = this->newsint_local + first * this->BetaProjector->get_num_nonloc_ions() *
                    ct.max_nl * ct.noncoll_factor;
        this->BetaProjector->project(this, newsint,
              first*ct.noncoll_factor, N * ct.noncoll_factor, weight);
        delete(RT1);

        if(ct.ldaU_mode != LDA_PLUS_U_NONE)
        {
            RmgTimer RTL("3-MgridSubspace: ldaUop x psi");
            LdaplusUxpsi(this, 0, this->nstates, this->orbitalsint_local);
        }

        {
            RT1 = new RmgTimer("3-MgridSubspace: AppNls");
            AppNls(this, newsint, this->Kstates[first].psi, this->nv, 
                   &this->ns[first * pbasis_noncoll],
                   0, std::min(block_size, N));
            delete(RT1);

        }

        // Here insert preconditioner
        rmg::sync_device();
        RmgTimer *RT2 = new RmgTimer("3-MgridSubspace: Mg_eig");
        this->mg_bridge->solve_block(
            h_psi_block, h_vnl_block, global_eigvals,
            first, N, ct.eig_parm.mg_timestep);
        delete RT2;
        rmg::sync_device();
        for (int is = 0; is < N; is++) Kstates[first + is].eig[0] = global_eigvals[first + is];

        RmgTimer RTO("3-MgridSubspace: ortho");
        int stop = std::min(first+block_size, this->nstates-first);
        MGOrtho.orthogonalize(first, stop, this->orbital_storage, true);

    }

    /* wavefunctions have changed, projectors have to be recalculated */
    {
        RT1 = new RmgTimer("3-MgridSubspace: Beta x psi");
        KpointType *newsint = this->newsint_local + first * this->BetaProjector->get_num_nonloc_ions() *
                    ct.max_nl * ct.noncoll_factor;
        this->BetaProjector->project(this, newsint,
              first * ct.noncoll_factor, N * ct.noncoll_factor, weight);
        delete(RT1);
    }

    if(ct.ldaU_mode != LDA_PLUS_U_NONE)
    {   
        RmgTimer RTL("3-MgridSubspace: ldaUop x psi"); 
        LdaplusUxpsi(this, 0, this->nstates, this->orbitalsint_local);
    }
}
