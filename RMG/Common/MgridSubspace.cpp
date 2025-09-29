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
#include "GlobalSums.h"
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
    bool potential_acceleration = (ct.potential_acceleration_constant_step > 0.0);
    int bs = ct.non_local_block_size;

    for(int is=0;is < this->nstates;is++) Kstates[is].dptr = NULL;

    int nb = this->nstates / bs;
    int irem = this->nstates % bs;
    for(int ib=0;ib < nb;ib++)
    {
        MgridSubspace (ib*bs, bs, bs, vtot_psi, vxc_psi);
    }
    if(irem)
    {
        MgridSubspace (nb*bs, irem, bs, vtot_psi, vxc_psi);
    }


    RmgTimer *RT1 = new RmgTimer("3-MgridSubspace: Diagonalization");
    this->Subdiag (vtot_psi, vxc_psi, ct.subdiag_driver);
//  To use BlockDiag comment out the Subdiag line above and uncomment
//  the line below. Block boundaries are set in BlockDiag.
//    this->BlockDiag(vtot_psi, vxc_psi);
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
    int my_pe_x, my_pe_y, my_pe_z;
    this->G->pe2xyz(pct.gridpe, &my_pe_x, &my_pe_y, &my_pe_z);
    int my_pe_offset = my_pe_x % pct.coalesce_factor;


    KpointType *weight = this->nl_weight;
#if HIP_ENABLED || CUDA_ENABLED
    weight = this->nl_weight_gpu;
#endif
    bool potential_acceleration = (ct.potential_acceleration_constant_step > 0.0);

    int pbasis_noncoll = pbasis * ct.noncoll_factor;
    double *nvtot_psi = vtot_psi;;
    if(pct.coalesce_factor > 1)
    {
        nvtot_psi = new double[pbasis * pct.coalesce_factor];
        GatherGrid(this->G, pbasis, vtot_psi, nvtot_psi);
    }

    // Set trade images coalesce_factor
    this->T->set_coalesce_factor(pct.coalesce_factor);

    int active_threads = ct.MG_THREADS_PER_NODE;
    if(ct.mpi_queue_mode && (active_threads > 1)) active_threads--;

    // When grid coalescing is enabled we need nstates to be an integral
    // multiple of (active_threads * pct.coalesce_factor) in order for the
    // coalescing routines to work correctly. So we pad nstates to satisfy
    // that condition where required (stored in mstates).
    int mstates = N / (active_threads * pct.coalesce_factor);
    if(N % (active_threads * pct.coalesce_factor)) mstates++;
    mstates = mstates * (active_threads * pct.coalesce_factor);

    // We adjust the block size here for threading and coalescing
    int block_size = bs;
    block_size = block_size / (active_threads * pct.coalesce_factor);
    block_size = block_size * (active_threads * pct.coalesce_factor);
    int nblocks = mstates / block_size;
    int irem = mstates % block_size;
    if(irem) nblocks++;

    for(int is=first;is < N+first;is++) Kstates[is].eig[1] = Kstates[is].eig[0];

    std::vector<double> deig(20);
    std::fill(deig.begin(), deig.end(), 0.0);
    if(ct.use_rmm_diis)
        for(int is=first;is < N+first;is++) 
            Kstates[is].dptr = new diis<KpointType>(4, pbasis_noncoll);

    for(int vcycle = 0;vcycle < ct.eig_parm.mucycles;vcycle++)
    {

        // Zero out dvh array if potential acceleration is enabled
        if(potential_acceleration)
        {
           int stop = ct.ndvh * pbasis * pct.coalesce_factor;
           for(int i=0;i < stop;i++) this->dvh[i] = 0.0;
           PotentialAccelerationReset(my_pe_offset*active_threads + ct.dvh_skip/pct.coalesce_factor);
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
            for(int st1=0;st1 < block_size;st1+=active_threads*pct.coalesce_factor)
            {
                SCF_THREAD_CONTROL thread_control;

                RT1 = new RmgTimer("3-MgridSubspace: Mg_eig");
                int istart = my_pe_offset*active_threads;
                int nthreads = active_threads;
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
                        if(pct.coalesce_factor > 1 && ct.coalesce_states)
                            thread_control.nv = this->nv - bofs * pbasis_noncoll;
                        else
                            thread_control.nv = (void *)&this->nv[(st1 + ist + istart) * pbasis_noncoll];
                        thread_control.ns = (void *)&this->ns[(first+sindex) * pbasis_noncoll];  // ns is not blocked!
                        thread_control.basetag = this->Kstates[first + sindex].istate;
                        thread_control.extratag1 = active_threads;
                        thread_control.extratag2 = bofs + st1;
                        thread_control.extratag3 = st1 + ist + istart;
                    }
                    QueueThreadTask(ist, thread_control);
                }

                // Thread tasks are set up so run them
                if(!ct.mpi_queue_mode && nthreads) T->run_thread_tasks(nthreads);
                if(ct.mpi_queue_mode) T->run_thread_tasks(active_threads, Rmg_Q);
                delete RT1;
            }
            if(ct.mpi_queue_mode) T->run_thread_tasks(active_threads, Rmg_Q);

        } // end for ib

        RT1 = new RmgTimer("3-MgridSubspace: Mg_eig");
        if(ct.mpi_queue_mode) T->run_thread_tasks(active_threads, Rmg_Q);
        delete RT1;

        for(int is=first;is < N+first;is++)
        {
            Kstates[is].eig[1] = Kstates[is].eig[0];
        }

        // Seems to be necessary for Broyden mixing in some cases.
        if(vcycle != (ct.eig_parm.mucycles-1))
        {
            RmgTimer RTO("3-MgridSubspace: ortho");
            ct.davidson_2stage_ortho=true;
            DavidsonOrtho(first, this->nstates-first,
                          pbasis_noncoll, this->orbital_storage);
        }
    }

    // Set trade images coalesce factor back to 1 for other routines.
    this->T->set_coalesce_factor(1);

    if(pct.coalesce_factor > 1)
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


