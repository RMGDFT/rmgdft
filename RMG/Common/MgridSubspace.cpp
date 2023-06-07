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
#include "blas.h"


// Solver that uses multigrid preconditioning and subspace rotations
// Part of Kpoint class

template void Kpoint<double>::MgridSubspace(double *, double *vxc_psi);
template void Kpoint<std::complex<double>>::MgridSubspace(double *, double *vxc_psi);


template <class KpointType> void Kpoint<KpointType>::MgridSubspace (double *vtot_psi, double *vxc_psi)
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
    double mean_occ_res = DBL_MAX;
    double mean_unocc_res = DBL_MAX;
    double max_occ_res = 0.0;
    double max_unocc_res = 0.0;
    double min_occ_res = DBL_MAX;
    double min_unocc_res = DBL_MAX;
    bool potential_acceleration = (ct.potential_acceleration_constant_step > 0.0);

    int pbasis_noncoll = pbasis * ct.noncoll_factor;
    double *nvtot_psi = vtot_psi;;
    double *coarse_vtot;
    if(pct.coalesce_factor > 1)
    {
        nvtot_psi = new double[pbasis * pct.coalesce_factor];
        coarse_vtot = new double[pbasis * pct.coalesce_factor];
        GatherGrid(this->G, pbasis, vtot_psi, nvtot_psi);
    }
    else
    {
        coarse_vtot = new double[pbasis];
    }

    // Set trade images coalesce_factor
    this->T->set_coalesce_factor(pct.coalesce_factor);

    // Restrict coarse vtot
    {
        Mgrid MG(this->L, this->T);
        int ixoff, iyoff, izoff;
        int dimx = pct.coalesce_factor * G->get_PX0_GRID(1);
        int dimy = G->get_PY0_GRID(1);
        int dimz = G->get_PZ0_GRID(1);
        int NX_GRID = G->get_NX_GRID(1);
        int NY_GRID = G->get_NY_GRID(1);
        int NZ_GRID = G->get_NZ_GRID(1);
        int sbasis = (dimx + 2)*(dimy + 2)*(dimz + 2);
        int dx2 = MG.MG_SIZE (dimx, 0, NX_GRID, G->get_PX_OFFSET(1), dimx, &ixoff, ct.boundaryflag);
        int dy2 = MG.MG_SIZE (dimy, 0, NY_GRID, G->get_PY_OFFSET(1), dimy, &iyoff, ct.boundaryflag);
        int dz2 = MG.MG_SIZE (dimz, 0, NZ_GRID, G->get_PZ_OFFSET(1), dimz, &izoff, ct.boundaryflag);
        double *work = new double[sbasis]();
        CPP_pack_ptos_convert (work, nvtot_psi, dimx, dimy, dimz);
        this->T->trade_images (work, dimx, dimy, dimz, FULL_TRADE);
        MG.mg_restrict (work, coarse_vtot, dimx, dimy, dimz, dx2, dy2, dz2, ixoff, iyoff, izoff);
        delete [] work;
    }
 
    int active_threads = ct.MG_THREADS_PER_NODE;
    if(ct.mpi_queue_mode && (active_threads > 1)) active_threads--;

    // When grid coalescing is enabled we need nstates to be an integral
    // multiple of (active_threads * pct.coalesce_factor) in order for the
    // coalescing routines to work correctly. So we pad nstates to satisfy
    // that condition where required (stored in mstates).
    int mstates = this->nstates / (active_threads * pct.coalesce_factor);
    if(this->nstates % (active_threads * pct.coalesce_factor)) mstates++;
    mstates = mstates * (active_threads * pct.coalesce_factor);

    // We adjust the block size here for threading and coalescing
    int block_size = ct.non_local_block_size;
    block_size = block_size / (active_threads * pct.coalesce_factor);
    block_size = block_size * (active_threads * pct.coalesce_factor);
    int nblocks = mstates / block_size;
    int irem = mstates % block_size;
    if(irem) nblocks++;

    for(int vcycle = 0;vcycle < ct.eig_parm.mucycles;vcycle++)
    {

        // Zero out dvh array if potential acceleration is enabled
        if(potential_acceleration)
        {
           int stop = this->ndvh * pbasis * pct.coalesce_factor;
           for(int i=0;i < stop;i++) this->dvh[i] = 0.0;
           PotentialAccelerationReset(my_pe_offset*active_threads + this->dvh_skip/pct.coalesce_factor);
        }

        // Update betaxpsi        
        RT1 = new RmgTimer("3-MgridSubspace: Beta x psi");
        this->BetaProjector->project(this, this->newsint_local, 0, mstates * ct.noncoll_factor, weight);
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
            AppNls(this, this->newsint_local, this->Kstates[bofs].psi, this->nv, 
                   &this->ns[bofs * pbasis_noncoll],
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
                        thread_control.job = HYBRID_SKIP;
                        if(!ct.mpi_queue_mode && nthreads == active_threads) nthreads = ist;
                    }
                    else
                    {
                        thread_control.job = HYBRID_EIG;
                        thread_control.vtot = nvtot_psi;
                        thread_control.vxc_psi = vxc_psi;
                        thread_control.coarse_vtot = coarse_vtot;
                        thread_control.vcycle = vcycle;
                        thread_control.sp = &this->Kstates[sindex];
                        thread_control.p3 = (void *)this;
// this is a hack since we are backing the nv array up before it's start in order to make
// the code in MgEigState work properly. Fix at some point.
                        if(pct.coalesce_factor > 1 && ct.coalesce_states)
                            thread_control.nv = this->nv - bofs * pbasis_noncoll;
                        else
                            thread_control.nv = (void *)&this->nv[(st1 + ist + istart) * pbasis_noncoll];
                        thread_control.ns = (void *)&this->ns[sindex * pbasis_noncoll];  // ns is not blocked!
                        thread_control.basetag = this->Kstates[sindex].istate;
                        thread_control.extratag1 = active_threads;
                        thread_control.extratag2 = bofs + st1;
                        thread_control.extratag3 = st1 + ist + istart;
                    }
                    QueueThreadTask(ist, thread_control);
                }

                // Thread tasks are set up so run them
                if(!ct.mpi_queue_mode && nthreads) T->run_thread_tasks(nthreads);
                delete RT1;
            }
            if(ct.mpi_queue_mode) T->run_thread_tasks(active_threads, Rmg_Q);

        } // end for ib

        RT1 = new RmgTimer("3-MgridSubspace: Mg_eig");
        if(ct.mpi_queue_mode) T->run_thread_tasks(active_threads, Rmg_Q);
        delete RT1;
    }

    // Set trade images coalesce factor back to 1 for other routines.
    this->T->set_coalesce_factor(1);

    if(pct.coalesce_factor > 1)
    {
        delete [] nvtot_psi;
        // Eigenvalues are not copied to all nodes in MgEigState when using coalesced grids.
        GatherEigs(this);
    }
    delete [] coarse_vtot;

    if(Verify ("freeze_occupied", true, this->ControlMap)) {

        // Orbital residual measures (used for some types of calculations
        this->max_unocc_res_index = (int)(ct.gw_residual_fraction * (double)this->nstates);
        this->mean_occ_res = 0.0;
        this->min_occ_res = DBL_MAX;
        this->max_occ_res = 0.0;
        this->mean_unocc_res = 0.0;
        this->min_unocc_res = DBL_MAX;
        this->max_unocc_res = 0.0;
        this->highest_occupied = 0;
        for(int istate = 0;istate < this->nstates;istate++) {
            if(this->Kstates[istate].occupation[0] > 0.0) {
                this->mean_occ_res += this->Kstates[istate].res;
                mean_occ_res += this->Kstates[istate].res;
                if(this->Kstates[istate].res >  this->max_occ_res)  this->max_occ_res = this->Kstates[istate].res;
                if(this->Kstates[istate].res <  this->min_occ_res)  this->min_occ_res = this->Kstates[istate].res;
                if(this->Kstates[istate].res >  max_occ_res)  max_occ_res = this->Kstates[istate].res;
                if(this->Kstates[istate].res <  min_occ_res)  min_occ_res = this->Kstates[istate].res;
                this->highest_occupied = istate;
            }
            else {
                if(istate <= this->max_unocc_res_index) {
                    this->mean_unocc_res += this->Kstates[istate].res;
                    mean_unocc_res += this->Kstates[istate].res;
                    if(this->Kstates[istate].res >  this->max_unocc_res)  this->max_unocc_res = this->Kstates[istate].res;
                    if(this->Kstates[istate].res <  this->min_unocc_res)  this->min_unocc_res = this->Kstates[istate].res;
                    if(this->Kstates[istate].res >  max_unocc_res)  max_unocc_res = this->Kstates[istate].res;
                    if(this->Kstates[istate].res <  min_unocc_res)  min_unocc_res = this->Kstates[istate].res;
                }
            }
        }
        this->mean_occ_res = this->mean_occ_res / (double)(this->highest_occupied + 1);
        this->mean_unocc_res = this->mean_unocc_res / (double)(this->max_unocc_res_index -(this->highest_occupied + 1));
        mean_occ_res = mean_occ_res / (double)(ct.num_kpts*(this->highest_occupied + 1));
        mean_unocc_res = mean_unocc_res / (double)(ct.num_kpts*this->max_unocc_res_index -(this->highest_occupied + 1));

        rmg_printf("Mean/Min/Max unoccupied wavefunction residual for kpoint %d  =  %10.5e  %10.5e  %10.5e\n", this->kidx, this->mean_unocc_res, this->min_unocc_res, this->max_unocc_res);

    }


    /* wavefunctions have changed, projectors have to be recalculated
     * but if we are using potential acceleration and not well converged yet
     * it is counterproductive to do so */
    if(!potential_acceleration || (potential_acceleration && (ct.rms <  5.0e-6)))
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


    RT1 = new RmgTimer("3-MgridSubspace: Diagonalization");
    this->Subdiag (vtot_psi, vxc_psi, ct.subdiag_driver);
    delete(RT1);

    // wavefunctions have changed, projectors have to be recalculated */
    RT1 = new RmgTimer("3-MgridSubspace: Beta x psi");
    this->BetaProjector->project(this, this->newsint_local, 0, nstates * ct.noncoll_factor, weight);
    delete(RT1);

}


