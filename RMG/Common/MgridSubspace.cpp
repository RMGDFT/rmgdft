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


// Solver that uses multigrid preconditioning and subspace rotations
// Part of Kpoint class

template void Kpoint<double>::MgridSubspace(double *);
template void Kpoint<std::complex<double>>::MgridSubspace(double *);


template <class KpointType> void Kpoint<KpointType>::MgridSubspace (double *vtot_psi)
{
    RmgTimer RT0("3-MgridSubspace"), *RT1;
    BaseThread *T = BaseThread::getBaseThread(0);
    int my_pe_x, my_pe_y, my_pe_z;
    this->G->pe2xyz(pct.gridpe, &my_pe_x, &my_pe_y, &my_pe_z);
    int my_pe_offset = my_pe_x % pct.coalesce_factor;

    double mean_occ_res = DBL_MAX;
    double mean_unocc_res = DBL_MAX;
    double max_occ_res = 0.0;
    double max_unocc_res = 0.0;
    double min_occ_res = DBL_MAX;
    double min_unocc_res = DBL_MAX;
    bool potential_acceleration = (ct.potential_acceleration_constant_step > 0.0);


    double *nvtot_psi = vtot_psi;;
    if(pct.coalesce_factor > 1)
    {
        nvtot_psi = new double[pbasis * pct.coalesce_factor];
        GatherGrid(this->G, pbasis, vtot_psi, nvtot_psi);
    }
    
    // Set trade images coalesce_factor
    this->T->set_coalesce_factor(pct.coalesce_factor);
    for(int vcycle = 0;vcycle < ct.eig_parm.mucycles;vcycle++)
    {

        int active_threads = ct.MG_THREADS_PER_NODE;
        if(ct.mpi_queue_mode && (active_threads > 1)) active_threads--;

        // Zero out dvh array if potential acceleration is enabled
        if(potential_acceleration)
        {
           int stop = this->ndvh * this->pbasis * pct.coalesce_factor;
           for(int i=0;i < stop;i++) this->dvh[i] = 0.0;
           PotentialAccelerationReset(my_pe_offset*active_threads + this->dvh_skip/pct.coalesce_factor);
        }

        // Update betaxpsi        
        RT1 = new RmgTimer("3-MgridSubspace: Beta x psi");
        this->BetaProjector->project(this, this->newsint_local, 0, nstates, this->nl_weight);
        delete(RT1);

        if(ct.ldaU_mode != LDA_PLUS_U_NONE)
        {
            RmgTimer RTL("3-MgridSubspace: ldaUop x psi");
            LdaplusUxpsi(this, 0, this->nstates, this->orbitalsint_local);
            this->ldaU->calc_ns_occ(this->orbitalsint_local, 0, this->nstates);
        }


        /* Update the wavefunctions */
        int istop = this->nstates / (active_threads * pct.coalesce_factor);
        istop = istop * active_threads * pct.coalesce_factor;

        // Apply the non-local operators to a block of orbitals
        RT1 = new RmgTimer("3-MgridSubspace: AppNls");
        AppNls(this, this->newsint_local, this->Kstates[0].psi, this->nv, this->ns, this->Bns,
               0, std::min(ct.non_local_block_size, this->nstates));
        delete(RT1);
        int first_nls = 0;

        int st1 = 0;
        while(st1 < this->nstates)
        {

            // Adjust thread count in case num_states is not evenly divisible by the number of threads
            while(active_threads > 1)
            {
                int icheck = st1 + active_threads*pct.coalesce_factor;
                if(icheck > this->nstates) 
                {
                    active_threads--;
                }
                else
                {
                    break;
                }
            }

            SCF_THREAD_CONTROL thread_control;

            // Make sure the non-local operators are applied for the next block if needed
            int check = first_nls + active_threads*pct.coalesce_factor;
            if(check > ct.non_local_block_size) 
            {
                RT1 = new RmgTimer("3-MgridSubspace: AppNls");
                AppNls(this, this->newsint_local, this->Kstates[st1].psi, this->nv, &this->ns[st1 * pbasis], this->Bns,
                       st1, std::min(ct.non_local_block_size, this->nstates - st1));
                first_nls = 0;
                delete(RT1);
            }
        
            RT1 = new RmgTimer("3-MgridSubspace: Mg_eig");
            int istart = my_pe_offset*active_threads;
            for(int ist = 0;ist < active_threads;ist++) {
                if((st1 + ist + istart) >= this->nstates) break;
                thread_control.job = HYBRID_EIG;
                thread_control.vtot = nvtot_psi;
                thread_control.vcycle = vcycle;
                thread_control.sp = &this->Kstates[st1 + ist + istart];
                thread_control.p3 = (void *)this;
                thread_control.nv = (void *)&this->nv[(first_nls + ist + istart) * pbasis];
                thread_control.ns = (void *)&this->ns[(st1 + ist + istart) * pbasis];  // ns is not blocked!
                thread_control.basetag = this->Kstates[st1 + ist + istart].istate;
                thread_control.extratag1 = active_threads;
                thread_control.extratag2 = st1;
                QueueThreadTask(ist, thread_control);
            }

            // Thread tasks are set up so run them
            if(!ct.mpi_queue_mode) T->run_thread_tasks(active_threads);
            if((check >= ct.non_local_block_size) && ct.mpi_queue_mode) T->run_thread_tasks(active_threads, Rmg_Q);

            delete RT1;

            // Increment index into non-local block and state index
            first_nls += active_threads*pct.coalesce_factor;
            st1+=active_threads*pct.coalesce_factor;
        }

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
        this->BetaProjector->project(this, this->newsint_local, 0, nstates, this->nl_weight);
        delete(RT1);

        if(ct.ldaU_mode != LDA_PLUS_U_NONE)
        {   
            RmgTimer RTL("3-MgridSubspace: ldaUop x psi"); 
            LdaplusUxpsi(this, 0, this->nstates, this->orbitalsint_local);
            this->ldaU->calc_ns_occ(this->orbitalsint_local, 0, this->nstates);
        }
    }


    RT1 = new RmgTimer("3-MgridSubspace: Diagonalization");
    this->Subdiag (vtot_psi, ct.subdiag_driver);
    delete(RT1);

    // wavefunctions have changed, projectors have to be recalculated */
    RT1 = new RmgTimer("3-MgridSubspace: Beta x psi");
    this->BetaProjector->project(this, this->newsint_local, 0, nstates, this->nl_weight);
    delete(RT1);

    if(ct.ldaU_mode != LDA_PLUS_U_NONE)
    {   
        RmgTimer RTL("3-MgridSubspace: ldaUop x psi"); 
        LdaplusUxpsi(this, 0, this->nstates, this->orbitalsint_local);
        this->ldaU->calc_ns_occ(this->orbitalsint_local, 0, this->nstates);
    }

}


