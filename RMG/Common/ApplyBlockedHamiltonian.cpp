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

template <typename KpointType>
void ApplyBlockedHamiltonian (Kpoint<KpointType> *kptr, int first, int N, int bs, KpointType *h_psi, double *vtot_psi, double *vxc_psi);


template <typename KpointType>
void ApplyBlockedHamiltonian(Kpoint<KpointType> *kptr, KpointType *h_psi, double *vtot_psi, double *vxc_psi)
{
    int bs = ct.non_local_block_size;
    int nb = kptr->nstates / bs;
    int irem = kptr->nstates % bs;
    for(int ib=0;ib < nb;ib++)
    {
        ApplyBlockedHamiltonian (kptr, ib*bs, bs, bs, h_psi, vtot_psi, vxc_psi);
    }
    if(irem)
    {
        ApplyBlockedHamiltonian (kptr, nb*bs, irem, bs,  h_psi,vtot_psi, vxc_psi);
    }


}

template <typename KpointType>
void ApplyBlockedHamiltonian (Kpoint<KpointType> *kptr, int first, int N, int bs, KpointType *h_psi, double *vtot_psi, double *vxc_psi)
{

    BaseThread *T = BaseThread::getBaseThread(0);
    int cfac = 1;
    if(ct.coalesce_states) cfac = pct.coalesce_factor;
    int my_pe_x, my_pe_y, my_pe_z;
    kptr->G->pe2xyz(pct.gridpe, &my_pe_x, &my_pe_y, &my_pe_z);
    int my_pe_offset = my_pe_x % cfac;

    int pbasis_noncoll = kptr->pbasis * ct.noncoll_factor;

    double *nvtot_psi = vtot_psi;;
    if(cfac > 1)
    {
        nvtot_psi = new double[kptr->pbasis * cfac];
        GatherGrid(kptr->G, kptr->pbasis, vtot_psi, nvtot_psi);
    }

    // Set trade images coalesce_factor
    kptr->T->set_coalesce_factor(cfac);

    int active_threads = ct.MG_THREADS_PER_NODE;
    if(ct.mpi_queue_mode && (active_threads > 1)) active_threads--;

    // When grid coalescing is enabled we need nstates to be an integral
    // multiple of (active_threads * pct.coalesce_factor) in order for the
    // coalescing routines to work correctly. So we pad nstates to satisfy
    // that condition where required (stored in mstates).
    int mstates = N / (active_threads * cfac);
    if(N % (active_threads * cfac)) mstates++;
    mstates = mstates * (active_threads * cfac);

    // We adjust the block size here for threading and coalescing
    int block_size = bs;
    block_size = block_size / (active_threads * cfac);
    block_size = block_size * (active_threads * cfac);
    int nblocks = mstates / block_size;
    int irem = mstates % block_size;
    if(irem) nblocks++;

        for(int ib = 0;ib < nblocks;ib++)
        {
            int bofs = ib * block_size;
            KpointType *newsint = kptr->newsint_local +
                    first * kptr->BetaProjector->get_num_nonloc_ions() *
                    ct.max_nl * ct.noncoll_factor;

            AppNls(kptr, newsint, kptr->Kstates[first + bofs].psi, kptr->nv, 
                   &kptr->ns[(first+bofs) * pbasis_noncoll],
                   bofs, std::min(block_size, mstates - bofs));
            for(int st1=0;st1 < block_size;st1+=active_threads*cfac)
            {
                SCF_THREAD_CONTROL thread_control;

                int istart = my_pe_offset*active_threads;
                for(int ist = 0;ist < active_threads;ist++) {
                    int sindex = bofs + st1 + ist + istart;
                    if(sindex >= mstates)
                    {
                        break;
                    }
                    else
                    {
                        thread_control.job = SUBDIAG_HAMILTONIAN;
                        thread_control.vtot = nvtot_psi;
                        thread_control.vxc_psi = vxc_psi;
                        thread_control.coarse_vtot = NULL;
                        thread_control.sp = &kptr->Kstates[first + sindex];
                        thread_control.p3 = (void *)kptr;
                        if(cfac > 1)
                            thread_control.nv = kptr->nv - (first+bofs) * pbasis_noncoll;
                        else
                            thread_control.nv = (void *)&kptr->nv[(st1 + ist + istart) * pbasis_noncoll];
                        thread_control.ns = (void *)&kptr->ns[(first+sindex) * pbasis_noncoll];  // ns is not blocked!
                        thread_control.basetag = kptr->Kstates[first + sindex].istate;
                        thread_control.extratag1 = active_threads;
                        thread_control.extratag2 = bofs + st1;
                        thread_control.extratag2 = first + bofs + st1;
                        thread_control.extratag3 = st1 + ist + istart;

                    }
                    QueueThreadTask(ist, thread_control);
                }

                // Thread tasks are set up so run them
                if(!ct.mpi_queue_mode && active_threads) T->run_thread_tasks(active_threads);
                if(ct.mpi_queue_mode) T->run_thread_tasks(active_threads, Rmg_Q);
            }
            if(!ct.mpi_queue_mode && active_threads) T->run_thread_tasks(active_threads);
            if(ct.mpi_queue_mode) T->run_thread_tasks(active_threads, Rmg_Q);

        } // end for ib

        if(!ct.mpi_queue_mode && active_threads) T->run_thread_tasks(active_threads);
        if(ct.mpi_queue_mode) T->run_thread_tasks(active_threads, Rmg_Q);


    // Set trade images coalesce factor back to 1 for other routines.
    kptr->T->set_coalesce_factor(1);

    if(cfac > 1) delete [] nvtot_psi;

}


template void ApplyBlockedHamiltonian (Kpoint<double> *kptr, int first, int N, int bs, double *h_psi, double *vtot_psi, double *vxc_psi);
template void ApplyBlockedHamiltonian(Kpoint<double> *kptr, double *h_psi, double *vtot_psi, double *vxc_psi);

template void ApplyBlockedHamiltonian (Kpoint<std::complex<double>> *kptr, int first, int N, int bs, std::complex<double> *h_psi, double *vtot_psi, double *vxc_psi);
template void ApplyBlockedHamiltonian(Kpoint<std::complex<double>> *kptr, std::complex<double> *h_psi, double *vtot_psi, double *vxc_psi);
