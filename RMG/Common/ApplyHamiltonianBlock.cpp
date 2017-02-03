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
#include "State.h"
#include "GlobalSums.h"
#include "BaseThread.h"
#include "TradeImages.h"
#include "RmgTimer.h"
#include "rmgthreads.h"
#include "RmgThread.h"
#include "Subdiag.h"
#include "Solvers.h"

#include "transition.h"


template double ApplyHamiltonianBlock<double>(Kpoint<double> *, int, int, double *, double *);
template double ApplyHamiltonianBlock<std::complex<double> >(Kpoint<std::complex<double>> *, int, int, std::complex<double> *, double *);


// Threaded routine that applies Hamiltonian operator to a block of orbitals of size num_states
// starting from first_state. 

template <typename KpointType>
double ApplyHamiltonianBlock (Kpoint<KpointType> *kptr, int first_state, int num_states, KpointType *h_psi, double *vtot)
{
    int pbasis = kptr->pbasis;
    BaseThread *T = BaseThread::getBaseThread(0);

    int active_threads = ct.MG_THREADS_PER_NODE;
    if(ct.mpi_queue_mode) active_threads--;
    if(active_threads < 1) active_threads = 1;

    int istop = num_states / active_threads;
    istop = istop * active_threads;

    // Apply the non-local operators to this block of orbitals
    AppNls(kptr, kptr->newsint_local, kptr->Kstates[first_state].psi, kptr->nv, &kptr->ns[first_state*pbasis], kptr->Bns,
           first_state, std::min(ct.non_local_block_size, num_states), false);

    int first_nls = 0;

    // Apply Hamiltonian to state 0 to get the diagonal from the finite diff operator. Work is repeated
    // in the thread loop below but that's not much extra work.
    double fd_diag = ApplyHamiltonian (kptr, kptr->Kstates[first_state].psi, &h_psi[first_state*pbasis], vtot, kptr->nv);

    for(int st1=first_state;st1 < first_state + istop;st1+=active_threads) {
        SCF_THREAD_CONTROL thread_control;

        // Make sure the non-local operators are applied for the next block if needed
        int check = first_nls + active_threads;
        if(check > ct.non_local_block_size) {
            AppNls(kptr, kptr->newsint_local, kptr->Kstates[st1].psi, kptr->nv, &kptr->ns[st1 * pbasis], kptr->Bns,
                   st1, std::min(ct.non_local_block_size, num_states + first_state - st1), false);
            first_nls = 0;
        }

        for(int ist = 0;ist < active_threads;ist++) {
            thread_control.job = HYBRID_APPLY_HAMILTONIAN;
            thread_control.vtot = vtot;
            thread_control.istate = st1 + ist;
            thread_control.sp = &kptr->Kstates[st1 + ist];
            thread_control.p1 = (void *)kptr->Kstates[st1 + ist].psi;
            thread_control.p2 = (void *)&h_psi[(st1 + ist) * pbasis];
            thread_control.p3 = (void *)kptr;
            thread_control.nv = (void *)&kptr->nv[(first_nls + ist) * pbasis];
            thread_control.ns = (void *)&kptr->ns[(st1 + ist) * pbasis];  // ns is not blocked!
            thread_control.basetag = kptr->Kstates[st1 + ist].istate;
            QueueThreadTask(ist, thread_control);
        }

        // Thread tasks are set up so run them
        if(!ct.mpi_queue_mode) T->run_thread_tasks(active_threads);

        // Increment index into non-local block
        first_nls += active_threads;
        
    }

    if(ct.mpi_queue_mode) T->run_thread_tasks(active_threads, Rmg_Q);

    // Process any remaining states in serial fashion
    for(int st1 = first_state + istop;st1 < first_state + num_states;st1++) {
         ApplyHamiltonian (kptr, kptr->Kstates[st1].psi, &h_psi[st1 * pbasis], vtot, &kptr->nv[first_nls * pbasis]);
         first_nls++;
    }
    
    return -0.5 * fd_diag;
}
