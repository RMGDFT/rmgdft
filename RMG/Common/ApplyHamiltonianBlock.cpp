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


template void ApplyHamiltonianBlock<double>(Kpoint<double> *, int, int, double *, double *);
template void ApplyHamiltonianBlock<std::complex<double> >(Kpoint<std::complex<double>> *, int, int, std::complex<double> *, double *);


// Threaded routine that applies Hamiltonian operator to a block of orbitals of size num_states
// starting from first_state. 

template <typename KpointType>
void ApplyHamiltonianBlock (Kpoint<KpointType> *kptr, int first_state, int num_states, KpointType *h_psi, double *vtot)
{
    int pbasis = kptr->pbasis;
    BaseThread *T = BaseThread::getBaseThread(0);

    int istop = num_states / T->get_threads_per_node();
    istop = istop * T->get_threads_per_node();
 
    // Apply the non-local operators to this block of orbitals
    AppNls(kptr, kptr->oldsint_local, kptr->Kstates[first_state].psi, kptr->nv, kptr->ns, kptr->Bns,
           first_state, std::min(ct.non_local_block_size, num_states));
    int first_nls = 0;

    for(int st1=first_state;st1 < first_state + istop;st1+=T->get_threads_per_node()) {
        SCF_THREAD_CONTROL thread_control[MAX_RMG_THREADS];

        // Make sure the non-local operators are applied for the next block if needed
        int check = first_nls + T->get_threads_per_node();
        if(check > ct.non_local_block_size) {
            AppNls(kptr, kptr->oldsint_local, kptr->Kstates[st1].psi, kptr->nv, &kptr->ns[st1 * pbasis], kptr->Bns,
                   st1, std::min(ct.non_local_block_size, num_states - st1));
            first_nls = 0;
        }

        for(int ist = 0;ist < T->get_threads_per_node();ist++) {
            thread_control[ist].job = HYBRID_APPLY_HAMILTONIAN;
            thread_control[ist].vtot = vtot;
            thread_control[ist].istate = st1 + ist;
            thread_control[ist].sp = &kptr->Kstates[st1 + ist];
            thread_control[ist].p1 = (void *)&kptr->Kstates[st1 + ist].psi;
            thread_control[ist].p2 = (void *)&h_psi[(st1 + ist) * pbasis];
            thread_control[ist].p3 = (void *)kptr;
            thread_control[ist].nv = (void *)&kptr->nv[(first_nls + ist) * pbasis];
            thread_control[ist].ns = (void *)&kptr->ns[(st1 + ist) * pbasis];  // ns is not blocked!
            T->set_pptr(ist, &thread_control[ist]);
        }

        // Thread tasks are set up so run them
        T->run_thread_tasks(T->get_threads_per_node());

        // Increment index into non-local block
        first_nls += T->get_threads_per_node();

        
    }

    // Process any remaining states in serial fashion
    for(int st1 = first_state + istop;st1 < num_states;st1++) {
         ApplyHamiltonian (kptr, kptr->Kstates[st1].psi, &h_psi[st1 * pbasis], vtot, &kptr->nv[first_nls * pbasis]);
         first_nls++;
    }
    
}
