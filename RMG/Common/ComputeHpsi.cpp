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
#include "rmg_error.h"
#include "rmgthreads.h"
#include "RmgTimer.h"
#include "RmgThread.h"
#include "GlobalSums.h"
#include "Kpoint.h"
#include "RmgGemm.h"
#include "Gpufuncs.h"
#include "Subdiag.h"
#include "GpuAlloc.h"
#include "ErrorFuncs.h"
#include "blas.h"
#include "Solvers.h"
#include "Functional.h"
#include "RmgMatrix.h"

#include "common_prototypes.h"
#include "common_prototypes1.h"
#include "transition.h"

#if CUDA_ENABLED
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cublas_v2.h>
#endif



template void Kpoint<double>::ComputeHpsi(double *, double *, double *);
template void Kpoint<std::complex<double>>::ComputeHpsi(double *, double *, 
        std::complex<double> * );

template <class KpointType> void Kpoint<KpointType>::ComputeHpsi (double *vtot_eig, double *vxc_psi, 
        KpointType *h_psi)
{
    RmgTimer RT0("Compute Hpsi");

    bool potential_acceleration = (ct.potential_acceleration_constant_step > 0.0) && (ct.scf_steps > 0);

    int pbasis_noncoll = pbasis * ct.noncoll_factor;
    BaseThread *T = BaseThread::getBaseThread(0);
    // State array is 4 * the number of states in length but memory above
    // the first set of nstates is unused in this routine so we can use it
    // as temporary space.
    KpointType *tmp_arrayT = Kstates[0].psi;
    tmp_arrayT += nstates * pbasis_noncoll ;

    // Apply operators on each wavefunction
    RmgTimer *RT1 = new RmgTimer("Compute Hpsi: apply operators");

    int active_threads = ct.MG_THREADS_PER_NODE;
    if(ct.mpi_queue_mode) active_threads--;
    if(active_threads < 1) active_threads = 1;

    // We adjust the block size here for threading
    int block_size = ct.non_local_block_size;
    block_size = block_size / active_threads;
    block_size = block_size * active_threads;
    int nblocks = this->nstates / block_size;
    int irem = this->nstates % block_size;
    if(irem) nblocks++;

    for(int ib = 0;ib < nblocks;ib++)
    {
        int bofs = ib * block_size;
        RmgTimer *RT3 = new RmgTimer("4-Diagonalization: AppNls");
        AppNls(this, this->newsint_local, this->Kstates[bofs].psi, this->nv, 
               &this->ns[bofs * pbasis_noncoll],
               bofs, std::min(block_size, this->nstates - bofs));
        delete(RT3);
        for(int st1=0;st1 < block_size;st1+=active_threads)
        {
            SCF_THREAD_CONTROL thread_control;

            int nthreads = active_threads;
            for(int ist = 0;ist < active_threads;ist++) {
                int sindex = bofs + st1 + ist;
                if(sindex >= this->nstates)
                {
                    thread_control.job = HYBRID_SKIP;
                    if(!ct.mpi_queue_mode && nthreads == active_threads) nthreads = ist;
                }
                else
                {
                    thread_control.job = HYBRID_APPLY_HAMILTONIAN;
                    thread_control.vtot = vtot_eig;
                    thread_control.vxc_psi = vxc_psi;
                    thread_control.extratag1 = potential_acceleration;
                    thread_control.istate = sindex;
                    thread_control.sp = &this->Kstates[sindex];
                    thread_control.p1 = (void *)Kstates[sindex].psi;
                    thread_control.p2 = (void *)&h_psi[sindex * pbasis_noncoll];
                    thread_control.p3 = (void *)this;
                    thread_control.nv = (void *)&this->nv[(st1 + ist) * pbasis_noncoll];
                    thread_control.ns = (void *)&this->ns[sindex * pbasis_noncoll];  // ns is not blocked!
                    thread_control.basetag = this->Kstates[sindex].istate;

                }
                QueueThreadTask(ist, thread_control);
            }

            // Thread tasks are set up so run them
            if(!ct.mpi_queue_mode && nthreads) T->run_thread_tasks(nthreads);
        }
        if(ct.mpi_queue_mode) T->run_thread_tasks(active_threads, Rmg_Q);

    } // end for ib

    delete(RT1);

    DeviceSynchronize();

}

