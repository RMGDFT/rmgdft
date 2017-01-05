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

#include "BaseThread.h"
#include "BaseThreadControl.h"
#include "RmgThread.h"
#include "rmg_error.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "rmgthreads.h"
#include "const.h"
#include "prototypes.h"
#include "transition.h"
#include "Kpoint.h"
#include "Subdiag.h"
#include "Solvers.h"
#ifdef USE_NUMA
    #include "numa.h"
#endif
#if GPU_ENABLED
#include <cuda.h>
#include <cuda_runtime.h>
#endif

#include <sys/types.h>
#include <pthread.h>
#include <sched.h>


// Main thread function specific to subprojects
void *run_threads(void *v) {

    Kpoint<double> *kptr_d;
    Kpoint<std::complex<double>> *kptr_c;
    BaseThreadControl *s;
    SCF_THREAD_CONTROL *ss;
    s = (BaseThreadControl *)v;
    BaseThread *T = BaseThread::getBaseThread(0);

#ifdef USE_NUMA
    struct bitmask *thread_cpumask = NULL;
    if(ct.use_numa) {
        thread_cpumask = numa_allocate_cpumask();

        //numa_set_localalloc();
        // For case with 1 MPI proc per host restrict threads to numa nodes
        if(pct.procs_per_host == 1) 
        {
            int t_tid = 0;
            for(unsigned int idx=0;idx<pct.cpumask->size;idx++)
            {
                if(numa_bitmask_isbitset(pct.cpumask, idx))
                {
                    if(s->tid == t_tid)
                    {
                        numa_bitmask_clearall(thread_cpumask);
                        numa_bitmask_setbit(thread_cpumask, idx);
                        numa_sched_setaffinity(0, thread_cpumask);
                        if(ct.verbose) printf("C1 Binding rank %d thread %d to cpu %d.\n", pct.local_rank, s->tid, idx);
                        break;
                    }
                    t_tid++;
                }
            }
        }
        else if(pct.ncpus == pct.procs_per_host) 
        {
            copy_bitmask_to_bitmask(pct.cpumask, thread_cpumask);
        }
        else if(pct.procs_per_host == pct.numa_nodes_per_host)
        {
            int t_tid = 0;
            for(unsigned int idx=0;idx<pct.cpumask->size;idx++)
            {
                if(numa_bitmask_isbitset(pct.cpumask, idx))
                {
                    if(s->tid == t_tid)
                    {
                        numa_bitmask_clearall(thread_cpumask);
                        numa_bitmask_setbit(thread_cpumask, idx);
                        numa_sched_setaffinity(0, thread_cpumask);
                        if(ct.verbose) printf("C2 Binding rank %d thread %d to cpu %d.\n", pct.local_rank, s->tid, idx);
                        break;
                    }
                    t_tid++;
                }
            }
        }
        else if((pct.procs_per_host > pct.numa_nodes_per_host) && (pct.ncpus != pct.procs_per_host))
        {
            // Multiple MPI procs per numa node so default is to bind threads to a specific node
            copy_bitmask_to_bitmask(pct.cpumask, thread_cpumask);

            // If threads*procs = cpus/node then do a one to one binding
            unsigned int cpus_per_node = numa_bitmask_weight(pct.cpumask);
            unsigned int procs_per_node = pct.procs_per_host / pct.numa_nodes_per_host;
            if(procs_per_node*ct.THREADS_PER_NODE == cpus_per_node && !(pct.procs_per_host % pct.numa_nodes_per_host))
            {
                unsigned int proc_rank_in_node = pct.local_rank % procs_per_node;
                unsigned int nid = pct.local_rank / (pct.procs_per_host / pct.numa_nodes_per_host);
                unsigned int t_tid = nid*cpus_per_node;
                for(unsigned int idx=0;idx<pct.cpumask->size;idx++)
                {
                    if(numa_bitmask_isbitset(pct.cpumask, idx))
                    {
                        if((s->tid + proc_rank_in_node * ct.THREADS_PER_NODE + nid*cpus_per_node) == t_tid)
                        {
                            numa_bitmask_clearall(thread_cpumask);
                            numa_bitmask_setbit(thread_cpumask, idx);
                            numa_sched_setaffinity(0, thread_cpumask);
                            if(ct.verbose) printf("C3 Binding rank %d thread %d to cpu %d.\n", pct.local_rank, s->tid, idx);
                            break;
                        }
                        t_tid++;
                    }
                }
            }
        }
        else
        {
            ct.use_numa = false;
        }

        numa_sched_setaffinity(0, thread_cpumask);
        numa_set_localalloc();

    }


#endif
    

#if GPU_ENABLED
    cudaError_t cuerr;
#endif

    // Set up thread local storage
    rmg_set_tsd(s);


    // Get the control structure
    ss = (SCF_THREAD_CONTROL *)s->pptr;

#if GPU_ENABLED
    cudaSetDevice(ct.cu_dev); 
#endif

    while(1) {

#ifdef USE_NUMA
        if(ct.use_numa) numa_sched_setaffinity(0, thread_cpumask);
#endif

        // We sleep until main wakes us up
        T->thread_sleep();

#ifdef USE_NUMA
        if(ct.use_numa) numa_sched_setaffinity(0, thread_cpumask);
        //T->set_cpu_affinity(s->tid, pct.procs_per_host, pct.local_rank);
#endif

        // Get the control structure
        ss = (SCF_THREAD_CONTROL *)s->pptr;

        // Switch that controls what we do
        switch(ss->job) {
            case HYBRID_EIG:       // Performs a single multigrid sweep over an orbital
               if(ct.is_gamma) {
                   kptr_d = (Kpoint<double> *)ss->p3;
                   if(ct.rms > ct.preconditioner_thr)
                       MgEigState<double,float> (kptr_d, (State<double> *)ss->sp, ss->vtot, (double *)ss->nv, (double *)ss->ns, ss->vcycle);
                   else
                       MgEigState<double,double> (kptr_d, (State<double> *)ss->sp, ss->vtot, (double *)ss->nv, (double *)ss->ns, ss->vcycle);
               }
               else {
                   kptr_c = (Kpoint<std::complex<double>> *)ss->p3;
                   if(ct.rms > ct.preconditioner_thr)
                       MgEigState<std::complex<double>, std::complex<float> > (kptr_c, (State<std::complex<double> > *)ss->sp, ss->vtot, (std::complex<double> *)ss->nv, (std::complex<double> *)ss->ns, ss->vcycle);
                   else
                       MgEigState<std::complex<double>, std::complex<double> > (kptr_c, (State<std::complex<double> > *)ss->sp, ss->vtot, (std::complex<double> *)ss->nv, (std::complex<double> *)ss->ns, ss->vcycle);
               }
               break;
            case HYBRID_SKIP:
               break;
            case HYBRID_SUBDIAG_APP_AB:
               if(ct.is_gamma) {
                   kptr_d = (Kpoint<double> *)ss->p3;
                   State<double> *spd = (State<double> *)ss->sp;
                   ApplyOperators<double> (kptr_d, spd->istate, (double *)ss->p1, (double *)ss->p2, ss->vtot, 
                                          (double *)ss->nv, (double *)ss->Bns);
               }
               else {
                   kptr_c = (Kpoint<std::complex<double>> *)ss->p3;
                   State<std::complex<double> > *spc = (State<std::complex<double> > *)ss->sp;
                   ApplyOperators<std::complex<double> > (kptr_c, spc->istate, (std::complex<double> *)ss->p1, (std::complex<double> *)ss->p2, ss->vtot,
                                                          (std::complex<double> *)ss->nv, (std::complex<double> *)ss->Bns);
               } 
               break;
            case HYBRID_APPLY_HAMILTONIAN:
               if(ct.is_gamma) {
                   kptr_d = (Kpoint<double> *)ss->p3;
                   ApplyHamiltonian<double> (kptr_d, (double *)ss->p1, (double *)ss->p2, ss->vtot, (double *)ss->nv);
               }
               else {
                   kptr_c = (Kpoint<std::complex<double>> *)ss->p3;
                   ApplyHamiltonian<std::complex<double> > (kptr_c, (std::complex<double> *)ss->p1, (std::complex<double> *)ss->p2, ss->vtot, 
                                          (std::complex<double> *)ss->nv);
               } 
               break;
            case HYBRID_DAV_PRECONDITIONER:
               if(ct.is_gamma) {
                   kptr_d = (Kpoint<double> *)ss->p1;
                   DavPreconditionerOne<double> (kptr_d, (double *)ss->p2, ss->fd_diag, ss->eig, ss->vtot, ss->avg_potential);
               }
               else {
               } 
               break;
            case HYBRID_BETAX_PSI1_CALCULATE:
//               betaxpsi1_calculate_one(ss->sp, ss->ion, ss->nion, ss->sintR, ss->sintI, ss->kpt, ss->weiptr);
               break;
            default:
               break;
        }

    }

    return NULL;
}


