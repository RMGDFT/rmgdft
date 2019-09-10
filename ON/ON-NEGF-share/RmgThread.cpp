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
#include "transition.h"
#include <boost/lockfree/queue.hpp>
#include <boost/lockfree/spsc_queue.hpp>

#include "prototypes_on.h"
#include "init_var.h"

extern std::vector<ORBITAL_PAIR> OrbitalPairs;


// Work queue(s). Not sure if thread specific queues or a global queue is best here. With a global queue we can ensure
// that if one thread is running much slower than the others we will still get the max throughput. Per thread queues
// will reduce contention on the queue though. Also have to consider how big to make the queue. If it becomes full
// main could pause and let the queue drain but that means that main must somehow start the threads running.
boost::lockfree::spsc_queue<SCF_THREAD_CONTROL, boost::lockfree::fixed_sized<false>, boost::lockfree::capacity<16000> > Tasks[MAX_RMG_THREADS];

// Called from main to setup thread tasks
void QueueThreadTask(int tid, SCF_THREAD_CONTROL &task)
{
    Tasks[tid].push(task);
}

// Called from threads to get what they are supposed to do
bool PopThreadTask(int tid, SCF_THREAD_CONTROL &task)
{
    bool ret = Tasks[tid].pop(task);
    return ret;
}

// Called from main to terminate all threads
void RmgTerminateThreads(void)
{
    BaseThread *T = BaseThread::getBaseThread(0);
    SCF_THREAD_CONTROL thread_controls[MAX_RMG_THREADS];
    for(int ist = 0;ist < MAX_RMG_THREADS;ist++)
    {
        thread_controls[ist].job = HYBRID_THREAD_EXIT;
        QueueThreadTask(ist, thread_controls[ist]);
    }
    T->run_thread_tasks(T->get_threads_per_node());
    if(ct.mpi_queue_mode) Rmg_Q->set_exitflag();
    if(ct.mpi_queue_mode) Rmg_Q->run_manager();
    T->thread_joinall();
}

#ifdef USE_NUMA
    #include "numa.h"

// Some versions of the numa library don't have the numa_bitmask_weight
// function so provide our own version to cover those cases.
unsigned int count_numamask_set_bits(const struct bitmask *mask)
{
    unsigned int count = 0;
    for(unsigned int idx = 0;idx < mask->size;idx++)
    {
        if(numa_bitmask_isbitset(mask, idx))count++;
    }
    return count;
}

#endif
#if GPU_ENABLED
#include <cuda.h>
#include <cuda_runtime.h>
#endif

#include <sys/types.h>
#include <sys/resource.h>
#include <sched.h>


// Main thread function specific to subprojects
void *run_threads(void *v) {

    BaseThreadControl *s;
    SCF_THREAD_CONTROL ss;
    s = (BaseThreadControl *)v;
    BaseThread *T = BaseThread::getBaseThread(0);
    int my_tid = s->tid;
 
#ifdef USE_NUMA
    struct bitmask *thread_cpumask = NULL;
    if(ct.use_numa) {
        thread_cpumask = numa_allocate_cpumask();

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
            numa_sched_setaffinity(0, thread_cpumask);
        }
        else if(pct.procs_per_host == pct.numa_nodes_per_host)
        {
            copy_bitmask_to_bitmask(pct.cpumask, thread_cpumask);
            numa_sched_setaffinity(0, thread_cpumask);
        }
        else if((pct.procs_per_host > pct.numa_nodes_per_host) && (pct.ncpus != pct.procs_per_host))
        {
            copy_bitmask_to_bitmask(pct.cpumask, thread_cpumask);
            numa_sched_setaffinity(0, thread_cpumask);
        }
        else
        {
            ct.use_numa = false;
        }

        numa_set_localalloc();

    }


#endif
    
#if GPU_ENABLED
    bool dev_set = false;
    cudaError_t cuerr;
#endif

    // Set up thread local storage
    rmg_set_tsd(s);


    while(1) {

        // See if a task if waiting. If not then sleep.
        if(!PopThreadTask(my_tid, ss)) 
        {
            T->thread_sleep();

            // When woken go to end of loop then circle around to pick up task
            continue;
#if GPU_ENABLED
            if(!dev_set) cudaSetDevice(ct.cu_dev); 
            dev_set = true;
#endif
        }

        // Set the project specific basetag into the main BaseThreadControl structure.
        T->set_thread_basetag(my_tid, ss.basetag);
        T->set_pptr(my_tid, (void *)&ss);

        // Switch that controls what we do
        switch(ss.job) {
            case HYBRID_ORBITALS_DOT_PRODUCT:       // Performs dot product for one pair of orbitals
                OrbitDotOrbitBlock(ss.basetag, ss.extratag2,(double *)ss.nv, (double *)ss.ns);
                break;
            case HYBRID_THETA_PHI:       // Performs Theta_ij * Phi_j
                ThetaPhiBlock(ss.basetag, ss.extratag2,(double *)ss.nv);
                break;
            case HYBRID_THREAD_EXIT:
                T->thread_exit();
                return NULL;

            default:
                break;
        }

    }

    return NULL;
}

