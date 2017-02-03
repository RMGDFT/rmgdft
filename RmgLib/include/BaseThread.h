/*
 *
 * Copyright (c) 2013, Emil Briggs
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the <organization> nor the
 *       names of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.

 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * 
*/

#ifndef RMG_BaseThread_H
#define RMG_BaseThread_H 1


#include "rmg_error.h"
#include "ThreadBarrier.h"

// Maximum number of Rmg threads. Adjust based on hardware resources.
#define MAX_RMG_THREADS 32


#ifdef __cplusplus

#include <thread>
#include <mutex>
#include <atomic>
#include <condition_variable>
#include <boost/thread.hpp>
#include <boost/thread/tss.hpp>
#include "BaseThreadControl.h"
#include <boost/lockfree/queue.hpp>
#include "MpiQueue.h"

// Singleton class used to manage a thread pool
class BaseThread {

private:

    // Threads to use on each MPI node
    int THREADS_PER_NODE;
    std::atomic<int> ACTIVE_THREADS_PER_NODE;

    // This is used when running with MPI_THREAD_SERIALIZED to ensure 
    // proper serialization
    static std::mutex mpi_mutex;

    std::atomic<int> jobs;
    std::atomic<bool> in_threaded_region;

    // Parent grid_comm returned from get_grid_comm when not in threaded region
    MPI_Comm parent_grid_comm;

    // Pool of communicators not tied to a specific thread. Pool size is set equal to the
    // number of threads. Typically used to ensure that electronic orbitals being processed
    // by parallel threads use a unique communicator across nodes.
    MPI_Comm *comm_pool;

    // Job queue. Only truly independent thread tasks should be queued. Mainly used for MgEigState
    // when using mpi_queue_mode.
    boost::lockfree::queue<BaseThreadControl *, boost::lockfree::fixed_sized<true>> *jobqueue;

    // Initialization flag
    static int init_flag;

    // Points to the single instance
    static BaseThread *instance;

    // Thread function pointer
    static void *(*funcptr)(void *);

    // Private constructor
    BaseThread(int nthreads);

public:

    static BaseThread *getBaseThread(int nthreads);

    void RegisterThreadFunction(void *(*funcptr)(void *), MPI_Comm &comm);

    // Condition variable and mutex for threads
    static std::mutex thread_mutex;
    static std::condition_variable thread_cv;

    // Condition variable and mutex for main
    static std::mutex main_mutex;
    static std::condition_variable main_cv;

    // Thread ID number assigned by us
    int tid;

    // Used to implement local barriers inside of the scf loops
    //boost::barrier *scf_barrier;
    Thread_barrier_t *barrier;

    void run_thread_tasks(int jobs, MpiQueue *Queue);
    void run_thread_tasks(int jobs);
    void thread_barrier_wait(bool spin);
    int get_thread_basetag(void);
    void set_thread_basetag(int tid, int tag);
    BaseThreadControl *get_thread_control(void);
    int get_thread_tid(void);

    // Called from thread returns grid_comm specific to that thread
    // which is a duplicate of parent_grid_comm with a different context
    MPI_Comm get_grid_comm(void);

    // For a given value of index returns a unique communicator across nodes for each
    // value of index.
    MPI_Comm get_unique_comm(int index);

    void set_cpu_affinity(int tid, int procs_per_node, int local_rank);
    void RMG_MPI_lock(void);
    void RMG_MPI_unlock(void);
    void set_pptr(int tid, void *p);
    void *get_pptr(int tid);
    int is_loop_over_states(void);
    int get_threads_per_node(void);
    int get_active_threads_per_node(void);
    void thread_sleep(void);
    void thread_exit(void);
    void wake_all_threads(void);
    void thread_joinall(void);

};

void rmg_set_tsd(BaseThreadControl *p);
void rmg_release_tsd(void);
BaseThreadControl *rmg_get_tsd(void);

#endif
#endif
