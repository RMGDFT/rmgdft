/*
 *
 * Copyright (c) 2016, Emil Briggs
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

#ifndef RMG_MpiQueue_H
#define RMG_MpiQueue_H 1

#include <mpi.h>
#include "BaseGrid.h"
#include "rmg_error.h"
#ifdef USE_NUMA
    #include <numa.h>
#endif

/* Type of async request passed to the mpi trade_images manager */
#define RMG_MPI_ISEND 1
#define RMG_MPI_IRECV 2
#define RMG_MPI_SUM 3

#if __cplusplus

//#include "BaseThread.h"
#include <atomic>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <boost/lockfree/queue.hpp>
#include <boost/thread/thread.hpp>
#include <boost/lockfree/spsc_queue.hpp>

typedef struct
{
    // Used by both manager and client threads
    std::atomic_bool *is_completed;

    // Async operation type, RMG_MPI_ISEND, RMG_MPI_IRECV, RMG_MPI_SUM
    int type;

    // Data type, not needed for RMG_MPI_ISEND and RMG_MPI_IRECV which use bytes
    MPI_Datatype datatype;

    // Used only by client threads
    bool is_unpacked;

    // Initialized by clients but never change
    void *buf;
    int buflen;
    int target;
    int target_index;

    // Actual tag passed to mpi routines.
    int mpi_tag;

    MPI_Comm comm;
    MPI_Request req;
} mpi_queue_item_t;


class MpiQueue {


private:

    boost::thread QueueManager;
    static void manager_thread(MpiQueue *Q);
    std::atomic<bool> running;
    int max_threads;

#ifdef USE_NUMA
    struct bitmask *cpumask{NULL};
#endif

public:

#ifdef USE_NUMA
    MpiQueue(int max_size, int max_threads, void *mask);
#else
    MpiQueue(int max_size, int max_threads);
#endif
    ~MpiQueue(void);

    bool push(int tid, mpi_queue_item_t &item);
    bool pop(int tid, mpi_queue_item_t &item);
    void cvwait(std::mutex &mut, std::condition_variable &cv, std::atomic_int &var);
    void cvwait(std::mutex &mut, std::condition_variable &cv, std::atomic_bool &var);
    void spinwaitall(std::atomic_bool *items, int n);
    void spinwait(int count);
    void run_manager(void);
    void stop_manager(void);
    boost::lockfree::spsc_queue<mpi_queue_item_t, boost::lockfree::fixed_sized<true>> **queue;
};

#endif
#endif
