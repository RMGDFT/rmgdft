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

#include <mutex>
#include <thread>
#include <condition_variable>
#include <queue>
#include <boost/lockfree/queue.hpp>
#include <boost/lockfree/spsc_queue.hpp>
#include <boost/thread/thread.hpp>
#include "BaseThread.h"
#include "MpiQueue.h"



std::condition_variable manager_cv;
std::mutex manager_mutex;

// Manager thread for queue mode. Works best with high thread
// per node counts when one CPU is dedicated to the manager thread.
void MpiQueue::manager_thread(MpiQueue *Q)
{
    printf("Manager thread started.\n");

#ifdef USE_NUMA
    if(Q->cpumask)
    {
        numa_sched_setaffinity(0, Q->cpumask);
        numa_set_localalloc();
    }
#endif
    unsigned qcount_large = 0;
    unsigned qcount_small = 0;
    mpi_queue_item_t qobj;


    boost::lockfree::spsc_queue<mpi_queue_item_t, boost::lockfree::fixed_sized<true>> *postedq_large = 
        new boost::lockfree::spsc_queue<mpi_queue_item_t, boost::lockfree::fixed_sized<true>>(Q->max_threads*54);
    boost::lockfree::spsc_queue<mpi_queue_item_t, boost::lockfree::fixed_sized<true>> *postedq_small = 
        new boost::lockfree::spsc_queue<mpi_queue_item_t, boost::lockfree::fixed_sized<true>>(Q->max_threads*54);



    // Here we loop and process requests
    while(1)
    {

        // Threads post requests into a boost multi-producer, multi-consumer queue.
        // We dequeue them and pass them on to MPI then put the request info into separate
        // small and large request queues that we then check for completion.
        for(int tid=0;tid < Q->max_threads;tid++)
        {
            while(Q->queue[tid]->pop(qobj))
            {
                if(qobj.type == RMG_MPI_IRECV)
                {
                    MPI_Irecv(qobj.buf, qobj.buflen, MPI_BYTE, qobj.target, qobj.mpi_tag, qobj.comm, &qobj.req);
                }
                else if(qobj.type == RMG_MPI_ISEND)
                {
                    MPI_Isend(qobj.buf, qobj.buflen, MPI_BYTE, qobj.target, qobj.mpi_tag, qobj.comm, &qobj.req);
                }
                else if(qobj.type == RMG_MPI_SUM)
                {
                    if(qobj.datatype == MPI_DOUBLE)
                        MPI_Iallreduce(MPI_IN_PLACE, qobj.buf, qobj.buflen, MPI_DOUBLE, MPI_SUM, qobj.comm, &qobj.req);
                    if(qobj.datatype == MPI_FLOAT)
                        MPI_Iallreduce(MPI_IN_PLACE, qobj.buf, qobj.buflen, MPI_FLOAT, MPI_SUM, qobj.comm, &qobj.req);
                    if(qobj.datatype == MPI_INT)
                        MPI_Iallreduce(MPI_IN_PLACE, qobj.buf, qobj.buflen, MPI_INT, MPI_SUM, qobj.comm, &qobj.req);
                }
                else 
                {
                    printf("Error: unknown MPI type.\n");fflush(NULL);exit(0);
                }
                // Now push it into our already posted queues which are only accessed by this thread so faster
                if(qobj.buflen>=256)
                {
                    if(!postedq_large->push(qobj)) {printf("Error: full queue.\n");fflush(NULL);exit(0);}
                    qcount_large++;
                }
                else
                {
                    if(!postedq_small->push(qobj)) {printf("Error: full queue.\n");fflush(NULL);exit(0);}
                    qcount_small++;
                }
            }
        }

        // Now check the posted items to see if any have completed        
        int breakcount = 0;
        while(postedq_large->pop(qobj))
        {
            qcount_large--;
            int flag;
            MPI_Test(&qobj.req, &flag, MPI_STATUS_IGNORE);
            if(flag)
            {
                qobj.is_completed->store(true, std::memory_order_release);
            }
            else
            {
                // Not complete so push it back on the queue and loop around again
                qcount_large++;
                postedq_large->push(qobj);
                breakcount++;
                if(breakcount >= 16)break;
            }
        }

        breakcount = 0;
        while(postedq_small->pop(qobj))
        {
            qcount_small--;
            int flag;
            MPI_Test(&qobj.req, &flag, MPI_STATUS_IGNORE);
            if(flag)
            {
                qobj.is_completed->store(true, std::memory_order_release);
            }
            else
            {
                // Not complete so push it back on the queue and loop around again
                qcount_small++;
                postedq_small->push(qobj);
                breakcount++;
                if(breakcount >= 16)break;
            }
        }

        // No pending requests so yield CPU if running flag is set. Otherwise sleep.
        if((qcount_large == 0) && (qcount_small == 0))
        {
            if(Q->running.load(std::memory_order_acquire))
            {
                //std::this_thread::yield();
                Q->spinwait(100);
            }
            else
            {
                Q->cvwait(manager_mutex, manager_cv, Q->running);
#ifdef USE_NUMA
                if(Q->cpumask) numa_sched_setaffinity(0, Q->cpumask);
#endif
            }
        }
    }
}


#ifdef USE_NUMA
MpiQueue::MpiQueue(int max_size, int max_threads, void *mask)
{
    this->cpumask = (struct bitmask *)mask;
#else
MpiQueue::MpiQueue(int max_size, int max_threads)
{
#endif
    this->max_threads = max_threads;
    this->queue = new boost::lockfree::spsc_queue<mpi_queue_item_t, boost::lockfree::fixed_sized<true>> *[max_threads];
    for(int tid = 0;tid < max_threads;tid++)
    {
        this->queue[tid] = new boost::lockfree::spsc_queue<mpi_queue_item_t, boost::lockfree::fixed_sized<true>>(max_size);
    }
    this->running.store(true, std::memory_order_release);
    this->QueueManager = boost::thread(&manager_thread, this);
    this->running.store(false, std::memory_order_release);
}

MpiQueue::~MpiQueue(void)
{
    for(int tid=0;tid < this->max_threads;tid++) delete this->queue[tid];
    delete [] this->queue;
}

bool MpiQueue::push(int tid, mpi_queue_item_t &item)
{
    return this->queue[tid]->push(item);
}

bool MpiQueue::pop(int tid, mpi_queue_item_t &item)
{
    return this->queue[tid]->pop(item);
}

void MpiQueue::run_manager(void)
{
    this->running.store(true, std::memory_order_release);
    manager_cv.notify_one();
}

void MpiQueue::stop_manager(void)
{
    this->running.store(false, std::memory_order_release);
}

void MpiQueue::cvwait(std::mutex &mut, std::condition_variable &cv, std::atomic_int &var)
{
    std::unique_lock<std::mutex> lk(mut);
    cv.wait(lk, [&] {return var.load(std::memory_order_acquire)==true;});
}

void MpiQueue::cvwait(std::mutex &mut, std::condition_variable &cv, std::atomic_bool &var)
{
    std::unique_lock<std::mutex> lk(mut);
    cv.wait(lk, [&] {return var.load(std::memory_order_acquire)==true;});
}

void MpiQueue::spinwaitall(std::atomic_bool *items, int n)
{
    for(int i=0;i<n;i++)
    {
        while(!items[i].load(std::memory_order_acquire)) this->spinwait(200);
    }
}

void MpiQueue::spinwait(int n)
{
#if __GNUC__
  for (int i = 0; i < n; i++) asm volatile ("nop");
#endif
}

