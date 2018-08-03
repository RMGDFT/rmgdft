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
#include <sys/time.h>
#include <sys/resource.h>
#include <boost/next_prior.hpp>
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
    //printf("Manager thread started.\n");

#ifdef USE_NUMA
    if(Q->cpumask)
    {
        numa_sched_setaffinity(0, Q->cpumask);
        numa_set_localalloc();
    }
#endif
    unsigned qcount = 0;
    mpi_queue_item_t qobj;


    boost::lockfree::spsc_queue<mpi_queue_item_t, boost::lockfree::fixed_sized<true>> *postedq = 
        new boost::lockfree::spsc_queue<mpi_queue_item_t, boost::lockfree::fixed_sized<true>>(Q->max_threads*54);



    // Here we loop and process requests
    while(1)
    {

        // Threads post requests into a boost single-producer, single-consumer queue.
        // Manager thread dequeues them and passes them on to MPI then we put the request info
        // into a separate queue that we then check for completion.
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
                if(!postedq->push(qobj)) {printf("Error: full queue.\n");fflush(NULL);exit(0);}
                qcount++;
            }
        }

        // Now check the posted items to see if any have completed        
        int breakcount = Q->max_threads;
        while(postedq->pop(qobj))
        {
            qcount--;
            int flag=false;
            MPI_Test(&qobj.req, &flag, MPI_STATUS_IGNORE);
            if(flag)
            {
                qobj.is_completed->store(true, std::memory_order_release);
                if((qobj.type == RMG_MPI_IRECV) || (qobj.type == RMG_MPI_ISEND))
                    qobj.group_count->fetch_sub(1, std::memory_order_release);
            }
            else
            {
                // Not complete so push it back on the queue and loop around again
                qcount++;
                postedq->push(qobj);
                breakcount--;
                if(breakcount == 0)break;
            }
            if(qcount == 0) break;
        }


        // No pending requests so spin CPU if running flag is set. Otherwise sleep.
        if(qcount == 0)
        {
            if(Q->exitflag.load(std::memory_order_acquire)) return;
            if(Q->running.load(std::memory_order_acquire))
            {
                if(Q->spin_manager) MpiQueue::spin(10);
            }
            else
            {
                Q->cvwait(manager_mutex, manager_cv, Q->running);
            }
        }
    }
}


MpiQueue::MpiQueue(int max_size, int max_threads, void *mask, void *newtopology, bool spin_manager_thread, bool spin_worker_threads)
{
#ifdef USE_NUMA
    this->cpumask = (struct bitmask *)mask;
#endif
#ifdef USE_HWLOC
    this->topology = (hwloc_topology_t *)newtopology;
#endif
    this->max_threads = max_threads;
    this->spin_manager = spin_manager_thread;
    this->spin_workers = spin_worker_threads;
    this->queue = new boost::lockfree::spsc_queue<mpi_queue_item_t, boost::lockfree::fixed_sized<true>> *[max_threads];
    for(int tid = 0;tid < max_threads;tid++)
    {
        this->queue[tid] = new boost::lockfree::spsc_queue<mpi_queue_item_t, boost::lockfree::fixed_sized<true>>(max_size);
    }
    this->exitflag.store(false, std::memory_order_release);
    this->running.store(true, std::memory_order_release);
    this->QueueManager = boost::thread(&manager_thread, this);
    this->QueueManager.detach();
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

void MpiQueue::set_exitflag(void)
{
    this->exitflag.store(true, std::memory_order_release);
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

void MpiQueue::waitall(std::atomic_bool *items, int n)
{
    if(this->spin_workers)
    {
        for(int i=0;i<n;i++)
        {
            while(!items[i].load(std::memory_order_acquire)) this->wait(10);
        }
    }
    else
    {
        for(int i=0;i<n;i++)
        {
            while(!items[i].load(std::memory_order_acquire)) std::this_thread::yield();
        }
    }
}

void MpiQueue::wait(int n)
{
    if(this->spin_workers)
    {
        MpiQueue::spin(n);
    }
    else
    {
        std::this_thread::yield();
    }
}

void MpiQueue::waitgroup(std::atomic_int &count)
{
    while(count.load(std::memory_order_acquire))
    {
        if(this->spin_workers)
        {
            MpiQueue::spin(5);
        }
        else
        {
            std::this_thread::yield();
        }
    }
}

void MpiQueue::spin(int n)
{
#if __GNUC__
    for (int i = 0; i < n; i++) asm volatile ("nop");
#endif
}
