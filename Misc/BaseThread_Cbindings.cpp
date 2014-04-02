#include "BaseThread.h"
#include "rmg_error.h"
#include "const.h"


// Non member functions used for handling thread specific data
extern "C" void run_thread_tasks(int jobs)
{
    BaseThread *B = BaseThread::getBaseThread(0);
    B->run_thread_tasks(jobs);
}

extern "C" void thread_barrier_wait(void) 
{
    BaseThread *B = BaseThread::getBaseThread(0);
    B->thread_barrier_wait();
}

extern "C" int get_thread_basetag(void)
{
    BaseThread *B = BaseThread::getBaseThread(0);
    return B->get_thread_basetag();
}

extern "C" void set_thread_basetag(int tid, int tag)
{
    BaseThread *B = BaseThread::getBaseThread(0);
    B->set_thread_basetag(tid, tag);
}

extern "C" int get_thread_tid(void)
{
    BaseThread *B = BaseThread::getBaseThread(0);
    return B->get_thread_tid();
}

extern "C" void init_HYBRID_MODEL(int nthreads)
{
    BaseThread *B = BaseThread::getBaseThread(nthreads);
    // This is not a leak. We want B to live forever
    //B = new BaseThread(nthreads);
}

extern "C" void set_cpu_affinity(int tid)
{
    BaseThread *B = BaseThread::getBaseThread(0);
    B->set_cpu_affinity(tid);
}

extern "C" void RMG_MPI_lock(void)
{
    BaseThread *B = BaseThread::getBaseThread(0);
    B->RMG_MPI_lock();
}

extern "C" void RMG_MPI_unlock(void)
{
    BaseThread *B = BaseThread::getBaseThread(0);
    B->RMG_MPI_unlock();
}

extern "C" int is_loop_over_states(void) 
{
    BaseThread *B = BaseThread::getBaseThread(0);
    return B->is_loop_over_states();
}

extern "C" void set_pptr(int tid, void *p)
{
    BaseThread *B = BaseThread::getBaseThread(0);
    B->set_pptr(tid, p);
}

