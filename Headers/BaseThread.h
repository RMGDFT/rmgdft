#ifndef RMG_BaseThread_H
#define RMG_BaseThread_H 1


#include <signal.h>
#include <sys/types.h>
#include <unistd.h>
#include <mpi.h>

#if GPU_ENABLED
#include <cuda.h>
#include <cuda_runtime.h>
#endif
#include "params.h"
#include "rmgtypes.h"
#include "rmg_error.h"

// Maximum number of Rmg threads. Adjust based on hardware resources.
#define MAX_RMG_THREADS 32


// Project specific thread function
void *run_threads(void *s);


#ifdef __cplusplus

#include <thread>
#include <mutex>
#include <atomic>
#include <condition_variable>
#include <boost/thread.hpp>
#include <boost/thread/tss.hpp>

class BaseThread {

private:

    // Threads to use on each MPI node
    static int THREADS_PER_NODE;

    // Used to implement a local barrier inside of the scf loops
    static boost::barrier *scf_barrier;

    // This is used when running with MPI_THREAD_SERIALIZED to ensure 
    // proper serialization
    static std::mutex mpi_mutex;

    static boost::thread *threads[MAX_RMG_THREADS];
    static std::atomic<int> jobs;
    static std::atomic<bool> in_threaded_region;

    // These are used to ensure thread ordering
    static volatile int mpi_thread_order_counter;
    static std::mutex thread_order_mutex;

    // Used for timing information
    static std::mutex timings_mutex;

    // Base tag
    int basetag;


public:

    // Initialization flag
    static int init_flag;

    // Condition variable and mutex for threads
    static std::mutex thread_mutex;
    static std::condition_variable thread_cv;

    // Condition variable and mutex for main
    static std::mutex main_mutex;
    static std::condition_variable main_cv;

    // Thread ID number assigned by us
    int tid;

    // Pointer to project specific data structure
    void *pptr;

#if GPU_ENABLED
    // Cuda device stream
    void *cstream;
#endif

    BaseThread(int nthreads);
    void run_thread_tasks(int jobs);
    void thread_barrier_wait(void);
    int get_thread_basetag(void);
    void set_thread_basetag(int tid, int tag);
    BaseThread *get_thread_control(void);
    int get_thread_tid(void);
    //cudaStream_t *get_thread_cstream(void);
    void set_cpu_affinity(int tid);
    void RMG_MPI_lock(void);
    void RMG_MPI_unlock(void);
    void rmg_timings (int what, rmg_double_t time);
    void set_pptr(int tid, void *p);
    int is_loop_over_states(void);
    int get_threads_per_node(void);
    void thread_sleep(void);

};

void rmg_set_tsd(BaseThread *p);
BaseThread *rmg_get_tsd(void);

#endif
#endif
