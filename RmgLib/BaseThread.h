#ifndef RMG_BaseThread_H
#define RMG_BaseThread_H 1


#include "rmg_error.h"

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

// Singleton class used to manage a thread pool
class BaseThread {

private:

    // Threads to use on each MPI node
    int THREADS_PER_NODE;

    // Used to implement a local barrier inside of the scf loops
    boost::barrier *scf_barrier;

    // This is used when running with MPI_THREAD_SERIALIZED to ensure 
    // proper serialization
    std::mutex mpi_mutex;

    boost::thread *threads[MAX_RMG_THREADS];
    std::atomic<int> jobs;
    std::atomic<bool> in_threaded_region;

    // These are used to ensure thread ordering
    volatile int mpi_thread_order_counter;
    std::mutex thread_order_mutex;

    // Base tag
    int basetag;

    // Initialization flag
    static int init_flag;

    // Points to the single instance
    static BaseThread *instance;

    // Thread function pointer
    static void *(*funcptr)(void *);

    // Private constructur
    BaseThread(int nthreads);

public:

    static BaseThread *getBaseThread(int nthreads);

    void RegisterThreadFunction(void *(*funcptr)(void *));

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

    void run_thread_tasks(int jobs);
    void thread_barrier_wait(void);
    int get_thread_basetag(void);
    void set_thread_basetag(int tid, int tag);
    BaseThreadControl *get_thread_control(void);
    int get_thread_tid(void);
    void set_cpu_affinity(int tid);
    void RMG_MPI_lock(void);
    void RMG_MPI_unlock(void);
    void set_pptr(int tid, void *p);
    int is_loop_over_states(void);
    int get_threads_per_node(void);
    void thread_sleep(void);

};

void rmg_set_tsd(BaseThreadControl *p);
BaseThreadControl *rmg_get_tsd(void);

#endif
#endif
