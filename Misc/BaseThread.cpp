#include "BaseThread.h"
#include "rmg_error.h"
#include "const.h"
using namespace std;


// Main thread control structure
BaseThread *thread_controls[MAX_RMG_THREADS];

// Condition variable and mutex for threads
std::mutex BaseThread::thread_mutex;
std::condition_variable BaseThread::thread_cv;

// Condition variable and mutex for main
std::mutex BaseThread::main_mutex;
std::condition_variable BaseThread::main_cv;


// Constructor
BaseThread::BaseThread(int nthreads)
{

    int thread, retval, ncpus;
    BaseThread *s;

    if(nthreads > MAX_RMG_THREADS)
        rmg_error_handler (__FILE__, __LINE__, "Too many threads requested. Change MAX_RMG_THREADS and recompile if needed.");

    if(!BaseThread::init_flag) {

        BaseThread::THREADS_PER_NODE = nthreads;

        // Should work on linux and AIX
        ncpus = sysconf( _SC_NPROCESSORS_ONLN );
        printf("Hybrid mode with %d threads and %d cores per node.\n", nthreads, ncpus);

        BaseThread::init_flag = 1;

        // Create a set of long lived threads
        for(thread = 0;thread < BaseThread::THREADS_PER_NODE;thread++) {

            thread_controls[thread] = new BaseThread(nthreads);
            thread_controls[thread]->tid = thread;
            threads[thread] = new boost::thread(run_threads, (void *)thread_controls[thread]);

        }

    }

}


// Wakes jobs sleeping threads starting from tid=0 and counting up
// jobs must be less than THREADS_PER_NODE
void BaseThread::run_thread_tasks(int jobs) {

    int thread;

    if(jobs > BaseThread::THREADS_PER_NODE) {
        // If this happens it is a bug
        rmg_error_handler (__FILE__, __LINE__, "More jobs than available threads scheduled\n");
    }

    BaseThread::scf_barrier = new boost::barrier(jobs);
    std::unique_lock<std::mutex> lk(BaseThread::thread_mutex);
    BaseThread::in_threaded_region = true;
    BaseThread::jobs = jobs;
    BaseThread::thread_cv.notify_all();
    BaseThread::main_cv.wait(lk, [] {return BaseThread::jobs == 0;});
    BaseThread::in_threaded_region = false;
    delete(BaseThread::scf_barrier);
 
}

void BaseThread::thread_sleep(void)
{
    std::unique_lock<std::mutex> lk(BaseThread::thread_mutex);
    if(BaseThread::jobs > 0) {
        BaseThread::jobs--;
        BaseThread::main_cv.notify_one();
    }
    BaseThread::thread_cv.wait(lk);
}

// Blocks all threads until nthreads specified in the init call have reached this point
void BaseThread::thread_barrier_wait(void) {
    int nt = BaseThread::THREADS_PER_NODE;
    if(!BaseThread::in_threaded_region) return;
    BaseThread::scf_barrier->wait();
}

// Reads the basetag from the thread specific data. Returns 0 if we are not in
// a parallel region. Meant to be called from a thread.
int BaseThread::get_thread_basetag(void) {

    BaseThread *ss;
    if(!BaseThread::in_threaded_region) return 0;
    ss = rmg_get_tsd();
    if(!ss) return 0;

    return ss->basetag;

}

// Sets a basetag value. Meant to be called from the main program.
void BaseThread::set_thread_basetag(int tid, int tag)
{
    thread_controls[tid]->basetag = tag;
}

// Gets the threads control structure pointer
BaseThread *BaseThread::get_thread_control(void) {
    BaseThread *ss;
    if(!BaseThread::in_threaded_region) return NULL;
    ss = rmg_get_tsd();
    if(!ss) return NULL;
    return ss;
}

// Reads the tid from the thread specific data. Returns -1 if we are not in 
// a parallel region
int BaseThread::get_thread_tid(void) {

    BaseThread *ss;

    if(!BaseThread::in_threaded_region) return -1;
    ss = rmg_get_tsd();
    if(!ss) return -1;

    return ss->tid;
}

#if GPU_ENABLED
// Gets thread cstream
//cudaStream_t *BaseThread::get_thread_cstream(void) {

//    BaseThread *ss;
//    if(!in_threaded_region) return &ct.cuda_stream;  // Return main process stream
//    ss = (BaseThread *)pthread_getspecific(scf_thread_control_key);
//    if(!ss) return NULL;
//    return &ss->cstream;
//}
#endif

// Used for positioning and setting processor affinity. For now assumes that
// THREADS_PER_NODE is an even multiple of ct.ncpus. If this is not true it
// does not attemp to schedule
void BaseThread::set_cpu_affinity(int tid)
{
#if __linux__
    int s;
    cpu_set_t cpuset;
    pthread_t thread;

    if(BaseThread::THREADS_PER_NODE % sysconf( _SC_NPROCESSORS_ONLN )) return;

    s = tid % BaseThread::THREADS_PER_NODE;


    // Set affinity mask
    CPU_ZERO(&cpuset);
    CPU_SET(s, &cpuset);

    thread = pthread_self();
    pthread_setaffinity_np(thread, sizeof(cpu_set_t), &cpuset);
#endif
}

void BaseThread::RMG_MPI_lock(void) {
    BaseThread::mpi_mutex.lock();
}

void BaseThread::RMG_MPI_unlock(void) {
    BaseThread::mpi_mutex.unlock();
}

extern double timings[LAST_TIME];
std::mutex BaseThread::timings_mutex;
void BaseThread::rmg_timings (int what, double time)
{
    BaseThread::timings_mutex.lock();
    if(BaseThread::in_threaded_region) {
        timings[what] += time / BaseThread::THREADS_PER_NODE;
    }
    else {
        timings[what] += time;
    }
    BaseThread::timings_mutex.unlock();
}                               /* end rmg_timings */

int BaseThread::is_loop_over_states(void)
{
    BaseThread *ss;
    if(!BaseThread::in_threaded_region) return 0;
    ss = rmg_get_tsd();
    if(!ss) return 0;
    return 1;
}

int BaseThread::get_threads_per_node(void)
{
    if(BaseThread::THREADS_PER_NODE == 0)
        rmg_error_handler (__FILE__, __LINE__, "Threads not initialized yet");
    return BaseThread::THREADS_PER_NODE;
}

void BaseThread::set_pptr(int tid, void *p)
{
    thread_controls[tid]->pptr = p;
}

// Init flag
int BaseThread::init_flag=0;

// Threads to use on each MPI node. Default is 1
int BaseThread::THREADS_PER_NODE=1;

// Thread ID number assigned by us
int tid;

// Pointer to project specific data structure
void *pptr;

#if GPU_ENABLED
    // Cuda device stream
    void *cstream;
#endif

// Used to implement a local barrier inside of the scf loops
boost::barrier *BaseThread::scf_barrier;

// This is used when running with MPI_THREAD_SERIALIZED to ensure 
// proper serialization
std::mutex BaseThread::mpi_mutex;

boost::thread *BaseThread::threads[MAX_RMG_THREADS];
std::atomic<int> BaseThread::jobs (0);
std::atomic<bool> BaseThread::in_threaded_region (false);

// These are used to ensure thread ordering
volatile int BaseThread::mpi_thread_order_counter = 0;
std::mutex BaseThread::thread_order_mutex;

// Basetag
int basetag;


// Non member functions used for handling thread specific data
static boost::thread_specific_ptr<BaseThread> my_ptr;
void rmg_set_tsd(BaseThread *p)
{
    if(!my_ptr.get()) {
       my_ptr.reset(p);
    }
}
BaseThread *rmg_get_tsd(void) {
    return my_ptr.get();
}

extern "C" void run_thread_tasks(int jobs)
{
    BaseThread B(0);
    B.run_thread_tasks(jobs);
}

extern "C" void thread_barrier_wait(void) 
{
    BaseThread B(0);
    B.thread_barrier_wait();
}

extern "C" int get_thread_basetag(void)
{
    BaseThread B(0);
    return B.get_thread_basetag();
}

extern "C" void set_thread_basetag(int tid, int tag)
{
    BaseThread B(0);
    B.set_thread_basetag(tid, tag);
}

extern "C" BaseThread *get_thread_control(void)
{
    BaseThread B(0);
    return B.get_thread_control();
}

extern "C" int get_thread_tid(void)
{
    BaseThread B(0);
    return B.get_thread_tid();
}


#if GPU_ENABLED
//extern "C" cudaStream_t *get_thread_cstream(void)
//{
//    BaseThread B(0);
//    return B.get_thread_cstream();
//}
#endif

extern "C" void init_HYBRID_MODEL(int nthreads)
{
    BaseThread B(nthreads);
    // This is not a leak. We want B to live forever
    //B = new BaseThread(nthreads);
}

extern "C" void set_cpu_affinity(int tid)
{
    BaseThread B(0);
    B.set_cpu_affinity(tid);
}

extern "C" void RMG_MPI_lock(void)
{
    BaseThread B(0);
    B.RMG_MPI_lock();
}

extern "C" void RMG_MPI_unlock(void)
{
    BaseThread B(0);
    B.RMG_MPI_unlock();
}


extern "C" void rmg_timings (int what, double time)
{
    if(BaseThread::init_flag)
    {
        BaseThread B(0);
        B.rmg_timings(what, time);
    }
}

extern "C" int is_loop_over_states(void) 
{
    BaseThread B(0);
    return B.is_loop_over_states();
}

extern "C" void set_pptr(int tid, void *p)
{
    BaseThread B(0);
    B.set_pptr(tid, p);
}

