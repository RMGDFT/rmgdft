#include "BaseThread.h"
#include "rmg_error.h"
#include "const.h"
using namespace std;


// Main thread control structure
BaseThread *thread_controls[MAX_SCF_THREADS];

// Constructor
BaseThread::BaseThread(int nthreads)
{

    int thread, retval, ncpus;
    BaseThread *s;

    if(!BaseThread::init_flag) {

        BaseThread::THREADS_PER_NODE = nthreads;

        // Should work on linux and AIX
        ncpus = sysconf( _SC_NPROCESSORS_ONLN );
        printf("Hybrid mode with %d threads and %d cores per node.\n", nthreads, ncpus);

        sem_init(&BaseThread::thread_sem, 0, 0);
        pthread_attr_init( &BaseThread::thread_attrs );
        pthread_attr_setscope( &BaseThread::thread_attrs, PTHREAD_SCOPE_SYSTEM );
        pthread_attr_setschedpolicy( &BaseThread::thread_attrs, SCHED_RR);

        // Create the thread specific data key
        pthread_key_create(&BaseThread::scf_thread_control_key, NULL);

        // Create the main sync barrier
        pthread_barrier_init(&BaseThread::run_barrier, NULL, BaseThread::THREADS_PER_NODE);

        BaseThread::init_flag = 1;

        // Create a set of long lived threads
        for(thread = 0;thread < BaseThread::THREADS_PER_NODE;thread++) {

            thread_controls[thread] = new BaseThread(nthreads);
            thread_controls[thread]->tid = thread;
            sem_init(&thread_controls[thread]->this_sync, 0, 0);
            retval = pthread_create(&threads[thread], &BaseThread::thread_attrs, 
                     &run_threads, (void *)thread_controls[thread]);

        }

    }

}

// Called when the main thread of execution is waiting for a set of threads to finish
void BaseThread::wait_for_threads(int jobs) {
    int idx;
    for(idx = 0;idx < jobs;idx++) {
        sem_wait(&BaseThread::thread_sem);
    }
}

// Wakes jobs sleeping threads starting from tid=0 and counting up
// jobs must be less than THREADS_PER_NODE
void BaseThread::wake_threads(int jobs) {

    int thread;

    if(jobs > BaseThread::THREADS_PER_NODE) {
        // If this happens it is a bug
        rmg_error_handler("More jobs than available threads scheduled\n");
    }

    pthread_mutex_lock(&BaseThread::job_mutex);
    pthread_mutex_unlock(&BaseThread::job_mutex);

    for(thread = 0;thread < jobs;thread++) {
        sem_post(&thread_controls[thread]->this_sync);
    }

}

// Initialization function for barriers
void BaseThread::scf_barrier_init(int nthreads) {
    pthread_barrier_init(&BaseThread::scf_barrier, NULL, nthreads);
}

// Blocks all threads until nthreads specified in the init call have reached this point
void BaseThread::scf_barrier_wait(void) {
    if(!BaseThread::in_threaded_region) return;
    pthread_barrier_wait(&BaseThread::scf_barrier);
}

// Termination function
void BaseThread::scf_barrier_destroy(void) {
    pthread_barrier_destroy(&scf_barrier);
}
// Initializes the key and sets a value into it

void BaseThread::scf_tsd_init(void) {
 pthread_key_create(&scf_thread_control_key, NULL);
}

// Sets a value into the key
void BaseThread::scf_tsd_set_value(void *s) {
     pthread_setspecific(scf_thread_control_key, s);
}

// Deletes the key
void BaseThread::scf_tsd_delete(void) {
 pthread_key_delete(BaseThread::scf_thread_control_key);
}


// Reads the basetag from the thread specific data. Returns 0 if we are not in
// a parallel region. Meant to be called from a thread.
int BaseThread::get_thread_basetag(void) {

    BaseThread *ss;
    if(!BaseThread::in_threaded_region) return 0;
    ss = (BaseThread *)pthread_getspecific(BaseThread::scf_thread_control_key);
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
    ss = (BaseThread *)pthread_getspecific(BaseThread::scf_thread_control_key);
    if(!ss) return NULL;
    return ss;
}

// Reads the tid from the thread specific data. Returns -1 if we are not in 
// a parallel region
int BaseThread::get_thread_tid(void) {

    BaseThread *ss;

    if(!BaseThread::in_threaded_region) return -1;
    ss = (BaseThread *)pthread_getspecific(BaseThread::scf_thread_control_key);
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

}

void BaseThread::enter_threaded_region(void) {
    BaseThread::in_threaded_region = 1;
}
void BaseThread::leave_threaded_region(void) {
    BaseThread::in_threaded_region = 0;
}

void BaseThread::RMG_MPI_lock(void) {
    pthread_mutex_lock(&BaseThread::mpi_mutex);
}

void BaseThread::RMG_MPI_unlock(void) {
    pthread_mutex_unlock(&BaseThread::mpi_mutex);
}

void BaseThread::RMG_MPI_thread_order_lock(void) {
   int tid, i1, ntid;
   tid = get_thread_tid();

   if(tid < 0) {
       rmg_error_handler("\nError in RMG_MPI_thread_order_lock. Terminating.\n");
   }

   while(1) {

       // Acquire the lock
       pthread_mutex_lock(&BaseThread::thread_order_mutex);

       // See if it's our turn
       i1 = BaseThread::mpi_thread_order_counter % BaseThread::THREADS_PER_NODE;
       if(i1 == tid) {
           // Raise priority of next thread
           ntid = i1 + 1;
           if(ntid < BaseThread::THREADS_PER_NODE) {
               pthread_setschedprio(BaseThread::threads[ntid], -19);
           }
           return;
       }

       pthread_mutex_unlock(&BaseThread::thread_order_mutex);
       sched_yield();

   }

}

void BaseThread::RMG_MPI_thread_order_unlock(void) {

  int tid;
  tid = get_thread_tid();

  BaseThread::mpi_thread_order_counter++;
  pthread_setschedprio(BaseThread::threads[tid], -1);
  pthread_mutex_unlock(&BaseThread::thread_order_mutex);

}

extern rmg_double_t timings[LAST_TIME];
pthread_mutex_t BaseThread::timings_mutex = PTHREAD_MUTEX_INITIALIZER;
void BaseThread::rmg_timings (int what, rmg_double_t time)
{
    pthread_mutex_lock(&BaseThread::timings_mutex);
    if(BaseThread::in_threaded_region) {
        timings[what] += time / BaseThread::THREADS_PER_NODE;
    }
    else {
        timings[what] += time;
    }
    pthread_mutex_unlock(&BaseThread::timings_mutex);
}                               /* end rmg_timings */


int BaseThread::is_loop_over_states(void)
{
    BaseThread *ss;
    if(!BaseThread::in_threaded_region) return 0;
    ss = (BaseThread *)pthread_getspecific(BaseThread::scf_thread_control_key);
    if(!ss) return 0;
    return 1;
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

// Pthread identifier
pthread_t pthread_tid;

// Pointer to project specific data structure
void *pptr;

// Used to implement a local barrier for all threads inside of the run_threads function
pthread_barrier_t BaseThread::run_barrier;

// Used to implement a local barrier inside of the scf loops
pthread_barrier_t BaseThread::scf_barrier;
pthread_mutex_t BaseThread::job_mutex = PTHREAD_MUTEX_INITIALIZER;

// These are used to synchronize the main process and the worker threads
sem_t BaseThread::thread_sem;

// This is used when running with MPI_THREAD_SERIALIZED to ensure 
// proper serialization
pthread_mutex_t BaseThread::mpi_mutex = PTHREAD_MUTEX_INITIALIZER;

pthread_attr_t BaseThread::thread_attrs;
pthread_t BaseThread::threads[MAX_SCF_THREADS];
volatile int BaseThread::in_threaded_region = 0;

// These are used to ensure thread ordering
volatile int BaseThread::mpi_thread_order_counter = 0;
pthread_mutex_t BaseThread::thread_order_mutex = PTHREAD_MUTEX_INITIALIZER;

// Used for accessing thread specific data
pthread_key_t BaseThread::scf_thread_control_key;

// Synchronization semaphore for this instance
sem_t this_sync;

// Basetag
int basetag;

extern "C" void wait_for_threads(int jobs)
{
    BaseThread B(0);
    B.wait_for_threads(jobs);
}

extern "C" void wake_threads(int jobs)
{
    BaseThread B(0);
    B.wake_threads(jobs);
}

extern "C" void scf_barrier_init(int nthreads) 
{
    BaseThread B(0);
    B.scf_barrier_init(nthreads);
}

extern "C" void scf_barrier_wait(void) 
{
    BaseThread B(0);
    B.scf_barrier_wait();
}

extern "C" void scf_barrier_destroy(void) 
{
    BaseThread B(0);
    B.scf_barrier_destroy();
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
    BaseThread *B;
    // This is not a leak. We want B to live forever
    B = new BaseThread(nthreads);
}

extern "C" void set_cpu_affinity(int tid)
{
    BaseThread B(0);
    B.set_cpu_affinity(tid);
}

extern "C" void enter_threaded_region(void)
{
    BaseThread B(0);
    B.enter_threaded_region();
}

extern "C" void leave_threaded_region(void)
{
    BaseThread B(0);
    B.leave_threaded_region();
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

extern "C" void RMG_MPI_thread_order_lock(void) 
{
   BaseThread B(0);
   B.RMG_MPI_thread_order_lock(); 
}

extern "C" void RMG_MPI_thread_order_unlock(void) 
{
   BaseThread B(0);
   B.RMG_MPI_thread_order_unlock(); 
}

extern "C" void rmg_timings (int what, rmg_double_t time)
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

