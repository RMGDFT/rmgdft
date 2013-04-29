/*** QMD-MGDFT/main.h *****
 * NAME
 *   Ab initio real space code with multigrid acceleration
 *   Quantum molecular dynamics package.
 *   Version: 2.1.5
 * COPYRIGHT
 *   Copyright (C) 1995  Emil Briggs
 *   Copyright (C) 1998  Emil Briggs, Charles Brabec, Mark Wensell, 
 *                       Dan Sullivan, Chris Rapcewicz, Jerzy Bernholc
 *   Copyright (C) 2001  Emil Briggs, Wenchang Lu,
 *                       Marco Buongiorno Nardelli,Charles Brabec, 
 *                       Mark Wensell,Dan Sullivan, Chris Rapcewicz,
 *                       Jerzy Bernholc
 */

/*


            hybrid.c

  Functions and private data structures required for hybrid MPI-threads
  parallel version of the code.


 */

#define _GNU_SOURCE

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <pthread.h>
#include <signal.h>
#include <sys/types.h>
#include <unistd.h>
#include <semaphore.h>

#include "main.h"
#include "hybrid.h"


#if PAPI_PERFMON
#undef kill
#include <papi.h>
volatile double PTHREAD_PAPI_COUNTERS[MAX_SCF_THREADS];
#endif

#if HYBRID_MODEL

// Main thread control structure
SCF_THREAD_CONTROL thread_control[MAX_SCF_THREADS];

// Used to implement a local barrier for all threads inside of the run_threads function
static pthread_barrier_t run_barrier;

// Used to implement a local barrier inside of the scf loops
static pthread_barrier_t scf_barrier;
static pthread_mutex_t job_mutex = PTHREAD_MUTEX_INITIALIZER;

// These are used to synchronize the main process and the worker threads
sem_t thread_sem;
volatile int job_count=0;
static pthread_mutex_t sync_mutex = PTHREAD_MUTEX_INITIALIZER;

// This is used when running with MPI_THREAD_SERIALIZED to ensure 
// proper serialization
static pthread_mutex_t mpi_mutex = PTHREAD_MUTEX_INITIALIZER;

static pthread_attr_t thread_attrs;
static pthread_t threads[MAX_SCF_THREADS];
volatile int in_threaded_region = 0;
static void run_threads(SCF_THREAD_CONTROL *s);

// These are used to ensure thread ordering
volatile int mpi_thread_order_counter = 0;
static pthread_mutex_t thread_order_mutex = PTHREAD_MUTEX_INITIALIZER;

// Used for accessing thread specific data
pthread_key_t scf_thread_control_key;

extern int TRADE_GRID_EDGES;
extern int GRID_MAX1;
extern int GRID_MAX2;

// Initialization function called by main process
void init_HYBRID_MODEL(void) {

    int thread, retval;
    SCF_THREAD_CONTROL *s;

    // Should work on linux and AIX
    ct.ncpus = sysconf( _SC_NPROCESSORS_ONLN );
    printf("Hybrid mode with %d threads and %d cores per node.\n", ct.THREADS_PER_NODE, ct.ncpus);

    sem_init(&thread_sem, 0, 0);

    pthread_attr_init( &thread_attrs );
    pthread_attr_setscope( &thread_attrs, PTHREAD_SCOPE_SYSTEM );
    pthread_attr_setschedpolicy( &thread_attrs, SCHED_RR);

    // Create the thread specific data key
    scf_tsd_init();

    // Create the main sync barrier
    pthread_barrier_init(&run_barrier, NULL, ct.THREADS_PER_NODE);

    // Here we create a set of long lived threads
    for(thread = 0;thread < ct.THREADS_PER_NODE;thread++) {

        s = &thread_control[thread];
        thread_control[thread].tid = thread;
        sem_init(&thread_control[thread].sync, 0, 0);
        retval = pthread_create(&threads[thread], &thread_attrs, (void *)run_threads, (void *)&thread_control[thread]);
    }

#if GPU_ENABLED
    // Pop context then any other threads need to push/pop to access GPU
//    cuCtxPopCurrent(&ct.cu_context);
#endif

    
}

// Main thread function
void run_threads(SCF_THREAD_CONTROL *s) {

    int retval, alloc;

#if GPU_ENABLED
    cudaError_t cuerr;
#endif

    set_cpu_affinity(s->tid);

#if PAPI_PERFMON
    int EventSet = PAPI_NULL;
    int PAPI_event;
    long long Papi_values[4];

    Papi_values[0] = 0;
    Papi_values[1] = 0;
    Papi_values[2] = 0;
    Papi_values[3] = 0;

    if (PAPI_create_eventset(&EventSet) != PAPI_OK) {
         error_handler ("Cannot create PAPI event set.\n");
    }
    if (PAPI_add_event(EventSet, PAPI_DP_OPS) != PAPI_OK) {
         error_handler ("Cannot add PAPI event.\n");
    }
   
    PAPI_start(EventSet); 

#endif

    s->pthread_tid = pthread_self();
    pthread_setspecific(scf_thread_control_key, (void *)s);

    alloc = (pct.PX0_GRID + 2*MAX_TRADE_IMAGES) * (pct.PY0_GRID + 2*MAX_TRADE_IMAGES) * (pct.PZ0_GRID + 2*MAX_TRADE_IMAGES);

    retval = MPI_Alloc_mem(sizeof(REAL) * (alloc + 64) , MPI_INFO_NULL, &s->trade_buf);
    if(retval != MPI_SUCCESS) {
         error_handler("Error in MPI_Alloc_mem.\n");
    }

#if GPU_ENABLED
    cudaSetDevice(ct.cu_dev); 
    if(cudaSuccess != (cuerr = cudaStreamCreate(&s->cstream))) {
        fprintf (stderr, "Error: cudaStreamCreate failed for: threads setup\n");
        exit(-1);
    }
    if( cudaSuccess != cudaMallocHost((void **)&s->gpu_host_temp1, (pct.PX0_GRID + 4) * (pct.PY0_GRID + 4) * (pct.PZ0_GRID + 4) * sizeof(REAL) )){
        fprintf (stderr, "Error: cudaMallocHost failed for: threads gpu_host_temp\n");
        exit(-1);
    }
    if( cudaSuccess != cudaMallocHost((void **)&s->gpu_host_temp2, (pct.PX0_GRID + 4) * (pct.PY0_GRID + 4) * (pct.PZ0_GRID + 4) * sizeof(REAL) )){
        fprintf (stderr, "Error: cudaMallocHost failed for: threads gpu_host_temp\n");
        exit(-1);
    }
    cuMemHostRegister ( s->trade_buf, sizeof(REAL) * (alloc + 64), CU_MEMHOSTREGISTER_PORTABLE);
#endif

    // Wait until everyone gets here
    pthread_barrier_wait(&run_barrier);

    while(1) {

        // We sleep forever or until we get a signal that wakes us up
        sem_wait(&s->sync);

        // Switch that controls what we do
        switch(s->job) {
            case HYBRID_EIG:       // Performs a single multigrid sweep over an orbital
               mg_eig_state_driver(s->sp, 0, s->vtot, ct.mg_eig_precision);
               break;
            case HYBRID_SKIP:
               break;
#if GAMMA_PT
            case HYBRID_SUBDIAG_APP_AB:
               subdiag_app_AB_one(s->sp, s->p1, s->p2, s->vtot);
               break;
            case HYBRID_SUBDIAG_APP_A:
               subdiag_app_A_one(s->sp, s->p1, s->p2, s->vtot);
               break;
            case HYBRID_SUBDIAG_APP_B:
               subdiag_app_B_one(s->sp, s->p1);
               break;
#endif 
            case HYBRID_BETAX_PSI1_CALCULATE:
               betaxpsi1_calculate_one(s->sp, s->ion, s->nion, s->sintR, s->sintI, s->kpt, s->weiptr);
               break;
#if PAPI_PERFMON
            case HYBRID_FINALIZE_PAPI:
               PAPI_stop(EventSet, Papi_values);
               PTHREAD_PAPI_COUNTERS[s->tid] = Papi_values[0] + Papi_values[1];
               break;
#endif
            default:
               break;
        }

        // Let the main thread know that we are done
        sem_post(&thread_sem);

    }
}

#if PAPI_PERFMON
// Called at the end to retreive thread FLOP counts
long long Papi_thread_flops(int tid) {
    return PTHREAD_PAPI_COUNTERS[tid];
}

void Papi_init_omp_threads(int ithread) {

  pthread_t tid;

  ct.OpenMpEventSet[ithread] = PAPI_NULL;

  ct.OpenMpPthreadId[ithread] = pthread_self(); 

  if (PAPI_create_eventset(&ct.OpenMpEventSet[ithread]) != PAPI_OK) {
     error_handler ("Cannot create OpenMP event set.\n");
  }
  if (PAPI_add_event(ct.OpenMpEventSet[ithread], PAPI_DP_OPS) != PAPI_OK) {
     error_handler ("Cannot add OpenMP event.\n");
  }

  PAPI_start(ct.OpenMpEventSet[ithread]);

}

void Papi_finalize_omp_threads(int dummy) {

    int ithread;
    long long Papi_values[4];

    Papi_values[0] = 0;
    Papi_values[1] = 0;
    Papi_values[2] = 0;
    Papi_values[3] = 0;


    for(ithread = 0;ithread < ct.THREADS_PER_NODE;ithread++) {
        if(pthread_equal(ct.OpenMpPthreadId[ithread], pthread_self())) {

            PAPI_stop(ct.OpenMpEventSet[ithread], Papi_values);
            ct.OpenMpFlopCount[ithread] = Papi_values[0] + Papi_values[1];
            return;

        }        
    }

}

#endif

// Called when the main thread of execution is waiting for a set of threads to finish
void wait_for_threads(int jobs) {
    int idx;
    for(idx = 0;idx < jobs;idx++) {
        sem_wait(&thread_sem);
    }
}

// Wakes jobs sleeping threads starting from tid=0 and counting up
// jobs must be less than THREADS_PER_NODE
void wake_threads(int jobs) {

    int thread;
    
    if(jobs > ct.THREADS_PER_NODE) {
        // If this happens it is a bug
        printf("More jobs than available threads scheduled\n");   
        exit(0);
    }

    pthread_mutex_lock(&job_mutex);
    job_count = jobs;
    pthread_mutex_unlock(&job_mutex);

    for(thread = 0;thread < jobs;thread++) {
        sem_post(&thread_control[thread].sync);
    }
 
}

// Initialization function for barriers
void scf_barrier_init(int nthreads) {
    pthread_barrier_init(&scf_barrier, NULL, nthreads);
}

// Blocks all threads until nthreads specified in the init call have reached this point
void scf_barrier_wait(void) {
    if(!in_threaded_region) return;
    pthread_barrier_wait(&scf_barrier);
}

// Termination function
void scf_barrier_destroy(void) {
    pthread_barrier_destroy(&scf_barrier);
}



// Initializes the key and sets a value into it
void scf_tsd_init(void) {
 pthread_key_create(&scf_thread_control_key, NULL);
}

// Sets a value into the key
void scf_tsd_set_value(void *s) {
     pthread_setspecific(scf_thread_control_key, s);
}

// Deletes the key
void scf_tsd_delete(void) {
 pthread_key_delete(scf_thread_control_key);
}


// Reads the basetag from the thread specific data. Returns 0 if we are not in
// a parallel region
int get_thread_basetag(void) {

    SCF_THREAD_CONTROL *ss;
    if(!in_threaded_region) return 0;
    ss = (SCF_THREAD_CONTROL *)pthread_getspecific(scf_thread_control_key);
    if(!ss) return 0;

    return ss->sp->istate;

}

// Gets the threads control structure pointer
SCF_THREAD_CONTROL *get_thread_control(void) {
    SCF_THREAD_CONTROL *ss;
    if(!in_threaded_region) return NULL;
    ss = (SCF_THREAD_CONTROL *)pthread_getspecific(scf_thread_control_key);
    if(!ss) return NULL;
    return ss;
}

// Reads the tid from the thread specific data. Returns -1 if we are not in 
// a parallel region
int get_thread_tid(void) {

    SCF_THREAD_CONTROL *ss;

    if(!in_threaded_region) return -1;
    ss = (SCF_THREAD_CONTROL *)pthread_getspecific(scf_thread_control_key);
    if(!ss) return -1;

    return ss->tid;
}


// Reads the threads MPI grid communicator from the thread specific data
MPI_Comm *get_thread_grid_comm(void) {

    SCF_THREAD_CONTROL *ss;

    if(!in_threaded_region) return &pct.grid_comm;
    ss = (SCF_THREAD_CONTROL *)pthread_getspecific(scf_thread_control_key);
    if(!ss) return &pct.grid_comm;

    return &ss->grid_comm;
}


#if GPU_ENABLED
// Gets thread cstream
cudaStream_t *get_thread_cstream(void) {

    SCF_THREAD_CONTROL *ss;
    if(!in_threaded_region) return &ct.cuda_stream;  // Return main process stream
    ss = (SCF_THREAD_CONTROL *)pthread_getspecific(scf_thread_control_key);
    if(!ss) return NULL;
    return &ss->cstream;
}
#endif

// Reads the pinned memory storage from the thread specific data
void *get_thread_trade_buf(void) {

    SCF_THREAD_CONTROL *ss;

    if(!in_threaded_region) return NULL;
    ss = (SCF_THREAD_CONTROL *)pthread_getspecific(scf_thread_control_key);
    if(!ss) return NULL;
    return ss->trade_buf;

}

// Used for positioning and setting processor affinity. For now assumes that
// THREADS_PER_NODE is an even multiple of ct.ncpus. If this is not true it
// does not attemp to schedule
void set_cpu_affinity(int tid)
{
    int s;
    cpu_set_t cpuset;
    pthread_t thread;

    if(ct.THREADS_PER_NODE % ct.ncpus) return;

    s = tid % ct.THREADS_PER_NODE;

    
    // Set affinity mask
    CPU_ZERO(&cpuset);
    CPU_SET(s, &cpuset);

    thread = pthread_self();
    pthread_setaffinity_np(thread, sizeof(cpu_set_t), &cpuset);

}

void enter_threaded_region(void) {
    in_threaded_region = 1;
}
void leave_threaded_region(void) {
    in_threaded_region = 0;
}

void RMG_MPI_lock(void) {
    pthread_mutex_lock(&mpi_mutex);
}

void RMG_MPI_unlock(void) {
    pthread_mutex_unlock(&mpi_mutex);
}

void RMG_MPI_thread_order_lock(void) {
   int tid, i1, ntid;
   tid = get_thread_tid();

   if(tid < 0) {
       printf("\nError in RMG_MPI_thread_order_lock. Terminating.\n");
       MPI_Finalize();
       exit(0);
   }

   while(1) {

       // Acquire the lock
       pthread_mutex_lock(&thread_order_mutex);

       // See if it's our turn
       i1 = mpi_thread_order_counter % ct.THREADS_PER_NODE;
       if(i1 == tid) {
           // Raise priority of next thread
           ntid = i1 + 1;
           if(ntid < ct.THREADS_PER_NODE) {
               pthread_setschedprio(threads[ntid], -19);
           }
           return;
       }

       pthread_mutex_unlock(&thread_order_mutex);
       sched_yield();

   }

}

void RMG_MPI_thread_order_unlock(void) {
    
  int tid;
  tid = get_thread_tid();

  mpi_thread_order_counter++;
  pthread_setschedprio(threads[tid], -1);
  pthread_mutex_unlock(&thread_order_mutex);

}


extern REAL timings[LAST_TIME];
static pthread_mutex_t timings_mutex = PTHREAD_MUTEX_INITIALIZER;
void rmg_timings (int what, REAL time)
{
    pthread_mutex_lock(&timings_mutex);
    if(in_threaded_region) {
        timings[what] += time / ct.THREADS_PER_NODE;
    }
    else {
        timings[what] += time;
    }
    pthread_mutex_unlock(&timings_mutex);
}                               /* end rmg_timings */



#endif

// Tells us if it's safe to use OMP in a region
int rmg_is_open_mp_safe(void)
{

#if !HYBRID_MODEL
    return 1;
#else
    if(in_threaded_region) return 0;
    return 1;
#endif
}

// Tells us if we are executing a parallel region that is a loop over orbitals
int is_loop_over_states(void)
{

#if !HYBRID_MODEL
    return 1;
#else
    SCF_THREAD_CONTROL *ss;
    if(!in_threaded_region) return 0;
    ss = (SCF_THREAD_CONTROL *)pthread_getspecific(scf_thread_control_key);
    if(!ss) return 0;
    if((ss->job == HYBRID_EIG) || (ss->job == HYBRID_SUBDIAG_APP_A) || (ss->job == HYBRID_SUBDIAG_APP_B))
        return 1;

    return 0;
#endif
    
}


// Some dummy function stubs so we don't have to put in #if blocks all over the place
#if !HYBRID_MODEL
void RMG_MPI_thread_order_lock(void) {
}
void RMG_MPI_thread_order_unlock(void) {
}
void scf_barrier_wait(void) {
}
void scf_barrier_init(int nthreads) {
}
void scf_barrier_destroy(void) {
}
#endif


#if RMG_OMP_THREADS
// Some of our routines can be called from serial or parallel regions. To avoid
// oversubscription of resources when using OMP in these routines use the following
// set of functions.
void RMG_set_omp_parallel(void)
{
    if(rmg_is_open_mp_safe()) {
        omp_set_num_threads(RMG_OMP_THREADS);
    }
    else {
        omp_set_num_threads(1);
    }
}

void RMG_set_omp_single(void)
{

}
#endif
