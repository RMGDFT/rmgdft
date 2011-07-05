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

#if HYBRID_MODEL

// Main thread control structure
SCF_THREAD_CONTROL thread_control[THREADS_PER_NODE];

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
static pthread_t threads[THREADS_PER_NODE];
volatile int in_threaded_region = 0;
static void run_threads(SCF_THREAD_CONTROL *s);

// These are used to ensure thread ordering
volatile int mpi_thread_order_counter = 0;
static pthread_mutex_t thread_order_mutex = PTHREAD_MUTEX_INITIALIZER;

// Used for accessing thread specific data
pthread_key_t scf_thread_control_key;

// Initialization function called by main process
void init_HYBRID_MODEL(void) {

    int thread, retval;

    sem_init(&thread_sem, 0, 0);
    ct.main_thread_pid = getpid();

    pthread_attr_init( &thread_attrs );
    pthread_attr_setschedpolicy( &thread_attrs, SCHED_RR);

    // Create the thread specific data key
    scf_tsd_init();

    // Create the main sync barrier
    pthread_barrier_init(&run_barrier, NULL, THREADS_PER_NODE);

    // Here we create a set of long lived threads
    for(thread = 0;thread < THREADS_PER_NODE;thread++) {

        thread_control[thread].tid = thread;
        sem_init(&thread_control[thread].sync, 0, 0);
        retval = pthread_create(&threads[thread], &thread_attrs, (void *)run_threads, (void *)&thread_control[thread]);

    }

    
}

// Main thread function
void run_threads(SCF_THREAD_CONTROL *s) {

    s->pthread_tid = pthread_self();
    pthread_setspecific(scf_thread_control_key, (void *)s);

    // Wait until everyone gets here
    pthread_barrier_wait(&run_barrier);

    while(1) {

        // We sleep forever or until we get a signal that wakes us up
        sem_wait(&s->sync);

        // Switch that controls what we do
        switch(s->job) {
            case HYBRID_EIG:       // Performs a single multigrid sweep over an orbital
               mg_eig_state(s->sp, 0, s->vtot);
               break;
            case HYBRID_SKIP:
               break;
            case HYBRID_SUBDIAG_APP_A:
               subdiag_app_A_one(s->sp, s->p1, s->p2, s->vtot);
               break;
            case HYBRID_SUBDIAG_APP_B:
               subdiag_app_B_one(s->sp, s->p1);
               break;
            case HYBRID_BETAX_PSI1_CALCULATE:
               betaxpsi1_calculate_one(s->sp, s->ion, s->nion, s->sintR, s->sintI);
               break;
            default:
               break;
        }

        // Let the main thread know that we are done
        sem_post(&thread_sem);

    }
}

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
    
    if(jobs > THREADS_PER_NODE) {
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

// Reads the tid from the thread specific data. Returns -1 if we are not in 
// a parallel region
int get_thread_tid(void) {

    SCF_THREAD_CONTROL *ss;

    if(!in_threaded_region) return -1;
    ss = (SCF_THREAD_CONTROL *)pthread_getspecific(scf_thread_control_key);
    if(!ss) return -1;

    return ss->tid;
}


// Used for positioning and setting processor affinity
void set_cpu_affinity(void)
{
    int s, j;
    cpu_set_t cpuset;
    pthread_t thread;

    thread = pthread_self();

    // Set affinity mask to include CPUs 0 up to 12 which is good for dual Istanbul

    CPU_ZERO(&cpuset);
    for(j=0;j<12;j++) CPU_SET(j, &cpuset);

    s = pthread_setaffinity_np(thread, sizeof(cpu_set_t), &cpuset);

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
       i1 = mpi_thread_order_counter % THREADS_PER_NODE;
       if(i1 == tid) {
           // Raise priority of next thread
           ntid = i1 + 1;
           if(ntid < THREADS_PER_NODE) {
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
        timings[what] += time / THREADS_PER_NODE;
    }
    else {
        timings[what] += time;
    }
    pthread_mutex_unlock(&timings_mutex);
}                               /* end rmg_timings */



#endif

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
#endif
