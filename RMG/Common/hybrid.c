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
#include "main.h"
#include "hybrid.h"

#if HYBRID_MODEL
void ATL_assert(void)
{
}

// Used to implement a local barrier inside of the scf loops
pthread_barrier_t scf_barrier;

volatile int in_threaded_region = 0;

// Initialization function
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



// Used for accessing thread specific data
pthread_key_t scf_thread_control_key;

// Initializes the key and sets a value into it
void scf_tsd_init(void) {
 pthread_key_create(&scf_thread_control_key, NULL);
 in_threaded_region = 1;
}

// Sets a value into the key
void scf_tsd_set_value(void *s) {
     pthread_setspecific(scf_thread_control_key, s);
}

// Deletes the key
void scf_tsd_delete(void) {
 pthread_key_delete(scf_thread_control_key);
 in_threaded_region = 0;
}


// Reads the basetag from the thread specific data. Returns 0 if we are not in
// a parallel region
int get_thread_basetag(void) {

    MG_THREAD_STRUCT *ss;
    if(!in_threaded_region) return 0;
    ss = (MG_THREAD_STRUCT *)pthread_getspecific(scf_thread_control_key);
    if(!ss) return 0;

    return ss->sp->istate;

}

// Reads the tid from the thread specific data. Returns -1 if we are not in 
// a parallel region
int get_thread_tid(void) {

    MG_THREAD_STRUCT *ss;

    if(!in_threaded_region) return -1;
    ss = (MG_THREAD_STRUCT *)pthread_getspecific(scf_thread_control_key);
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
#endif
