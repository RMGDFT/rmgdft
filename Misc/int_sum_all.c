/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

/****f* QMD-MGDFT/real_sum_all.c *****
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
 * FUNCTION
 *   REAL real_sum_all(REAL x)
 *   Performs a scalar sum over all processors.
 * INPUTS
 *   x: defined in each processor
 * OUTPUT
 *   sum over all processors is returned to each processor
 * PARENTS
 *   get_ke.c get_rho.c get_te.c get_vh.c getpoi_bc.c init_nuc.c
 *   lforce.c mg_eig_state.c norm_psi.c 
 * CHILDREN
 *   nothing
 * SOURCE
 */


#include <float.h>
#include <math.h>
#include "main.h"


#if HYBRID_MODEL
    #include <hybrid.h>
    #include <pthread.h>
    volatile int int_sum_all_vector[THREADS_PER_NODE];
    volatile int recvbuf[THREADS_PER_NODE];
    volatile int int_sum_all_vector_state = 0;
    pthread_mutex_t int_sum_all_vector_lock = PTHREAD_MUTEX_INITIALIZER;
    static int int_sum_all_threaded(int x, int tid, MPI_Comm comm);
#endif



int int_sum_all (int x, MPI_Comm comm)
{

    int inreg;
    int outreg;
    int tid;
#if MD_TIMERS
    REAL time0;

    time0 = my_crtc ();
#endif
	
#if HYBRID_MODEL
    // If we are not in a pthreads parallel region get_thread_tid will return a negative
    // value so in that case we just fall through to the regular routine.
    tid = get_thread_tid();
    if(tid >= 0) {
        return int_sum_all_threaded(x, tid, comm);
    }
#endif

    inreg = x;

    MPI_Allreduce (&inreg, &outreg, 1, MPI_INT, MPI_SUM, comm);

#if MD_TIMERS
    rmg_timings (REAL_SUM_ALL_TIME, my_crtc () - time0);
#endif


    return outreg;



}                               /* end real_sum_all */


#if HYBRID_MODEL

// Used to sum a block of data from a set of threads operating in parallel
int int_sum_all_threaded(int x, int tid,  MPI_Comm comm) {

#if MD_TIMERS
  REAL time0;
  time0 = my_crtc ();
#endif

  // First load the data in the array. If int_sum_all_vector_state is 0 set it to 1
  pthread_mutex_lock(&int_sum_all_vector_lock);
      int_sum_all_vector_state = 1;
      int_sum_all_vector[tid] = x;
  pthread_mutex_unlock(&int_sum_all_vector_lock);


  // Wait until everyone gets here
  scf_barrier_wait();

  // Data is all loaded now and we only want one thread to do the MPI call
  // Might have some contention here for high core counts
  pthread_mutex_lock(&int_sum_all_vector_lock);
      if(int_sum_all_vector_state == 1) {
          MPI_Allreduce(int_sum_all_vector, recvbuf, THREADS_PER_NODE, MPI_INT, MPI_SUM, comm);
          int_sum_all_vector_state = 0;
      }
  pthread_mutex_unlock(&int_sum_all_vector_lock);

#if MD_TIMERS
   rmg_timings (REAL_SUM_ALL_TIME, my_crtc () - time0);
#endif

  return recvbuf[tid];

}

#endif
