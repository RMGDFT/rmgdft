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
 *   rmg_double_t real_sum_all(rmg_double_t x)
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


#include <hybrid.h>
#include <pthread.h>
volatile rmg_double_t real_sum_all_vector[MAX_SCF_THREADS];
volatile  rmg_double_t recvbuf[MAX_SCF_THREADS];
volatile int real_sum_all_vector_state = 0;
pthread_mutex_t real_sum_all_vector_lock = PTHREAD_MUTEX_INITIALIZER;
static rmg_double_t real_sum_all_threaded(rmg_double_t x, int tid, MPI_Comm comm);



rmg_double_t real_sum_all (rmg_double_t x, MPI_Comm comm)
{

    rmg_double_t inreg;
    rmg_double_t outreg;
    int tid;
#if MD_TIMERS
    rmg_double_t time0;

    time0 = my_crtc ();
#endif
	
    // If we are not in a pthreads parallel region get_thread_tid will return a negative
    // value so in that case we just fall through to the regular routine.
    tid = get_thread_tid();
    if(tid >= 0) {
        return real_sum_all_threaded(x, tid, comm);
    }

    inreg = x;

    MPI_Allreduce (&inreg, &outreg, 1, MPI_DOUBLE, MPI_SUM, comm);

#if MD_TIMERS
    rmg_timings (REAL_SUM_ALL_TIME, my_crtc () - time0);
#endif


    return outreg;



}                               /* end real_sum_all */



// Used to sum a block of data from a set of threads operating in parallel
rmg_double_t real_sum_all_threaded(rmg_double_t x, int tid,  MPI_Comm comm) {

#if MD_TIMERS
  rmg_double_t time0;
  time0 = my_crtc ();
#endif

  // First load the data in the array. If real_sum_all_vector_state is 0 set it to 1
  pthread_mutex_lock(&real_sum_all_vector_lock);
      real_sum_all_vector_state = 1;
      real_sum_all_vector[tid] = x;
  pthread_mutex_unlock(&real_sum_all_vector_lock);


  // Wait until everyone gets here
  thread_barrier_wait();

  // Data is all loaded now and we only want one thread to do the MPI call
  // Might have some contention here for high core counts
  pthread_mutex_lock(&real_sum_all_vector_lock);
      if(real_sum_all_vector_state == 1) {
          MPI_Allreduce(real_sum_all_vector, recvbuf, ct.THREADS_PER_NODE, MPI_DOUBLE, MPI_SUM, comm);
          real_sum_all_vector_state = 0;
      }
  pthread_mutex_unlock(&real_sum_all_vector_lock);

#if MD_TIMERS
   rmg_timings (REAL_SUM_ALL_TIME, my_crtc () - time0);
#endif

  return recvbuf[tid];

}

