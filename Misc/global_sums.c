/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

/****f* QMD-MGDFT/global_sums.c *****
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
 *   void global_sums (rmg_double_t *vect, int *length)
 *   Sums an array over all processors. For serial machines it just returns.
 * INPUTS
 *   vect: a vector to be sumed. Each processor has his own value
 *   length: length of the vector
 * OUTPUT
 *   vect: Each processor gets the sum.
 * PARENTS
 *   app_nl.c mg_eig_state.c nlforce_d.c nlforce_p.c
 *   nlforce_s.c rft.c scf.c subdiag_mpi.c symmetry.c write_avgd.c 
 *   write_avgv.c write_zstates.c
 * CHILDREN
 *   Mpi_allreduce is a MPI routine
 * SOURCE
 */


#include "main.h"
#include "common_prototypes.h"

static rmg_double_t *fixed_vector1 = NULL;
static rmg_double_t *fixed_vector2 = NULL;
#define MAX_FIXED_VECTOR 512

#if MPI

#if HYBRID_MODEL

#include <hybrid.h>
#include <pthread.h>
volatile rmg_double_t *global_sums_vector, *tvector;
volatile int global_sums_vector_state = 0;
pthread_mutex_t global_sums_vector_lock = PTHREAD_MUTEX_INITIALIZER;
static void global_sums_threaded (rmg_double_t *vect, int *length, int tid, MPI_Comm comm);



void global_sums_threaded (rmg_double_t *vect, int *length, int tid, MPI_Comm comm)
{

  rmg_double_t *rptr, *rptr1;

  if(*length < MAX_FIXED_VECTOR) {

      QMD_dcopy(*length, vect, 1, &fixed_vector1[*length * tid], 1);
      scf_barrier_wait();
      if(tid == 0)
          MPI_Allreduce(fixed_vector1, fixed_vector2, *length * ct.THREADS_PER_NODE, MPI_DOUBLE, MPI_SUM, comm);
      scf_barrier_wait();
      QMD_dcopy(*length,  &fixed_vector2[*length * tid], 1, vect, 1);

      return;
  }

  scf_barrier_wait();
  pthread_mutex_lock(&global_sums_vector_lock);
      if(global_sums_vector_state == 0) {
          my_malloc (global_sums_vector, *length * ct.THREADS_PER_NODE, rmg_double_t);
          my_malloc (tvector, *length * ct.THREADS_PER_NODE, rmg_double_t);
      }
      global_sums_vector_state = 1;
  pthread_mutex_unlock(&global_sums_vector_lock);

  QMD_dcopy(*length, vect, 1, &global_sums_vector[*length * tid], 1);
  
  // Wait until everyone gets here
  scf_barrier_wait();

  pthread_mutex_lock(&global_sums_vector_lock);
      if(global_sums_vector_state == 1) {
          MPI_Allreduce(global_sums_vector, tvector, *length * ct.THREADS_PER_NODE, MPI_DOUBLE, MPI_SUM, comm);
          global_sums_vector_state = 0;
      }
  pthread_mutex_unlock(&global_sums_vector_lock);
  QMD_dcopy(*length,  &tvector[*length * tid], 1, vect, 1);

  // Must wait until all threads have copied the data to vect before freeing memory
  scf_barrier_wait();

  pthread_mutex_lock(&global_sums_vector_lock);
      // ensures that the memory is only freed once
      if(tvector != NULL) {
          my_free(tvector);
          my_free(global_sums_vector);
          tvector = NULL;
          global_sums_vector = NULL;
      }
  pthread_mutex_unlock(&global_sums_vector_lock);

}

#endif


void init_global_sums(void) {
    int retval;
    retval = MPI_Alloc_mem(sizeof(rmg_double_t) * ct.THREADS_PER_NODE * MAX_FIXED_VECTOR , MPI_INFO_NULL, &fixed_vector1);
    if(retval != MPI_SUCCESS) {
        error_handler("Error in MPI_Alloc_mem.\n");
    }
    retval = MPI_Alloc_mem(sizeof(rmg_double_t) * ct.THREADS_PER_NODE * MAX_FIXED_VECTOR , MPI_INFO_NULL, &fixed_vector2);
    if(retval != MPI_SUCCESS) {
        error_handler("Error in MPI_Alloc_mem.\n");
    }
}


void global_sums (rmg_double_t * vect, int *length, MPI_Comm comm)
{
    int sizr, steps, blocks, newsize, tid;
    rmg_double_t *rptr, *rptr1;
    rmg_double_t rptr2[100];
#if MD_TIMERS
    rmg_double_t time0;

    time0 = my_crtc ();
#endif

#if HYBRID_MODEL
    tid = get_thread_tid();
    if(tid >= 0) {
        global_sums_threaded(vect, length, tid, comm);
        return;
    }
#endif

    /* Check for small vector case and handle on stack */
    if (*length < 100)
    {
        sizr = *length;
        QMD_dcopy (sizr, vect, 1, rptr2, 1);
       	MPI_Allreduce (rptr2, vect, sizr, MPI_DOUBLE, MPI_SUM, comm);

		
#if MD_TIMERS
        rmg_timings (GLOBAL_SUMS_TIME, my_crtc () - time0);
#endif
        return;
    }

    my_malloc (rptr, MAX_PWRK, rmg_double_t);
    newsize = MAX_PWRK;
    blocks = *length / newsize;
    sizr = (*length % newsize);

    rptr1 = vect;

    for (steps = 0; steps < blocks; steps++)
    {
        QMD_dcopy (newsize, rptr1, 1, rptr, 1);
        MPI_Allreduce (rptr, rptr1, newsize, MPI_DOUBLE, MPI_SUM, comm);  

        rptr1 += newsize;
    }

    if (sizr)
    {
        QMD_dcopy (sizr, rptr1, 1, rptr, 1);
        MPI_Allreduce (rptr, rptr1, sizr, MPI_DOUBLE, MPI_SUM, comm);
    }

    my_free (rptr);

#if MD_TIMERS
    rmg_timings (GLOBAL_SUMS_TIME, my_crtc () - time0);
#endif
}                               /* end global_sums */



#else



void global_sums (rmg_double_t * vect, int *length, MPI_Comm comm)
{
    return;
}

#endif


/******/
