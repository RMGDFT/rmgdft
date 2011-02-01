/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

/****f* QMD-MGDFT/threads.c *****
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
 *   Thread utility routines for SMP and NUMA machines.
 *   void create_threads(STATE *states)
 *   void start_threads(int action)
 *   void thread_dispatch(SCF_THREAD_CONTROL *s, int job)
 *   void wait_for_threads(void)
 * SOURCE
 */


#include "main.h"

#if SMP
#include <errno.h>


#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <math.h>
#include <stdio.h>



/* This is used to assign orbitals to threads. It is set to zero before a multigrid
 * sweep is started. Then as each state takes a thread it increments the count.
 */
volatile int state_queue;
pthread_mutex_t lock_mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t lock1_mutex = PTHREAD_MUTEX_INITIALIZER;

volatile int done_flag;
pthread_mutex_t lock2_mutex = PTHREAD_MUTEX_INITIALIZER;

/* If the pthreads package does not allow us to explicitly
 * control concurrency we do it ourselves. */
pthread_mutex_t concurrency_mutex = PTHREAD_MUTEX_INITIALIZER;
volatile int concurrency_count;


SCF_THREAD_CONTROL thread_control[MAX_THREADS + 1];
static pthread_t scf_threads[MAX_THREADS + 1];
static pthread_t main_thread;
QMD_sem_t main_sem[MAX_THREADS];
QMD_sem_t thread_sem;



static pthread_attr_t thread_attrs;



void create_threads (STATE * states)
{

    int thread, size, sizr;
    int offset;
    char tbuf[100];
    int maxm, retval;


    /* First set up the threads for the multigrid */
    pthread_attr_init (&thread_attrs);

    /* set contention scope to system for SMPs */
    maxm = ct.thread_concurrency / 2;

#if IRIX
    if (!maxm)
        maxm++;
    sprintf (tbuf, "threads %d", ct.num_threads + 1);
    dplace_line (tbuf);
    /* Fix block or cyclic distribution */
    sprintf (tbuf, "memories %d in topology cluster", maxm);
    dplace_line (tbuf);
    sprintf (tbuf, "distribute threads 1:%d across memories", ct.num_threads);
    dplace_line (tbuf);
#endif

    for (thread = 0; thread < ct.num_threads; thread++)
    {

        thread_control[thread].my_states = states;
        thread_control[thread].start = 0;
        thread_control[thread].tid = thread;
        QMD_sem_init (&main_sem[thread]);

    }                           /* end for */


    QMD_sem_init (&thread_sem);


    /* Now we set up the threads for orthogonalization */
    size = ct.psi_nbasis / ct.num_threads;
    sizr = ct.psi_nbasis % ct.num_threads;

    for (thread = 0; thread < ct.num_threads; thread++)
    {

        thread_control[thread].numpt = size;

    }                           /* end for */


    for (thread = 0; thread < sizr; thread++)
    {

        thread_control[thread].numpt++;

    }                           /* end for */


    /* Set up offsets for gather_psi */
    offset = 0;
    for (thread = 0; thread < ct.num_threads; thread++)
    {

        thread_control[thread].offset = offset;
        offset += thread_control[thread].numpt;

    }                           /* end for */


    for (thread = 0; thread < ct.num_threads; thread++)
    {

        retval = pthread_create (&scf_threads[thread],
                                 &thread_attrs, thread_scheduler, &thread_control[thread]);
        if (retval)
        {
            printf ("ERROR = %d\n", retval);
            error_handler ("Unable to create threads");
        }
    }                           /* end for */


    wait_for_threads ();


    retval = pthread_create (&main_thread, &thread_attrs, run, NULL);
    if (retval)
    {
        printf ("ERROR = %d\n", retval);
        error_handler ("Unable to create threads");
    }


#if (IRIX || AIX)
    pthread_setconcurrency (ct.thread_concurrency);
#endif
}                               /* end create_threads */





/* Starts the threads up. Eventually action should be an array so we can
 * assign different tasks to different threads.
 */
void start_threads (int action)
{

    int thread, st1;
    state_queue = 0;
    done_flag = 0;

    /* Set to zero */
    concurrency_count = 0;


    /* Assign the action */
    for (thread = 0; thread < ct.num_threads; thread++)
        thread_control[thread].start = action;


    /* Wake them up */
    for (thread = 0; thread < ct.num_threads; thread++)
    {

        if (action != SMP_ORTHO1)
        {
            while (concurrency_count >= ct.thread_concurrency)
                usleep (10000);
            pthread_mutex_lock (&concurrency_mutex);
            concurrency_count++;
            pthread_mutex_unlock (&concurrency_mutex);
            QMD_sem_post (&main_sem[thread]);
        }
        else
        {
            QMD_sem_post (&main_sem[thread]);
        }                       /* end if */

    }

}                               /* end start_threads */



/* This is called by the main thread when it is waiting for the auxiliary
 * threads to finish processing. I suppose we could do something a bit more
 * clever here for the synchronization but this works fine since the
 * auxiliary threads are typically working for periods of time on the order of
 * a few seconds.
 */
void wait_for_threads (void)
{

    QMD_sem_wait (&thread_sem);

}                               /* end wait_for_threads */





/* Thread scheduling routine */
void *thread_scheduler (void *set)
{

    int job, cnum, clines, cpad;
    SCF_THREAD_CONTROL *s;

    s = (SCF_THREAD_CONTROL *) set;


    /* The offset of is to prevent cache line thrashing in the ortho
     * routines. The number probably needs to be tuned for different machines
     * but 112 seems to work quite well on an Origin 2000 with R10K processors.
     */
    cnum = CACHE_LINE_SIZE / sizeof (REAL);
    clines = s->numpt / cnum;
    clines = 11 - (clines % 11);
    cpad = clines * cnum;

    s->lda = s->numpt + cpad;


    /* We allocate our memory here when the thread starts since on the O2K
     * this insures that the memory is local to the node on which the thread
     * is running.
     */
    my_calloc (s->rho, s->lda, REAL);
#if GAMMA_PT
    my_calloc (s->base_mem, (ct.num_states + 1) * s->lda + 112, REAL);
    my_calloc (s->darr, ct.num_states * ct.num_states, REAL);
#else
    my_calloc (s->base_mem, ct.num_kpts * 2 * (ct.num_states + 1) * s->lda + 112, REAL);
    my_calloc (s->darr, 2 * ct.num_states * ct.num_states, REAL);
#endif
    my_calloc (s->rhocore, s->lda, REAL);
    my_calloc (s->scratch1, 2 * s->lda, REAL);

#if IRIX
    my_calloc (s->vtot, P0_BASIS, REAL);
#else
    if (!s->tid)
        my_calloc (s->vtot, P0_BASIS, REAL);
#endif

    my_malloc (s->eigs, ct.num_states, REAL);
    my_malloc (s->occs, ct.num_states, REAL);

    /* Initialization is done so wake up the main thread */
    pthread_mutex_lock (&lock1_mutex);
    done_flag++;
    if (done_flag == ct.num_threads)
        QMD_sem_post (&thread_sem);
    pthread_mutex_unlock (&lock1_mutex);


    /* Each thread sits in this loop waiting to be woken up by the main
     * thread. The main thread tells this thread what to do by setting the
     * start field in the thread control structure to a particular action
     * code.
     */
    while (1)
    {

        QMD_sem_wait (&main_sem[s->tid]);

        job = s->start;
        s->start = 0;

        thread_dispatch (s, job);

        pthread_mutex_lock (&lock1_mutex);
        done_flag++;
        if (done_flag == ct.num_threads)
            QMD_sem_post (&thread_sem);
        pthread_mutex_unlock (&lock1_mutex);

    }                           /* end while */


    return NULL;

}                               /* end thread_scheduler */




/* Thread dispatch routine */
void thread_dispatch (SCF_THREAD_CONTROL * s, int job)
{

    int idx, threadb, offset, size;
    int thread, st1, ione = 1;

    switch (job)
    {

    case SMP_EIG:
        /* Process the orbitals that are assigned to me */
#if IRIX

        offset = 0;
        for (thread = 0; thread < ct.num_threads; thread++)
        {

            threadb = s->tid + thread;
            if (threadb > (ct.num_threads - 1))
                threadb -= ct.num_threads;
            offset = thread_control[threadb].offset;
            size = thread_control[threadb].numpt;
            QMD_scopy (size, thread_control[threadb].scratch1, ione,
                       &thread_control[s->tid].vtot[offset], ione);
            offset += size;

        }                       /* end for */

#else

        s->vtot = thread_control[0].vtot;

#endif

        /* Each thread increments state_queue to indicate that it is
         * processing a particular thread.
         */
        while (1)
        {

            pthread_mutex_lock (&lock_mutex);
            if (state_queue == ct.num_kpts * ct.num_states)
            {
                pthread_mutex_unlock (&lock_mutex);
                break;
            }                   /* end if */
            st1 = state_queue;
            state_queue++;
            pthread_mutex_unlock (&lock_mutex);
            mg_eig_state (&s->my_states[st1], s->tid, s->vtot);

        }                       /* end while */

        break;

    case SMP_ORTHO1:
        ortho1_smp (s);
        break;

    case SMP_ORTHO2:
        ortho2_smp (s);
        break;

    case SMP_GET_RHO:
        get_rho_smp (s);
        break;

    case SMP_SORT_PSI:
        sort_psi_smp (s);
        break;

    case SMP_DIAG1:
        subdiag1_smp (s);
        break;

    case SMP_DIAG2:
        subdiag2_smp (s);
        break;

    case SMP_NLFORCE:
        /* Zero out the force component array */
        for (idx = 0; idx < ct.num_ions; idx++)
        {

            s->force[idx][0] = 0.0;
            s->force[idx][1] = 0.0;
            s->force[idx][2] = 0.0;

        }                       /* end for */

        nlforce (s->my_states, s->tid);

        break;

    case SMP_GETNLOP:
        get_nlop_smp (s->tid);
        break;

    case SMP_SKIP:
        /* This is here in case the main thread wishes to assign tasks to
         * some threads but not all. In that case it sets s->start to skip
         * for the threads that don't need to do anything.
         */
        break;

    default:
        break;

    }                           /* end switch */

    if (job != SMP_ORTHO1)
    {
        pthread_mutex_lock (&concurrency_mutex);
        concurrency_count--;
        pthread_mutex_unlock (&concurrency_mutex);
    }                           /* end if */

}                               /* end thread_dispatch */


void QMD_thread_barrier (QMD_thread_barrier_struct * bs)
{

    pthread_mutex_lock (&bs->mutex);
    bs->count++;

    if (bs->count == ct.num_threads)
    {

        pthread_cond_broadcast (&bs->cond);

    }
    else
    {

        while (bs->count != ct.num_threads)
        {

            pthread_cond_wait (&bs->cond, &bs->mutex);

        }                       /* end while */

    }                           /* end if */

    pthread_mutex_unlock (&bs->mutex);

}                               /* end QMD_thread_barrier */



#endif
/******/
