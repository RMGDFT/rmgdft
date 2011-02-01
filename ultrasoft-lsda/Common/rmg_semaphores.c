/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

/****f* QMD-MGDFT/rmg_semaphores.c *****
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
 *   Semapore routines for systems that don't have them (AIX).
 * SOURCE
 */

#include	"main.h"


#if SMP
#if (AIX || LINUX || XT3)

void QMD_sem_init (QMD_sem_t * sem)
{
    pthread_mutex_init (&sem->lock, NULL);
    pthread_cond_init (&sem->cond, NULL);
    sem->count = 0;
}

void QMD_sem_destroy (QMD_sem_t * sem)
{
    pthread_mutex_destroy (&sem->lock);
    pthread_cond_destroy (&sem->cond);
}

void p_operation_cleanup (void *arg)
{
    QMD_sem_t *sem;
    sem = (QMD_sem_t *) arg;
    pthread_mutex_unlock (&sem->lock);
}

void QMD_sem_wait (QMD_sem_t * sem)
{
    pthread_mutex_lock (&sem->lock);
    pthread_cleanup_push (p_operation_cleanup, sem);
    while (sem->count == 0)
        pthread_cond_wait (&sem->cond, &sem->lock);
    sem->count++;
    /*
     *  Note that the pthread_cleanup_pop subroutine will
     *  execute the p_operation_cleanup routine
     */
    pthread_cleanup_pop (1);
}


void QMD_sem_post (QMD_sem_t * sem)
{
    pthread_mutex_lock (&sem->lock);
    sem->count--;
    if (sem->count <= 0)
        pthread_cond_signal (&sem->cond);
    pthread_mutex_unlock (&sem->lock);
}
#else



void QMD_sem_init (QMD_sem_t * sem)
{
    sem_init (&sem->sem, 0, 0);
}

void QMD_sem_destroy (QMD_sem_t * sem)
{
    sem_destroy (&sem->sem);
}

void QMD_sem_wait (QMD_sem_t * sem)
{
    sem_wait (&sem->sem);
}

void QMD_sem_post (QMD_sem_t * sem)
{
    sem_post (&sem->sem);
}

#endif

#endif

/******/
