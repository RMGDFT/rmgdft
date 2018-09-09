/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "main.h"
#include "prototypes_on.h"
#include "init_var.h"
#include "BaseThread.h"
#include "rmgthreads.h"
#include "RmgThread.h"

void precond(double *x)
{

    BaseThread *T = BaseThread::getBaseThread(0);
    int active_threads = ct.MG_THREADS_PER_NODE;

    SCF_THREAD_CONTROL thread_control;

    double *work1, *work2;
    double *sp;
    double gamma;
    int size;
    int idx;
    int istate;

    int nx = ct.max_orbit_nx;

    my_malloc_init( work1, 4 * ct.max_orbit_size, double );
    my_malloc_init( work2, 4 * ct.max_orbit_size, double );

    gamma = get_gamma_precond(vtot_c, states[0].eig[0]);


    int num_orbitals_thispe = ct.state_end - ct.state_begin;

    if( num_orbitals_thispe < active_threads) active_threads = num_orbitals_thispe;

    int task_start[active_threads];
    int task_end[active_threads];

    DistributeTasks(active_threads, num_orbitals_thispe, task_start, task_end);

    for(int ist = 0;ist < active_threads;ist++) {

        thread_control.job = HYBRID_ORBITALS_DOT_PRODUCT;
        thread_control.sp = (void *)states1;
        thread_control.p1 = (void *)states;
        thread_control.p2 = (void *)states;

        thread_control.nv = (void *)Aij;
        thread_control.ns = (void *)Bij;

        thread_control.basetag = task_start[ist];
        thread_control.extratag1 = active_threads;
        thread_control.extratag2 = task_end[ist];

        QueueThreadTask(ist, thread_control);

        //DotProductOrbitOrbit(&states1[st1], &states[st2], &states[st1],  H, S, onepair);
    }

    // Thread tasks are set up so run them
    T->run_thread_tasks(active_threads);





    size = 0;
    for (istate = ct.state_begin; istate < ct.state_end; istate++)
    {
        sp = x + size;
        size += states[istate].size;
        int ixx = states[istate].ixmax - states[istate].ixmin + 1;
        int iyy = states[istate].iymax - states[istate].iymin + 1;
        int izz = states[istate].izmax - states[istate].izmin + 1;

        for (idx = 0; idx < states[istate].size; idx++)
        {
            work1[idx] = gamma * sp[idx];
        }
        /* compute the preconditioned steepest descent direction
         * -> work1 */

//        precond_mg_c(sp, work1, work2, istate);
        precond_mg(sp, work1, work2, istate);
          app_mask(istate, work1, 0);
//          ZeroBoundaryC(work1, ixx, iyy, izz);


        for (idx = 0; idx < states[istate].size; idx++)
        {
            sp[idx] = 0.5 * work1[idx];
        }

    }

    my_free(work1);
    my_free(work2);

}

