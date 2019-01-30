/************************** SVN Revision Information **************************
 **    $Id: orbit_dot_orbit.c 4341 2018-04-12 01:20:36Z emil $    **
 ******************************************************************************/

/*


   orbit_dot_orbit.c

   work_matrix_row(i,j) = <states[i].psiR in this pe | states1[j].psiR>


 */

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
#include "transition.h"

extern std::vector<ORBITAL_PAIR> OrbitalPairs;

void OrbitDotOrbit(STATE * states, STATE * states1, double *Aij, double *Bij)
{

    //BaseThread *T = BaseThread::getBaseThread(0);
    int active_threads = ct.MG_THREADS_PER_NODE;

    if( OrbitalPairs.size() == 0) return;

//    printf("\n %d num pairs ", (int)OrbitalPairs.size());

    if( (int)OrbitalPairs.size() < active_threads) active_threads = (int)OrbitalPairs.size();

    int pair_start[active_threads];
    int pair_end[active_threads];

    DistributeTasks(active_threads, (int)OrbitalPairs.size(), pair_start, pair_end);
    int ist;
#pragma omp parallel private(ist)
    {
#pragma omp for schedule(static,1) nowait
        for(int ist = 0;ist < active_threads;ist++) {
            OrbitDotOrbitBlock(pair_start[ist], pair_end[ist], Aij, Bij); 
        }


    }
#if 0

    SCF_THREAD_CONTROL thread_control;
    //for(int ist = 0;ist < active_threads;ist++) {
    for(int ist = 0;ist < (int)OrbitalPairs.size();ist++) {

        thread_control.job = HYBRID_ORBITALS_DOT_PRODUCT;
        thread_control.sp = (void *)states1;
        thread_control.p1 = (void *)states;
        thread_control.p2 = (void *)states;

        thread_control.nv = (void *)Aij;
        thread_control.ns = (void *)Bij;

        thread_control.basetag = ist;
        thread_control.extratag1 = active_threads;
        thread_control.extratag2 = ist+1;

        QueueThreadTask(ist%active_threads, thread_control);

        //DotProductOrbitOrbit(&states1[st1], &states[st2], &states[st1],  H, S, onepair);
    }

    // Thread tasks are set up so run them
    T->run_thread_tasks(active_threads, Rmg_Q);
#endif


}

