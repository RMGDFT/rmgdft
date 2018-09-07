/************************** SVN Revision Information **************************
 **    $Id: orbit_dot_orbit.c 4341 2018-04-12 01:20:36Z emil $    **
 ******************************************************************************/


#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>


void DistributeTasks(int active_threads, int num_tasks, int *task_start, int *task_end)
{

    int tasks_per_thread = num_tasks/active_threads;
    int tasks_remain = num_tasks % active_threads;


    task_start[0] = 0;

    for(int ist = 0;ist < active_threads-1;ist++) {
        task_end[ist] = task_start[ist]  + tasks_per_thread;
        if (ist < tasks_remain) task_end[ist] += 1;
        task_start[ist+1] = task_end[ist];
    }
    task_end[active_threads-1] = num_tasks;
}

