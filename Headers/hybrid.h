/*
 *
 *  Ab initio real space code with multigrid acceleration
 *  Quantum molecular dynamics package.
 *  Version: 2.1.5
 *  COPYRIGHT
 *  Copyright (C) 1995  Emil Briggs
 *  Copyright (C) 1998  Emil Briggs, Charles Brabec, Mark Wensell, 
 *  Dan Sullivan, Chris Rapcewicz, Jerzy Bernholc
 *  Copyright (C) 2001  Emil Briggs, Wenchang Lu,
 *  Marco Buongiorno Nardelli,Charles Brabec, 
 *  Mark Wensell,Dan Sullivan, Chris Rapcewicz,
 *  Jerzy Bernholc
 */

#if !(defined(_WIN32) || defined(_WIN64))
#include <semaphore.h>
#endif
#include <mpi.h>
#include "typedefs.h"
#include "rmgtypes.h"
#include "rmgtypedefs.h"
#include "rmgthreads.h"


#define HYBRID_EIG 1
#define HYBRID_SKIP 2
#define HYBRID_SUBDIAG_APP_A 3
#define HYBRID_SUBDIAG_APP_B 4
#define HYBRID_BETAX_PSI1_CALCULATE 5
#define HYBRID_SUBDIAG_APP_AB 6


void thread_barrier_wait(void);


int get_thread_basetag(void);
int get_thread_tid(void);
SCF_THREAD_CONTROL *get_thread_control(void);
void run_thread_tasks(int jobs);
void init_HYBRID_MODEL(int npes, int thispe, int nthreads, MPI_Comm comm);
int is_loop_over_states(void);

