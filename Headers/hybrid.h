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

#include <semaphore.h>
#include <mpi.h>
#include "typedefs.h"
#include "rmgtypes.h"
#include "rmgtypedefs.h"
#include "rmgthreads.h"

#if HYBRID_MODEL

#define HYBRID_EIG 1
#define HYBRID_SKIP 2
#define HYBRID_SUBDIAG_APP_A 3
#define HYBRID_SUBDIAG_APP_B 4
#define HYBRID_BETAX_PSI1_CALCULATE 5
#define HYBRID_SUBDIAG_APP_AB 6


void scf_barrier_init(int nthreads);
void scf_barrier_wait(void);
void scf_barrier_destroy(void);


void scf_tsd_init(void);
void scf_tsd_set_value(void *s);
void scf_tsd_delete(void);
int get_thread_basetag(void);
int get_thread_tid(void);
SCF_THREAD_CONTROL *get_thread_control(void);
void set_cpu_affinity(int tid);
void wait_for_threads(int jobs);
void wake_threads(int jobs);
void init_HYBRID_MODEL(void);
void enter_threaded_region(void);
void leave_threaded_region(void);
int is_loop_over_states(void);

#if GPU_ENABLED
cudaStream_t *get_thread_cstream(void);
#endif
#endif
