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

#if HYBRID_MODEL

#define HYBRID_EIG 1
#define HYBRID_SKIP 2
#define HYBRID_SUBDIAG_APP_A 3
#define HYBRID_SUBDIAG_APP_B 4

typedef struct {
  STATE *sp;
  REAL *vtot;
  int basetag;
  int tid;
  REAL *p1;    // p1 and p2 are only used in some parallel regions
  REAL *p2;
} MG_THREAD_STRUCT;

void scf_barrier_init(int nthreads);
void scf_barrier_wait(void);
void scf_barrier_destroy(void);


void scf_tsd_init(void);
void scf_tsd_set_value(void *s);
void scf_tsd_delete(void);
int get_thread_basetag(void);
int get_thread_tid(void);
void mg_eig_state_threaded(MG_THREAD_STRUCT *ss);
void set_cpu_affinity(void);
void wait_for_threads(int jobs);
void wake_threads(int jobs);
void init_HYBRID_MODEL(void);
void enter_threaded_region(void);
void leave_threaded_region(void);
int is_loop_over_states(void);
void RMG_MPI_lock(void);
void RMG_MPI_unlock(void);
void RMG_MPI_thread_order_lock(void);
void RMG_MPI_thread_order_unlock(void);
#endif
