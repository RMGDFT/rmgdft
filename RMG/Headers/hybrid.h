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
void global_sums_threaded (REAL *vect, int *length, int tid);
void mg_eig_state_threaded(MG_THREAD_STRUCT *ss);
void set_cpu_affinity(void);
