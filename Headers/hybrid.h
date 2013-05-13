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
#define HYBRID_BETAX_PSI1_CALCULATE 5
#define HYBRID_FINALIZE_PAPI 6
#define HYBRID_SUBDIAG_APP_AB 7



typedef struct {
  STATE *sp;
  rmg_double_t *vtot;
  int basetag;
  int tid;
  rmg_double_t *p1;    // p1 and p2 are only used in some parallel regions
  rmg_double_t *p2;
} MG_THREAD_STRUCT;

void scf_barrier_init(int nthreads);
void scf_barrier_wait(void);
void scf_barrier_destroy(void);


void scf_tsd_init(void);
void scf_tsd_set_value(void *s);
void scf_tsd_delete(void);
int get_thread_basetag(void);
int get_thread_tid(void);
void *get_thread_trade_buf(void);
SCF_THREAD_CONTROL *get_thread_control(void);
void mg_eig_state_threaded(MG_THREAD_STRUCT *ss);
void set_cpu_affinity(int tid);
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
void trade_images_manager(void *s);
long long Papi_thread_flops(int tid);
MPI_Comm *get_thread_grid_comm(void);

#if GPU_ENABLED
cudaStream_t *get_thread_cstream(void);
#endif
#endif
