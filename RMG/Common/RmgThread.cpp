#include "BaseThread.h"
#include "RmgThread.h"
#include "rmg_error.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "rmgthreads.h"
#include "const.h"
#include "prototypes.h"
using namespace std;


// Main thread function specific to subprojects
void *run_threads(void *v) {

    int retval;
    BaseThread *s;
    SCF_THREAD_CONTROL *ss;
    s = (BaseThread *)v;
    

#if GPU_ENABLED
    cudaError_t cuerr;
#endif

    s->set_cpu_affinity(s->tid);

    s->pthread_tid = pthread_self();
    pthread_setspecific(s->scf_thread_control_key, (void *)s);

#if GPU_ENABLED
    cudaSetDevice(ct.cu_dev); 
    if(cudaSuccess != (cuerr = cudaStreamCreate(&s->cstream))) {
        fprintf (stderr, "Error: cudaStreamCreate failed for: threads setup\n");
        exit(-1);
    }
    if( cudaSuccess != cudaMallocHost((void **)&s->gpu_host_temp1, (get_PX0_GRID() + 4) * (get_PY0_GRID() + 4) * (get_PZ0_GRID() + 4) * sizeof(rmg_double_t) )){
        fprintf (stderr, "Error: cudaMallocHost failed for: threads gpu_host_temp\n");
        exit(-1);
    }
    if( cudaSuccess != cudaMallocHost((void **)&s->gpu_host_temp2, (get_PX0_GRID() + 4) * (get_PY0_GRID() + 4) * (get_PZ0_GRID() + 4) * sizeof(rmg_double_t) )){
        fprintf (stderr, "Error: cudaMallocHost failed for: threads gpu_host_temp\n");
        exit(-1);
    }
    cuMemHostRegister ( s->trade_buf, sizeof(rmg_double_t) * (alloc + 64), CU_MEMHOSTREGISTER_PORTABLE);
#endif

    // Wait until everyone gets here
    pthread_barrier_wait(&s->run_barrier);

    while(1) {

        // We sleep forever or until we get a signal that wakes us up
        sem_wait(&s->this_sync);

        // Get the control structure
        ss = (SCF_THREAD_CONTROL *)s->pptr;

        // Switch that controls what we do
        switch(ss->job) {
            case HYBRID_EIG:       // Performs a single multigrid sweep over an orbital
               mg_eig_state_driver(ss->sp, 0, ss->vtot, ct.mg_eig_precision);
               break;
            case HYBRID_SKIP:
               break;
#if GAMMA_PT
            case HYBRID_SUBDIAG_APP_AB:
               subdiag_app_AB_one(ss->sp, (rmg_double_t *)ss->p1, (rmg_double_t *)ss->p2, ss->vtot);
               break;
            case HYBRID_SUBDIAG_APP_A:
               subdiag_app_A_one(ss->sp, (rmg_double_t *)ss->p1, (rmg_double_t *)ss->p2, ss->vtot);
               break;
            case HYBRID_SUBDIAG_APP_B:
               subdiag_app_B_one(ss->sp, (rmg_double_t *)ss->p1);
               break;
#endif 
            case HYBRID_BETAX_PSI1_CALCULATE:
               betaxpsi1_calculate_one(ss->sp, ss->ion, ss->nion, ss->sintR, ss->sintI, ss->kpt, ss->weiptr);
               break;
            default:
               break;
        }

        // Let the main thread know that we are done
        sem_post(&s->thread_sem);

    }

}


