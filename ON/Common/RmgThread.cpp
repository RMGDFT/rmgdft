#include "BaseThread.h"
#include "RmgThread.h"
#include "rmg_error.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "rmgthreads.h"
#include "const.h"
//#include "prototypes.h"
#if GPU_ENABLED
#include <cuda.h>
#include <cuda_runtime.h>
#endif


// Main thread function specific to subprojects
void *run_threads(void *v) {

    int retval;
    BaseThreadControl *s;
    SCF_THREAD_CONTROL *ss;
    s = (BaseThreadControl *)v;
    BaseThread *T = BaseThread::getBaseThread(0);


#if GPU_ENABLED
    cudaError_t cuerr;
#endif

    T->set_cpu_affinity(s->tid, pct.procs_per_host, pct.local_rank);

    // Set up thread local storage
    rmg_set_tsd(s);

    // Get the control structure
    ss = (SCF_THREAD_CONTROL *)s->pptr;

#if GPU_ENABLED
    cudaSetDevice(ct.cu_dev); 
#endif

    // Quick hack fix
    sleep(1000000);
    T->thread_sleep();
#if 0
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
               subdiag_app_AB_one(ss->sp, (double *)ss->p1, (double *)ss->p2, ss->vtot);
               break;
            case HYBRID_SUBDIAG_APP_A:
               subdiag_app_A_one(ss->sp, (double *)ss->p1, (double *)ss->p2, ss->vtot);
               break;
            case HYBRID_SUBDIAG_APP_B:
               subdiag_app_B_one(ss->sp, (double *)ss->p1);
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

#endif
}


