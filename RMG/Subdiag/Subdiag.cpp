#include <complex>
#include "FiniteDiff.h"
#include "const.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "rmg_alloc.h"
#include "rmg_error.h"
#include "rmgthreads.h"
#include "RmgTimer.h"
#include "RmgThread.h"
#include "GlobalSums.h"
#include "Subdiag.h"

#include "prototypes.h"
#include "common_prototypes.h"
#include "common_prototypes1.h"
#include "transition.h"

#if GPU_ENABLED
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cublas_v2.h>
#endif


static double *tmp_arrayR = NULL;
static double *tmp_array2R = NULL;


template void Subdiag<double>(Kpoint<double> *, double *, double *, double *, int);
template void Subdiag<std::complex<double> >(Kpoint<std::complex<double>> *, double *, double *, double *, int);

template <typename KpointType>
void Subdiag (Kpoint<KpointType> *kptr, double *vh, double *vnuc, double *vxc, int subdiag_driver)
{

    RmgTimer RT0("Diagonalization");

    BaseThread *T = BaseThread::getBaseThread(0);


    // First time through allocate pinned memory for buffers
    if(!tmp_arrayR) {

        int retval1 = MPI_Alloc_mem(kptr->pbasis * kptr->nstates * sizeof(KpointType) , MPI_INFO_NULL, &tmp_arrayR);
        int retval2 = MPI_Alloc_mem(kptr->pbasis * kptr->nstates * sizeof(KpointType) , MPI_INFO_NULL, &tmp_array2R);
        if((retval1 != MPI_SUCCESS) || (retval2 != MPI_SUCCESS)) {
            rmg_error_handler (__FILE__, __LINE__, "Memory allocation failure in Subdiag");
        }

        #if GPU_ENABLED
            cudaHostRegister( tmp_arrayR, kptr->pbasis * kptr->nstates * sizeof(KpointType), cudaHostRegisterPortable);
            cudaHostRegister( tmp_array2R, kptr->pbasis * kptr->nstates * sizeof(KpointType), cudaHostRegisterPortable);
        #endif

    }


    // Get vtot on fine grid 
    int FP0_BASIS = kptr->G->get_P0_BASIS(kptr->G->get_default_FG_RATIO());
    double *vtot = new double[FP0_BASIS];
    double *vtot_eig = new double[kptr->pbasis];
    for (int idx = 0; idx < FP0_BASIS; idx++) {
        vtot[idx] = vh[idx] + vxc[idx] + vnuc[idx];
    }

    // vtot holds potential on fine grid so restrict it to the orbital grid and store
    // result in vtot_eig
    get_vtot_psi (vtot_eig, vtot, kptr->G->get_default_FG_RATIO());


    // Apply AB operator on each wavefunction
    RmgTimer *RT1 = new RmgTimer("Diagonalization: apply operators");

    int dimx = kptr->G->get_PX0_GRID(1);
    int dimy = kptr->G->get_PY0_GRID(1);
    int dimz = kptr->G->get_PZ0_GRID(1);

    // Trade images on coarse grid and store result back in vtot
    if(ct.kohn_sham_fd_order == APP_CI_FOURTH)
        kptr->T->trade_imagesx (vtot_eig, vtot, dimx, dimy, dimz, 1, FULL_TRADE);

    if(ct.kohn_sham_fd_order == APP_CI_SIXTH)
        kptr->T->trade_imagesx (vtot_eig, vtot, dimx, dimy, dimz, 2, FULL_TRADE);
    

    AppNls(kptr, pct.newsintR_local, pct.newsintI_local);
    

    // Each thread applies the operator to one wavefunction
    KpointType *a_psi = (KpointType *)tmp_arrayR;
    KpointType *b_psi = (KpointType *)tmp_array2R;

    int istop = kptr->nstates / T->get_threads_per_node();
    istop = istop * T->get_threads_per_node();

    for(int st1=0;st1 < istop;st1 += T->get_threads_per_node()) {
        SCF_THREAD_CONTROL thread_control[MAX_RMG_THREADS];
        for(int ist = 0;ist < ct.THREADS_PER_NODE;ist++) {
            thread_control[ist].job = HYBRID_SUBDIAG_APP_AB;
            thread_control[ist].sp = &kptr->kstates[st1 + ist];
            thread_control[ist].p1 = &a_psi[(st1 + ist) * kptr->pbasis];
            thread_control[ist].p2 = &b_psi[(st1 + ist) * kptr->pbasis];
            thread_control[ist].vtot = vtot;
            T->set_pptr(ist, &thread_control[ist]);
        }

        // Thread tasks are set up so wake them
        T->run_thread_tasks(T->get_threads_per_node());

    }

    // Process any remaining orbitals serially
    for(int st1 = istop;st1 < kptr->nstates;st1++) {
        subdiag_app_AB_one (&kptr->kstates[st1], &a_psi[st1 * kptr->pbasis], &b_psi[st1 *  kptr->pbasis], vtot);
    }


    delete(RT1);

    
    // Dispatch to correct subroutine    
    switch(subdiag_driver) {

        case SUBDIAG_LAPACK:
            break;
        default:
            rmg_error_handler(__FILE__, __LINE__, "Invalid subdiag_driver type");

    } // end switch


    // free memory
    delete [] vtot_eig;
    delete [] vtot;
}

