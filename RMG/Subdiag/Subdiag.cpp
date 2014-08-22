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
#include "Kpoint.h"
#include "Subdiag.h"
#include "blas.h"

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
static double *global_matrix_storage = NULL;


template void Subdiag<double>(Kpoint<double> *, double *, double *, double *, int);
template void Subdiag<std::complex<double> >(Kpoint<std::complex<double>> *, double *, double *, double *, int);

template <typename KpointType>
void Subdiag (Kpoint<KpointType> *kptr, double *vh, double *vnuc, double *vxc, int subdiag_driver)
{

    //RmgTimer RT0("Diagonalization");
    BaseGrid *G = kptr->G;
    Lattice *L = kptr->L;
    rmg_printf("\nSUBSPACE DIAGONALIZATION\n");
    KpointType *global_matrix;

    const char *trans_t = "t";
    const char *trans_n = "n";
    const char *trans_c = "c";

    int num_states = kptr->nstates;
    int pbasis = kptr->pbasis;

    double vel = L->get_omega() / ((double)(G->get_NX_GRID(1) * G->get_NY_GRID(1) * G->get_NZ_GRID(1)));

    // For MPI routines
    int factor = 1;
    if(!ct.is_gamma) factor = 2;


    BaseThread *T = BaseThread::getBaseThread(0);


    // First time through allocate pinned memory for buffers
    if(!tmp_arrayR) {

        int retval1 = MPI_Alloc_mem(kptr->pbasis * kptr->nstates * sizeof(KpointType) , MPI_INFO_NULL, &tmp_arrayR);
        int retval2 = MPI_Alloc_mem(kptr->pbasis * kptr->nstates * sizeof(KpointType) , MPI_INFO_NULL, &tmp_array2R);
        int retval3 = MPI_Alloc_mem(kptr->nstates * kptr->nstates * sizeof(KpointType) , MPI_INFO_NULL, &global_matrix_storage);

        if((retval1 != MPI_SUCCESS) || (retval2 != MPI_SUCCESS) || (retval3 != MPI_SUCCESS) ) {
            rmg_error_handler (__FILE__, __LINE__, "Memory allocation failure in Subdiag");
        }

        #if GPU_ENABLED
            cudaHostRegister( tmp_arrayR, kptr->pbasis * kptr->nstates * sizeof(KpointType), cudaHostRegisterPortable);
            cudaHostRegister( tmp_array2R, kptr->pbasis * kptr->nstates * sizeof(KpointType), cudaHostRegisterPortable);
            cudaHostRegister( global_matrix, kptr->nstates * kptr->nstates * sizeof(KpointType), cudaHostRegisterPortable);
        #endif

    }
    global_matrix = (KpointType *)global_matrix_storage;

    // Get vtot on fine grid 
    int FP0_BASIS = kptr->G->get_P0_BASIS(kptr->G->get_default_FG_RATIO());
    double *vtot = new double[4*FP0_BASIS];
    double *vtot_eig = new double[kptr->pbasis];
    KpointType *Aij = new KpointType[kptr->nstates * kptr->nstates];
    KpointType *Bij = new KpointType[kptr->nstates * kptr->nstates];
    KpointType *Sij = new KpointType[kptr->nstates * kptr->nstates];
    double *eigs = new double[kptr->nstates];

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


    // Apply Nls
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
            thread_control[ist].p1 = (void *)&a_psi[(st1 + ist) * kptr->pbasis];
            thread_control[ist].p2 = (void *)&b_psi[(st1 + ist) * kptr->pbasis];
            thread_control[ist].p3 = (void *)kptr;
            thread_control[ist].vtot = vtot_eig;
            T->set_pptr(ist, &thread_control[ist]);
        }

        // Thread tasks are set up so wake them
        T->run_thread_tasks(T->get_threads_per_node());

    }

    // Process any remaining orbitals serially
    for(int st1 = istop;st1 < kptr->nstates;st1++) {
        ApplyOperators (kptr, st1, &a_psi[st1 * kptr->pbasis], &b_psi[st1 * kptr->pbasis], vtot_eig);
    }

    delete(RT1);
    /* Operators applied and we now have
         tmp_arrayR:  A|psi> + BV|psi> + B|beta>dnm<beta|psi>
         tmp_array2R:  B|psi> + B|beta>qnm<beta|psi> */


    // Compute A matrix
    RT1 = new RmgTimer("Diagonalization: matrix setup");

    if(ct.is_gamma) {

        double alpha(1.0);
        double beta(0.0);
        dgemm (trans_t, trans_n, &num_states, &num_states, &pbasis, &alpha, (double *)kptr->Kstates[0].psi, &pbasis,
                                tmp_arrayR, &pbasis, &beta, (double *)global_matrix, &num_states);

    } 
    else {

        KpointType alpha(1.0);
        KpointType beta(0.0);
        zgemm (trans_c, trans_n, &num_states, &num_states, &pbasis, (double *)&alpha, (double *)kptr->Kstates[0].psi, &pbasis,
                                (double *)tmp_arrayR, &pbasis, (double *)&beta, (double *)global_matrix, &num_states);

    }

    delete(RT1);

    // Reduce matrix and store copy in Aij
    RT1 = new RmgTimer("Diagonalization: MPI_Allreduce");
    MPI_Allreduce(MPI_IN_PLACE, (double *)global_matrix, num_states * num_states * factor, MPI_DOUBLE, MPI_SUM, pct.grid_comm);

    for(int idx = 0;idx < num_states*num_states;idx++) Aij[idx] = global_matrix[idx];
    delete(RT1);



    // Compute S matrix
    RT1 = new RmgTimer("Diagonalization: matrix setup");
    if(ct.is_gamma) {

        double alpha(vel);
        double beta(0.0);
        dgemm (trans_t, trans_n, &num_states, &num_states, &pbasis, &alpha, (double *)kptr->Kstates[0].psi, &pbasis,
                                (double *)kptr->ns, &pbasis, &beta, (double *)global_matrix, &num_states);

    }
    else {

        KpointType alpha(vel);
        KpointType beta(0.0);
        zgemm (trans_c, trans_n, &num_states, &num_states, &pbasis, (double *)&alpha, (double *)kptr->Kstates[0].psi, &pbasis,
                                (double *)kptr->ns, &pbasis, (double *)&beta, (double *)global_matrix, &num_states);

    }
    delete(RT1);

    // Reduce matrix and store copy in Sij
    RT1 = new RmgTimer("Diagonalization: MPI_Allreduce");
    MPI_Allreduce(MPI_IN_PLACE, (double *)global_matrix, num_states * num_states * factor, MPI_DOUBLE, MPI_SUM, pct.grid_comm);
    for(int idx = 0;idx < num_states*num_states;idx++) Sij[idx] = global_matrix[idx];
    delete(RT1);



    // Compute B matrix
    RT1 = new RmgTimer("Diagonalization: matrix setup");
    if(ct.is_gamma) {

        double alpha(vel);
        double beta(0.0);
        dgemm (trans_t, trans_n, &num_states, &num_states, &pbasis, &alpha, (double *)kptr->Kstates[0].psi, &pbasis,
                                tmp_array2R, &pbasis, &beta, (double *)global_matrix, &num_states);

    } 
    else {

        KpointType alpha(vel);
        KpointType beta(0.0);
        zgemm (trans_c, trans_n, &num_states, &num_states, &pbasis, (double *)&alpha, (double *)kptr->Kstates[0].psi, &pbasis,
                                (double *)tmp_array2R, &pbasis, (double *)&beta, (double *)global_matrix, &num_states);

    }
    delete(RT1);

    // Reduce matrix and store copy in Bij
    RT1 = new RmgTimer("Diagonalization: MPI_Allreduce");
    MPI_Allreduce(MPI_IN_PLACE, (double *)global_matrix, num_states * num_states * factor, MPI_DOUBLE, MPI_SUM, pct.grid_comm);
    for(int idx = 0;idx < num_states*num_states;idx++) Bij[idx] = global_matrix[idx];
    delete(RT1);


    // global_matrix holds Bij now, we store a copy in Bij as well and pass Bij to the driver routine in globalmatrix as well

    // Dispatch to correct subroutine, eigs will hold eigenvalues on return and global_matrix will hold the eigenvectors
    switch(subdiag_driver) {

        case SUBDIAG_LAPACK:
            Subdiag_Lapack (kptr, (KpointType *)Aij, (KpointType *)Bij, (KpointType *)Sij, eigs, (KpointType *)global_matrix);
            break;
        default:
            rmg_error_handler(__FILE__, __LINE__, "Invalid subdiag_driver type");

    } // end switch


    // If subspace diagonalization is used everystep, use eigenvalues obtained here 
    // as the correct eigenvalues
    if (ct.diag == 1) {
        for (int st1 = 0; st1 < num_states; st1++) {
            kptr->kstates[st1].eig[0] = eigs[st1];
            kptr->Kstates[st1].eig[0] = eigs[st1];
        }
    }


    // Update the orbitals
    RT1 = new RmgTimer("Diagonalization: Update orbitals");
    if(ct.is_gamma) {

        double alpha(1.0);
        double beta(0.0);
        dgemm("n", "n", &pbasis, &num_states, &num_states,
                            (double *)&alpha,
                            (double *)kptr->orbital_storage, &pbasis,
                            (double *)global_matrix, &num_states,
                            (double *)&beta, (double *)tmp_arrayR, &pbasis );
    }
    else {

        KpointType alpha(1.0);
        KpointType beta(0.0);
        zgemm("n", "n", &pbasis, &num_states, &num_states,
                            (double *)&alpha,
                            (double *)kptr->orbital_storage, &pbasis,
                            (double *)global_matrix, &num_states,
                            (double *)&beta, (double *)tmp_arrayR, &pbasis );

    }
    
    KpointType *tptr = (KpointType *)tmp_arrayR;
    for(int idx = 0;idx < num_states * pbasis;idx++) kptr->orbital_storage[idx] = tptr[idx];
    delete(RT1);

    // free memory
    delete [] Sij;
    delete [] Bij;
    delete [] Aij;
    delete [] vtot_eig;
    delete [] vtot;

}

