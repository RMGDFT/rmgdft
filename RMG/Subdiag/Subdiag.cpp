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
#include "RmgGemm.h"
#include "Subdiag.h"
#include "GpuAlloc.h"
#include "ErrorFuncs.h"
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




template void Subdiag<double>(Kpoint<double> *, double *, double *, double *, int);
template void Subdiag<std::complex<double> >(Kpoint<std::complex<double>> *, double *, double *, double *, int);

template <typename KpointType>
void Subdiag (Kpoint<KpointType> *kptr, double *vh, double *vnuc, double *vxc, int subdiag_driver)
{
    RmgTimer RT0("Diagonalization");
    rmg_printf("\nSUBSPACE DIAGONALIZATION\n");

    BaseGrid *G = kptr->G;
    Lattice *L = kptr->L;

    int num_states = kptr->nstates;
    int pbasis = kptr->pbasis;
    double vel = L->get_omega() / ((double)(G->get_NX_GRID(1) * G->get_NY_GRID(1) * G->get_NZ_GRID(1)));

    // For MPI routines
    int factor = 1;
    if(!ct.is_gamma) factor = 2;

    BaseThread *T = BaseThread::getBaseThread(0);

    static KpointType *tmp_arrayT = NULL;
    static KpointType *tmp_array2T = NULL;
    static KpointType *global_matrix = NULL;

    KpointType *Agpu = NULL;
    KpointType *gpu_eigvectors = NULL;
    KpointType *NULLptr = NULL;

#if GPU_ENABLED

    cublasStatus_t custat;
    static KpointType *Aij = NULL;
    static KpointType *Bij = NULL;
    static KpointType *Sij = NULL;

    // Grab some gpu memory
    gpu_eigvectors = (KpointType *)GpuMalloc(num_states * num_states * sizeof( KpointType ));
    Agpu = (KpointType *)GpuMalloc(pbasis * num_states * sizeof( KpointType ));

    // Start wavefunctions transferring to the GPU
    RmgCudaError(__FILE__,__LINE__, cublasSetVector(pbasis * num_states, sizeof( KpointType ), kptr->orbital_storage, 1, Agpu, 1 ) , "Problem transferring orbitals to GPU");


#else

    KpointType *Aij = new KpointType[kptr->nstates * kptr->nstates];
    KpointType *Bij = new KpointType[kptr->nstates * kptr->nstates];
    KpointType *Sij = new KpointType[kptr->nstates * kptr->nstates];

#endif



    char *trans_t = "t";
    char *trans_n = "n";
    char *trans_c = "c";
    char *trans_m;
    if(typeid(KpointType) == typeid(std::complex<double>)) {
         trans_m = trans_c;
    }
    else {
        trans_m = trans_t;
    }   


    // First time through allocate pinned memory for buffers
    if(!tmp_arrayT) {

        int retval1 = MPI_Alloc_mem(kptr->pbasis * kptr->nstates * sizeof(KpointType) , MPI_INFO_NULL, &tmp_arrayT);
        int retval2 = MPI_Alloc_mem(kptr->pbasis * kptr->nstates * sizeof(KpointType) , MPI_INFO_NULL, &tmp_array2T);
        int retval3 = MPI_Alloc_mem(kptr->nstates * kptr->nstates * sizeof(KpointType) , MPI_INFO_NULL, &global_matrix);

        if((retval1 != MPI_SUCCESS) || (retval2 != MPI_SUCCESS) || (retval3 != MPI_SUCCESS) ) {
            rmg_error_handler (__FILE__, __LINE__, "Memory allocation failure in Subdiag");
        }

        #if GPU_ENABLED
            RmgCudaError(__FILE__, __LINE__, cudaHostRegister( tmp_arrayT, kptr->pbasis * kptr->nstates * sizeof(KpointType), cudaHostRegisterPortable), "Error registering memory.\n");
            RmgCudaError(__FILE__, __LINE__, cudaHostRegister( tmp_array2T, kptr->pbasis * kptr->nstates * sizeof(KpointType), cudaHostRegisterPortable), "Error registering memory.\n");
            RmgCudaError(__FILE__, __LINE__, cudaHostRegister( global_matrix, kptr->nstates * kptr->nstates * sizeof(KpointType), cudaHostRegisterPortable), "Error registering memory.\n");
            RmgCudaError(__FILE__, __LINE__, cudaMallocHost((void **)&Aij, kptr->nstates * kptr->nstates * sizeof(KpointType)), "Error allocating memory.\n");
            RmgCudaError(__FILE__, __LINE__, cudaMallocHost((void **)&Bij, kptr->nstates * kptr->nstates * sizeof(KpointType)), "Error allocating memory.\n");
            RmgCudaError(__FILE__, __LINE__, cudaMallocHost((void **)&Sij, kptr->nstates * kptr->nstates * sizeof(KpointType)), "Error allocating memory.\n");
        #endif

    }


    // Get vtot on fine grid 
    int FP0_BASIS = kptr->G->get_P0_BASIS(kptr->G->get_default_FG_RATIO());
    double *vtot = new double[4*FP0_BASIS];
    double *vtot_eig = new double[kptr->pbasis];
    double *eigs = new double[2*kptr->nstates];

    for (int idx = 0; idx < FP0_BASIS; idx++) {
        vtot[idx] = vh[idx] + vxc[idx] + vnuc[idx];
    }

    // vtot holds potential on fine grid so restrict it to the orbital grid and store
    // result in vtot_eig
    get_vtot_psi (vtot_eig, vtot, kptr->G->get_default_FG_RATIO());


    // Apply operators on each wavefunction
    RmgTimer *RT1 = new RmgTimer("Diagonalization: apply operators");


    // Apply Nls
    AppNls(kptr, kptr->newsint_local);


    // Each thread applies the operator to one wavefunction
    KpointType *a_psi = (KpointType *)tmp_arrayT;
    KpointType *b_psi = (KpointType *)tmp_array2T;

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
         tmp_arrayT:  A|psi> + BV|psi> + B|beta>dnm<beta|psi>
         tmp_array2T:  B|psi> + B|beta>qnm<beta|psi> */


    // Compute A matrix
    RT1 = new RmgTimer("Diagonalization: matrix setup");
    KpointType alpha(1.0);
    KpointType beta(0.0);
    RmgGemm(trans_m, trans_n, num_states, num_states, pbasis, alpha, kptr->orbital_storage, pbasis, tmp_arrayT, pbasis, beta, global_matrix, num_states, Agpu, NULLptr, NULLptr, false, true, false, true);
    delete(RT1);

    // Reduce matrix and store copy in Aij
    RT1 = new RmgTimer("Diagonalization: MPI_Allreduce");
    MPI_Allreduce(MPI_IN_PLACE, (double *)global_matrix, num_states * num_states * factor, MPI_DOUBLE, MPI_SUM, pct.grid_comm);

    for(int idx = 0;idx < num_states*num_states;idx++) Aij[idx] = global_matrix[idx];
    delete(RT1);



    // Compute S matrix
    RT1 = new RmgTimer("Diagonalization: matrix setup");
    KpointType alpha1(vel);
    RmgGemm (trans_m, trans_n, num_states, num_states, pbasis, alpha1, kptr->orbital_storage, pbasis, kptr->ns, pbasis, beta, global_matrix, num_states, Agpu, NULLptr, NULLptr, false, true, false, true);
    delete(RT1);

    // Reduce matrix and store copy in Sij
    RT1 = new RmgTimer("Diagonalization: MPI_Allreduce");
    MPI_Allreduce(MPI_IN_PLACE, (double *)global_matrix, num_states * num_states * factor, MPI_DOUBLE, MPI_SUM, pct.grid_comm);
    for(int idx = 0;idx < num_states*num_states;idx++) Sij[idx] = global_matrix[idx];
    delete(RT1);


    // We need B matrix for US pseudopotentials and/or MEHRSTELLEN_DISCRETIZATION
    if(!ct.norm_conserving_pp || (ct.norm_conserving_pp && ct.discretization == MEHRSTELLEN_DISCRETIZATION)) {

        // Compute B matrix
        RT1 = new RmgTimer("Diagonalization: matrix setup");
        RmgGemm (trans_m, trans_n, num_states, num_states, pbasis, alpha1, kptr->orbital_storage, pbasis, tmp_array2T, pbasis, beta, global_matrix, num_states, Agpu, NULLptr, NULLptr, false, true, false, true);
        delete(RT1);

        // Reduce matrix and store copy in Bij
        RT1 = new RmgTimer("Diagonalization: MPI_Allreduce");
        MPI_Allreduce(MPI_IN_PLACE, (double *)global_matrix, num_states * num_states * factor, MPI_DOUBLE, MPI_SUM, pct.grid_comm);
        for(int idx = 0;idx < num_states*num_states;idx++) Bij[idx] = global_matrix[idx];
        delete(RT1);

    }


    // global_matrix holds Bij now, we store a copy in Bij as well and pass Bij to the driver routine in globalmatrix as well

    // Dispatch to correct subroutine, eigs will hold eigenvalues on return and global_matrix will hold the eigenvectors.
    // If the MAGMA driver is selected the eigenvectors will be stored on the GPU in gpu_global_matrix
    RT1 = new RmgTimer("Diagonalization: Eigensolver");
    switch(subdiag_driver) {

        case SUBDIAG_LAPACK:
            Subdiag_Lapack (kptr, (KpointType *)Aij, (KpointType *)Bij, (KpointType *)Sij, eigs, (KpointType *)global_matrix);
            break;
        case SUBDIAG_SCALAPACK:
            Subdiag_Scalapack (kptr, (KpointType *)Aij, (KpointType *)Bij, (KpointType *)Sij, eigs, (KpointType *)global_matrix);
            break;
        case SUBDIAG_MAGMA:
            Subdiag_Magma (kptr, (KpointType *)Aij, (KpointType *)Bij, (KpointType *)Sij, eigs, (KpointType *)global_matrix, (KpointType *)gpu_eigvectors);
            break;
        default:
            rmg_error_handler(__FILE__, __LINE__, "Invalid subdiag_driver type");

    } // end switch
    delete(RT1);


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
    if(subdiag_driver == SUBDIAG_MAGMA) {

        // The magma driver leaves the eigenvectors on the gpu in gpu_eigvectors and the current
        // orbitals are already there in Agpu. So we set C to kptr->orbital_storage and have the
        // gpu transfer the rotated orbitals directly back there.
        // The magma driver also brings the eigvectors back from the gpu and puts them in global_matrix
        // so we can rotate the betaxpsi as well.
        RmgGemm(trans_n, trans_n, pbasis, num_states, num_states, alpha, 
                NULLptr, pbasis, NULLptr, num_states, beta, kptr->orbital_storage, pbasis, 
                Agpu, gpu_eigvectors, NULLptr, false, false, false, true);
    }
    else {

        RmgGemm(trans_n, trans_n, pbasis, num_states, num_states, alpha, 
                kptr->orbital_storage, pbasis, global_matrix, num_states, beta, tmp_arrayT, pbasis, 
                Agpu, NULLptr, NULLptr, false, true, false, true);

        // And finally copy them back
        for(int idx = 0;idx < num_states * pbasis;idx++) kptr->orbital_storage[idx] = tmp_arrayT[idx];
    }

    delete(RT1);

    // Rotate new betaxpsi
    int size = kptr->sint_size;
    KpointType *sint_ptr = kptr->newsint_local;
    KpointType *work = new KpointType[num_states];
    KpointType *tmp_sint_ptr = new KpointType[size];
    KpointType ZERO_t(0.0);

    for(int ion=0;ion < pct.num_nonloc_ions;ion++) {

        for(int ip = 0; ip < ct.max_nl; ip++) {
            for(int st1=0;st1 < num_states;st1++) {
                work[st1] = ZERO_t;
                for(int st2=0;st2 < num_states;st2++) {
                    work[st1] += global_matrix[st1*num_states + st2] * sint_ptr[ion*num_states*ct.max_nl  + st2*ct.max_nl + ip];
                }
            }
            for(int st1=0;st1 < num_states;st1++) {
                sint_ptr[ion*num_states*ct.max_nl + st1*ct.max_nl + ip] = work[st1];
            }
        }

    }
    delete [] tmp_sint_ptr;
    delete [] work;

#if GPU_ENABLED
    GpuFree(Agpu);
    GpuFree(gpu_eigvectors);
#endif

    // free memory
    delete [] eigs;
#if !GPU_ENABLED
    delete [] Sij;
    delete [] Bij;
    delete [] Aij;
#endif
    delete [] vtot_eig;
    delete [] vtot;

}

