#include <complex>
#include "FiniteDiff.h"
#include "const.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "rmg_error.h"
#include "rmgthreads.h"
#include "RmgTimer.h"
#include "RmgThread.h"
#include "GlobalSums.h"
#include "Kpoint.h"
#include "RmgGemm.h"
#include "Gpufuncs.h"
#include "Subdiag.h"
#include "GpuAlloc.h"
#include "ErrorFuncs.h"
#include "blas.h"
#include "Solvers.h"

#include "common_prototypes.h"
#include "common_prototypes1.h"
#include "transition.h"
#include "prototypes_tddft.h"

#if GPU_ENABLED
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cublas_v2.h>
#endif


template void HSmatrix<double>(Kpoint<double> *, double *, double *, double *, double *);
template void HSmatrix<std::complex<double> >(Kpoint<std::complex<double>> *, double *,double *,  std::complex<double> *, std::complex<double> *); 

    template <typename KpointType>
void HSmatrix (Kpoint<KpointType> *kptr, double *vtot_eig,double *vxc_psi,  KpointType *Hmat, KpointType *Smat) 
{
    RmgTimer RT0("4-Diagonalization");
    bool potential_acceleration = (ct.potential_acceleration_constant_step > 0.0) && (ct.scf_steps > 0);

    double vel = kptr->L->get_omega() / ((double)(kptr->G->get_NX_GRID(1) * kptr->G->get_NY_GRID(1) * kptr->G->get_NZ_GRID(1)));

    int nstates = kptr->nstates;
    MPI_Comm grid_comm = kptr->grid_comm;
    KpointType *orbital_storage = kptr->orbital_storage;
    KpointType *ns = kptr->ns;
    KpointType *nv = kptr->nv;
    KpointType *newsint_local = kptr->newsint_local;
    int pbasis_noncoll = kptr->pbasis * ct.noncoll_factor;

    // For MPI routines
    int factor = 1;
    if(!ct.is_gamma) factor = 2;

    
    BaseThread *T = BaseThread::getBaseThread(0);
    // State array is 4 * the number of states in length but memory above
    // the first set of nstates is unused in this routine so we can use it
    // as temporary space.
    KpointType *tmp_arrayT = kptr->Kstates[0].psi;
    tmp_arrayT += kptr->nstates * pbasis_noncoll;

    static KpointType *global_matrix1;

#if GPU_ENABLED

    KpointType *Aij = (KpointType *)GpuMallocManaged(nstates * nstates * sizeof(KpointType));
    KpointType *Sij = (KpointType *)GpuMallocManaged(nstates * nstates * sizeof(KpointType));
    if(!global_matrix1) global_matrix1 = (KpointType *)GpuMallocManaged(nstates * nstates * sizeof(KpointType));     
    GpuFill((double *)Aij, factor*nstates * nstates, 0.0);
    GpuFill((double *)Sij, factor*nstates * nstates, 0.0);

#else
    KpointType *Aij = new KpointType[nstates * nstates]();
    KpointType *Sij = new KpointType[nstates * nstates];
    if(!global_matrix1) {
        int retval1 = MPI_Alloc_mem(ct.max_states * ct.max_states * sizeof(KpointType) , MPI_INFO_NULL, &global_matrix1);
        if(retval1 != MPI_SUCCESS) rmg_error_handler (__FILE__, __LINE__, "Memory allocation failure in Subdiag");
    }
    for(int ix=0;ix < nstates*nstates;ix++)global_matrix1[ix] = 0.0;
#endif


    char *trans_t = "t";
    char *trans_n = "n";
    char *trans_c = "c";
    char *trans_a = trans_t;
    if(typeid(KpointType) == typeid(std::complex<double>)) trans_a = trans_c;


    // Apply operators on each wavefunction
    RmgTimer *RT1 = new RmgTimer("4-Diagonalization: apply operators");
    RmgTimer *RT2 = new RmgTimer("4-Diagonalization: AppNls");

    // Apply Nls
    AppNls(kptr, newsint_local, kptr->Kstates[0].psi, nv, ns, 0, std::min(ct.non_local_block_size, nstates));
    delete RT2;
    int first_nls = 0;

    // Each thread applies the operator to one wavefunction
    KpointType *h_psi = (KpointType *)tmp_arrayT;

#if GPU_ENABLED
    // Until the finite difference operators are being applied on the GPU it's faster
    // to make sure that the result arrays are present on the cpu side.
    int device = -1;
    cudaMemPrefetchAsync ( h_psi, nstates*pbasis_noncoll*sizeof(KpointType), device, NULL);
    cudaDeviceSynchronize();
#endif

    int active_threads = ct.MG_THREADS_PER_NODE;
    if(ct.mpi_queue_mode) active_threads--;
    if(active_threads < 1) active_threads = 1;

    int istop = nstates / active_threads;
    istop = istop * active_threads;
    for(int st1=0;st1 < istop;st1 += active_threads) {
        SCF_THREAD_CONTROL thread_control;
        // Make sure the non-local operators are applied for the next block if needed
         int check = first_nls + active_threads;
         if(check > ct.non_local_block_size) {
             RmgTimer *RT3 = new RmgTimer("4-Diagonalization: apply operators: AppNls");
#if GPU_ENABLED
             cudaDeviceSynchronize();
#endif
             AppNls(kptr, newsint_local, kptr->Kstates[st1].psi, nv, &ns[st1 * pbasis_noncoll],
                    st1, std::min(ct.non_local_block_size, nstates - st1));
#if GPU_ENABLED
             cudaDeviceSynchronize();
#endif
             first_nls = 0;
             delete RT3;
         }

        for(int ist = 0;ist < active_threads;ist++) {
            thread_control.job = HYBRID_APPLY_HAMILTONIAN;
            thread_control.vtot = vtot_eig;
            thread_control.vxc_psi = vxc_psi;
            thread_control.extratag1 = potential_acceleration;
            thread_control.istate = st1 + ist;
            thread_control.sp = &kptr->Kstates[st1 + ist];
            thread_control.p1 = (void *)kptr->Kstates[st1 + ist].psi;
            thread_control.p2 = (void *)&h_psi[(st1 + ist) * pbasis_noncoll];
            thread_control.p3 = (void *)kptr;
            thread_control.nv = (void *)&nv[(first_nls + ist) * pbasis_noncoll];
            thread_control.ns = (void *)&ns[(st1 + ist) * pbasis_noncoll];  // ns is not blocked!
            thread_control.basetag = kptr->Kstates[st1 + ist].istate;
            QueueThreadTask(ist, thread_control);
        }

        // Thread tasks are set up so wake them
        if(!ct.mpi_queue_mode) T->run_thread_tasks(active_threads);
        if((check >= ct.non_local_block_size) && ct.mpi_queue_mode) T->run_thread_tasks(active_threads, Rmg_Q);

        // Increment index into non-local block
        first_nls += active_threads;

    }

    if(ct.mpi_queue_mode) T->run_thread_tasks(active_threads, Rmg_Q);

    // Process any remaining orbitals serially
    for(int st1 = istop;st1 < nstates;st1++) {
        // Make sure the non-local operators are applied for the next state if needed
         int check = first_nls + 1;
         if(check > ct.non_local_block_size) {
             RmgTimer *RT3 = new RmgTimer("4-Diagonalization: apply operators: AppNls");
#if GPU_ENABLED
             cudaDeviceSynchronize();
#endif
             AppNls(kptr, newsint_local, kptr->Kstates[st1].psi, nv, &ns[st1 * pbasis_noncoll], st1, std::min(ct.non_local_block_size, nstates - st1));
#if GPU_ENABLED
             cudaDeviceSynchronize();
#endif
             first_nls = 0;
             delete RT3;
         }
        ApplyHamiltonian<KpointType, KpointType> (kptr, st1, kptr->Kstates[st1].psi, &h_psi[st1 * pbasis_noncoll], 
                      vtot_eig, vxc_psi, &nv[first_nls * pbasis_noncoll], potential_acceleration);

        first_nls++;
    }
    delete(RT1);
    /* Operators applied and we now have
         tmp_arrayT:  A|psi> + BV|psi> + B|beta>dnm<beta|psi>
         tmp_array2T:  B|psi> + B|beta>qnm<beta|psi> */

#if GPU_ENABLED
    cudaDeviceSynchronize();
#endif

    // Compute A matrix
    RT1 = new RmgTimer("4-Diagonalization: matrix setup/reduce");
    KpointType alpha(vel);
    KpointType beta(0.0);

    if(ct.is_gamma)
        RmgSyrkx("L", "T", nstates, pbasis_noncoll, alpha, orbital_storage, pbasis_noncoll, tmp_arrayT, pbasis_noncoll, beta, Aij, nstates);
    else
        RmgGemm(trans_a, trans_n, nstates, nstates, pbasis_noncoll, alpha, orbital_storage, pbasis_noncoll, tmp_arrayT, pbasis_noncoll, beta, Aij, nstates);

#if HAVE_ASYNC_ALLREDUCE
    // Asynchronously reduce it
    MPI_Request MPI_reqAij;
    MPI_Request MPI_reqSij;
    if(ct.use_async_allreduce)
       MPI_Iallreduce(MPI_IN_PLACE, (double *)Aij, nstates * nstates * factor, MPI_DOUBLE, MPI_SUM, grid_comm, &MPI_reqAij);
    else
       MPI_Allreduce(MPI_IN_PLACE, (double *)Aij, nstates * nstates * factor, MPI_DOUBLE, MPI_SUM, grid_comm);
#else
    MPI_Allreduce(MPI_IN_PLACE, (double *)Aij, nstates * nstates * factor, MPI_DOUBLE, MPI_SUM, grid_comm);
#endif

    // Compute S matrix
    KpointType alpha1(vel);
    if(ct.norm_conserving_pp)
    {
        RmgSyrkx("L", "T", nstates, pbasis_noncoll, alpha1, orbital_storage, pbasis_noncoll,  orbital_storage, pbasis_noncoll, beta, Sij, nstates);
    }
    else
    {
        RmgGemm (trans_a, trans_n, nstates, nstates, pbasis_noncoll, alpha1, orbital_storage, pbasis_noncoll, ns, pbasis_noncoll, beta, Sij, nstates);
    }


#if HAVE_ASYNC_ALLREDUCE
    // Asynchronously reduce Sij request
    if(ct.use_async_allreduce)
        MPI_Iallreduce(MPI_IN_PLACE, (double *)Sij, nstates * nstates * factor, MPI_DOUBLE, MPI_SUM, grid_comm, &MPI_reqSij);
    else
        MPI_Allreduce(MPI_IN_PLACE, (double *)Sij, nstates * nstates * factor, MPI_DOUBLE, MPI_SUM, grid_comm);
#else
    MPI_Allreduce(MPI_IN_PLACE, (double *)Sij, nstates * nstates * factor, MPI_DOUBLE, MPI_SUM, grid_comm);
#endif


#if HAVE_ASYNC_ALLREDUCE
    // Wait for S request to finish and when done store copy in Sij
    if(ct.use_async_allreduce) MPI_Wait(&MPI_reqAij, MPI_STATUS_IGNORE);
    if(ct.use_async_allreduce) MPI_Wait(&MPI_reqSij, MPI_STATUS_IGNORE);
#endif

    delete(RT1);

    for(int i = 0; i < nstates * nstates; i++)
    {
        Hmat[i] = Aij[i];
        Smat[i] = Sij[i];
    }
    // free memory
#if GPU_ENABLED
    GpuFreeManaged(Sij);
    GpuFreeManaged(Aij);
#else
    delete [] Sij;
    delete [] Aij;
#endif



}

