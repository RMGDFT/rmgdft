/*
 *
 * Copyright 2014 The RMG Project Developers. See the COPYRIGHT file 
 * at the top-level directory of this distribution or in the current
 * directory.
 * 
 * This file is part of RMG. 
 * RMG is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * any later version.
 *
 * RMG is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

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
#include "RmgParallelFft.h"

#include "../../RMG/Headers/prototypes.h"
#include "common_prototypes.h"
#include "common_prototypes1.h"
#include "prototypes_tddft.h"
#include "transition.h"

#if GPU_ENABLED
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cublas_v2.h>
#endif


template void HSmatrix<double>(Kpoint<double> *, double *, double *, double *);
template void HSmatrix<std::complex<double> >(Kpoint<std::complex<double>> *, double *, std::complex<double> *, std::complex<double> *); 

    template <typename KpointType>
void HSmatrix (Kpoint<KpointType> *kptr, double *vtot_eig, KpointType *Aij, KpointType *Sij) 
{
    //rmg_printf("\nSUBSPACE DIAGONALIZATION\n");

    BaseGrid *G = kptr->G;
    Lattice *L = kptr->L;

    int num_states = kptr->nstates;
    int pbasis = kptr->pbasis;
    double vel = L->get_omega() / ((double)(G->get_NX_GRID(1) * G->get_NY_GRID(1) * G->get_NZ_GRID(1)));

    // For MPI routines
    int factor = 1;
    if(!ct.is_gamma) factor = 2;

    BaseThread *T = BaseThread::getBaseThread(0);

    static KpointType *tmp_arrayT;
    static KpointType *global_matrix1;
    static KpointType *global_matrix2;

    KpointType *Agpu = NULL;
    KpointType *NULLptr = NULL;

#if GPU_ENABLED

    cublasStatus_t custat;
    KpointType *tmp_array2T = (KpointType *)GpuMallocHost(pbasis * kptr->nstates * sizeof(KpointType));     

#else

    KpointType *tmp_array2T = new KpointType[pbasis * kptr->nstates];

#endif



    char *trans_t = "t";
    char *trans_n = "n";
    char *trans_c = "c";
    char *trans_a;
    if(typeid(KpointType) == typeid(std::complex<double>)) {
        trans_a = trans_c;
    }
    else {
        trans_a = trans_t;
    }   


    // First time through allocate pinned memory for buffers
    if(!tmp_arrayT) {

        int retval1 = MPI_Alloc_mem(pbasis * ct.max_states * sizeof(KpointType) , MPI_INFO_NULL, &tmp_arrayT);
        int retval2 = MPI_Alloc_mem(ct.max_states * ct.max_states * sizeof(KpointType) , MPI_INFO_NULL, &global_matrix1);
        int retval3 = MPI_Alloc_mem(ct.max_states * ct.max_states * sizeof(KpointType) , MPI_INFO_NULL, &global_matrix2);

        if((retval1 != MPI_SUCCESS) || (retval2 != MPI_SUCCESS) || (retval3 != MPI_SUCCESS)) {
            rmg_error_handler (__FILE__, __LINE__, "Memory allocation failure in Subdiag");
        }

#if GPU_ENABLED
        RmgCudaError(__FILE__, __LINE__, cudaHostRegister( tmp_arrayT, pbasis * ct.max_states * sizeof(KpointType), cudaHostRegisterPortable), "Error registering memory.\n");
        RmgCudaError(__FILE__, __LINE__, cudaHostRegister( global_matrix1, ct.max_states * ct.max_states * sizeof(KpointType), cudaHostRegisterPortable), "Error registering memory.\n");
        RmgCudaError(__FILE__, __LINE__, cudaHostRegister( global_matrix2, ct.max_states * ct.max_states * sizeof(KpointType), cudaHostRegisterPortable), "Error registering memory.\n");
#endif

    }


    // Apply operators on each wavefunction
    RmgTimer *RT1 = new RmgTimer("Diagonalization: apply operators");

    // Apply Nls
    AppNls(kptr, kptr->newsint_local, kptr->Kstates[0].psi, kptr->nv, kptr->ns, kptr->Bns,
            0, std::min(ct.non_local_block_size, kptr->nstates));
    int first_nls = 0;


    // Each thread applies the operator to one wavefunction
    KpointType *a_psi = (KpointType *)tmp_arrayT;
    KpointType *b_psi = (KpointType *)tmp_array2T;

    int istop = kptr->nstates / T->get_threads_per_node();
    istop = istop * T->get_threads_per_node();
    for(int st1=0;st1 < istop;st1 += T->get_threads_per_node()) {
        SCF_THREAD_CONTROL thread_control[MAX_RMG_THREADS];
        // Make sure the non-local operators are applied for the next block if needed
        int check = first_nls + T->get_threads_per_node();
        if(check > ct.non_local_block_size) {
            AppNls(kptr, kptr->newsint_local, kptr->Kstates[st1].psi, kptr->nv, &kptr->ns[st1 * pbasis], kptr->Bns,
                    st1, std::min(ct.non_local_block_size, kptr->nstates - st1));
            first_nls = 0;
        }

        for(int ist = 0;ist < ct.THREADS_PER_NODE;ist++) {
            thread_control[ist].job = HYBRID_SUBDIAG_APP_AB;
            thread_control[ist].sp = &kptr->Kstates[st1 + ist];
            thread_control[ist].p1 = (void *)&a_psi[(st1 + ist) * pbasis];
            thread_control[ist].p2 = (void *)&b_psi[(st1 + ist) * pbasis];
            thread_control[ist].p3 = (void *)kptr;
            thread_control[ist].vtot = vtot_eig;
            thread_control[ist].nv = (void *)&kptr->nv[(first_nls + ist) * pbasis];
            thread_control[ist].Bns = (void *)&kptr->Bns[(first_nls + ist) * pbasis];

            T->set_pptr(ist, &thread_control[ist]);
        }

        // Thread tasks are set up so wake them
        T->run_thread_tasks(T->get_threads_per_node());

        // Increment index into non-local block
        first_nls += T->get_threads_per_node();

    }

    // Process any remaining orbitals serially
    for(int st1 = istop;st1 < kptr->nstates;st1++) {
        ApplyOperators (kptr, st1, &a_psi[st1 * pbasis], &b_psi[st1 * pbasis], vtot_eig, 
                &kptr->nv[first_nls * pbasis], &kptr->Bns[first_nls * pbasis]);
        first_nls++;
    }

    delete(RT1);
    /* Operators applied and we now have
tmp_arrayT:  A|psi> + BV|psi> + B|beta>dnm<beta|psi>
tmp_array2T:  B|psi> + B|beta>qnm<beta|psi> */

    // Compute A matrix
    RT1 = new RmgTimer("Diagonalization: matrix setup/reduce");
    KpointType alpha(1.0);
    KpointType beta(0.0);
    RmgGemm(trans_a, trans_n, num_states, num_states, pbasis, alpha, kptr->orbital_storage, pbasis, tmp_arrayT, pbasis, beta, global_matrix1, num_states, Agpu, NULLptr, NULLptr, false, true, false, true);

#if HAVE_ASYNC_ALLREDUCE
    // Asynchronously reduce it
    MPI_Request MPI_reqAij;
    MPI_Iallreduce(MPI_IN_PLACE, (double *)global_matrix1, num_states * num_states * factor, MPI_DOUBLE, MPI_SUM, pct.grid_comm, &MPI_reqAij);
#else
    MPI_Allreduce(MPI_IN_PLACE, (double *)global_matrix1, num_states * num_states * factor, MPI_DOUBLE, MPI_SUM, pct.grid_comm);
#endif

    // Compute S matrix
    KpointType alpha1(vel);
    RmgGemm (trans_a, trans_n, num_states, num_states, pbasis, alpha1, kptr->orbital_storage, pbasis, kptr->ns, pbasis, beta, global_matrix2, num_states, Agpu, NULLptr, NULLptr, false, true, false, true);

#if HAVE_ASYNC_ALLREDUCE
    // Wait for Aij request to finish
    MPI_Wait(&MPI_reqAij, MPI_STATUS_IGNORE);
#endif

#if HAVE_ASYNC_ALLREDUCE
    // Asynchronously reduce Sij request
    MPI_Request MPI_reqSij;
    MPI_Iallreduce(MPI_IN_PLACE, (double *)global_matrix2, num_states * num_states * factor, MPI_DOUBLE, MPI_SUM, pct.grid_comm, &MPI_reqSij);
#else
    MPI_Allreduce(MPI_IN_PLACE, (double *)global_matrix2, num_states * num_states * factor, MPI_DOUBLE, MPI_SUM, pct.grid_comm);
#endif

    // Store reduced Aij back in Aij matrix
    for(int idx = 0;idx < num_states*num_states;idx++) Aij[idx] = global_matrix1[idx];


#if HAVE_ASYNC_ALLREDUCE
    // Wait for S request to finish and when done store copy in Sij
    MPI_Wait(&MPI_reqSij, MPI_STATUS_IGNORE);
#endif
    for(int idx = 0;idx < num_states*num_states;idx++) Sij[idx] = global_matrix2[idx];
    delete(RT1);

    // Free up tmp_array2T
#if GPU_ENABLED
    GpuFreeHost(tmp_array2T);
#else
    delete [] tmp_array2T;
#endif


}

