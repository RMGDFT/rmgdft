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

#include "prototypes.h"
#include "common_prototypes.h"
#include "common_prototypes1.h"
#include "transition.h"

#if GPU_ENABLED
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cublas_v2.h>
#endif


template void Subdiag<double>(Kpoint<double> *, double *, int);
template void Subdiag<std::complex<double> >(Kpoint<std::complex<double>> *, double *, int);

// Temporary hack to enable the use of these arrays in GeneralDiag.
void *GMatrix1;
void *GMatrix2;

template <typename KpointType>
void Subdiag (Kpoint<KpointType> *kptr, double *vtot_eig, int subdiag_driver)
{
    RmgTimer RT0("4-Diagonalization");
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

    // State array is 4 * the number of states in length but memory above
    // the first set of nstates is unused in this routine so we can use it
    // as temporary space.
    KpointType *tmp_arrayT = kptr->Kstates[0].psi;
    tmp_arrayT += num_states * pbasis;

    static KpointType *global_matrix1;
    static KpointType *global_matrix2;

    KpointType *Agpu = NULL;
    KpointType *NULLptr = NULL;

#if GPU_ENABLED

    cublasStatus_t custat;
    static KpointType *Aij;
    static KpointType *Bij;
    static KpointType *Sij;
    KpointType *tmp_array2T = (KpointType *)GpuMallocHost(pbasis * kptr->nstates * sizeof(KpointType));     

#else

    KpointType *Aij = new KpointType[kptr->nstates * kptr->nstates];
    KpointType *Bij = new KpointType[kptr->nstates * kptr->nstates];
    KpointType *Sij = new KpointType[kptr->nstates * kptr->nstates];
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
    if(!global_matrix1) {

        int retval1 = MPI_Alloc_mem(ct.max_states * ct.max_states * sizeof(KpointType) , MPI_INFO_NULL, &global_matrix1);
        int retval2 = MPI_Alloc_mem(ct.max_states * ct.max_states * sizeof(KpointType) , MPI_INFO_NULL, &global_matrix2);

        if((retval1 != MPI_SUCCESS) || (retval2 != MPI_SUCCESS)) {
            rmg_error_handler (__FILE__, __LINE__, "Memory allocation failure in Subdiag");
        }

        GMatrix1 = (void *)global_matrix1;
        GMatrix2 = (void *)global_matrix2;

        #if GPU_ENABLED
            RmgCudaError(__FILE__, __LINE__, cudaHostRegister( global_matrix1, ct.max_states * ct.max_states * sizeof(KpointType), cudaHostRegisterPortable), "Error registering memory.\n");
            RmgCudaError(__FILE__, __LINE__, cudaHostRegister( global_matrix2, ct.max_states * ct.max_states * sizeof(KpointType), cudaHostRegisterPortable), "Error registering memory.\n");
            RmgCudaError(__FILE__, __LINE__, cudaMallocHost((void **)&Aij, ct.max_states * ct.max_states * sizeof(KpointType)), "Error allocating memory.\n");
            RmgCudaError(__FILE__, __LINE__, cudaMallocHost((void **)&Bij, ct.max_states * ct.max_states * sizeof(KpointType)), "Error allocating memory.\n");
            RmgCudaError(__FILE__, __LINE__, cudaMallocHost((void **)&Sij, ct.max_states * ct.max_states * sizeof(KpointType)), "Error allocating memory.\n");
        #endif

    }

    
    double *eigs = new double[2*kptr->nstates];


    // Apply operators on each wavefunction
    RmgTimer *RT1 = new RmgTimer("4-Diagonalization: apply operators");
    RmgTimer *RT2 = new RmgTimer("4-Diagonalization: AppNls");

    // Apply Nls
    AppNls(kptr, kptr->newsint_local, kptr->Kstates[0].psi, kptr->nv, kptr->ns, kptr->Bns,
           0, std::min(ct.non_local_block_size, kptr->nstates));
    delete RT2;
    int first_nls = 0;

    // Each thread applies the operator to one wavefunction
    KpointType *a_psi = (KpointType *)tmp_arrayT;
    KpointType *b_psi = (KpointType *)tmp_array2T;

    int active_threads = ct.MG_THREADS_PER_NODE;
    if(ct.mpi_queue_mode) active_threads--;
    if(active_threads < 1) active_threads = 1;

    int istop = kptr->nstates / active_threads;
    istop = istop * active_threads;
    for(int st1=0;st1 < istop;st1 += active_threads) {
        SCF_THREAD_CONTROL thread_control;
        // Make sure the non-local operators are applied for the next block if needed
         int check = first_nls + active_threads;
         if(check > ct.non_local_block_size) {
             RmgTimer *RT3 = new RmgTimer("4-Diagonalization: apply operators: AppNls");

             AppNls(kptr, kptr->newsint_local, kptr->Kstates[st1].psi, kptr->nv, &kptr->ns[st1 * pbasis], kptr->Bns,
                    st1, std::min(ct.non_local_block_size, kptr->nstates - st1));
             first_nls = 0;
             delete RT3;
         }

        for(int ist = 0;ist < active_threads;ist++) {
            thread_control.job = HYBRID_SUBDIAG_APP_AB;
            thread_control.sp = &kptr->Kstates[st1 + ist];
            thread_control.p1 = (void *)&a_psi[(st1 + ist) * pbasis];
            thread_control.p2 = (void *)&b_psi[(st1 + ist) * pbasis];
            thread_control.p3 = (void *)kptr;
            thread_control.vtot = vtot_eig;
            thread_control.nv = (void *)&kptr->nv[(first_nls + ist) * pbasis];
            thread_control.Bns = (void *)&kptr->Bns[(first_nls + ist) * pbasis];
            thread_control.basetag = kptr->Kstates[st1 + ist].istate;
            QueueThreadTask(ist, thread_control);
        }

        // Thread tasks are set up so wake them
        if(!ct.mpi_queue_mode) T->run_thread_tasks(active_threads);

        // Increment index into non-local block
        first_nls += active_threads;

    }

    if(ct.mpi_queue_mode) T->run_thread_tasks(active_threads, Rmg_Q);

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
    RT1 = new RmgTimer("4-Diagonalization: matrix setup/reduce");
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

    // We need B matrix for US pseudopotentials and/or MEHRSTELLEN_DISCRETIZATION
    if(!ct.norm_conserving_pp || (ct.norm_conserving_pp && ct.discretization == MEHRSTELLEN_DISCRETIZATION)) {

        // Compute B matrix
        RmgGemm (trans_a, trans_n, num_states, num_states, pbasis, alpha1, kptr->orbital_storage, pbasis, tmp_array2T, pbasis, beta, global_matrix1, num_states, Agpu, NULLptr, NULLptr, false, true, false, true);

        // Reduce matrix and store copy in Bij
        MPI_Allreduce(MPI_IN_PLACE, (double *)global_matrix1, num_states * num_states * factor, MPI_DOUBLE, MPI_SUM, pct.grid_comm);
        for(int idx = 0;idx < num_states*num_states;idx++) Bij[idx] = global_matrix1[idx];

    }

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

    // global_matrix1 holds Bij now, we store a copy in Bij as well and pass Bij to the driver routine in globalmatrix as well

    // Dispatch to correct subroutine, eigs will hold eigenvalues on return and global_matrix1 will hold the eigenvectors.
    // The eigenvectors may be stored in row-major or column-major format depending on the type of diagonaliztion method
    // used. This is handled during the rotation of the orbitals by trans_b which is set by the driver routine.
    RT1 = new RmgTimer("4-Diagonalization: Eigensolver");
    char *trans_b = "n";
    switch(subdiag_driver) {

        case SUBDIAG_LAPACK:
            trans_b = Subdiag_Lapack (kptr, (KpointType *)Aij, (KpointType *)Bij, (KpointType *)Sij, eigs, (KpointType *)global_matrix1);
            break;
        case SUBDIAG_SCALAPACK:
            trans_b = Subdiag_Scalapack (kptr, (KpointType *)Aij, (KpointType *)Bij, (KpointType *)Sij, eigs, (KpointType *)global_matrix1);
            break;
        case SUBDIAG_ELPA:
            trans_b = Subdiag_Elpa (kptr, (KpointType *)Aij, (KpointType *)Bij, (KpointType *)Sij, eigs, (KpointType *)global_matrix1);
            break;
        case SUBDIAG_MAGMA:
#if GPU_ENABLED && MAGMA_LIBS
            trans_b = Subdiag_Magma (kptr, (KpointType *)Aij, (KpointType *)Bij, (KpointType *)Sij, eigs, (KpointType *)global_matrix1);
#else
            trans_b = Subdiag_Lapack (kptr, (KpointType *)Aij, (KpointType *)Bij, (KpointType *)Sij, eigs, (KpointType *)global_matrix1);
#endif
            break;
        default:
            rmg_error_handler(__FILE__, __LINE__, "Invalid subdiag_driver type");

    } // end switch
    delete(RT1);


    // If subspace diagonalization is used everystep, use eigenvalues obtained here 
    // as the correct eigenvalues
    if (ct.diag == 1) {
        for (int st1 = 0; st1 < num_states; st1++) {
            kptr->Kstates[st1].eig[0] = eigs[st1];
        }
    }

    // Update the orbitals
    RT1 = new RmgTimer("4-Diagonalization: Update orbitals");

    RmgGemm(trans_n, trans_b, pbasis, num_states, num_states, alpha, 
            kptr->orbital_storage, pbasis, global_matrix1, num_states, beta, tmp_arrayT, pbasis, 
            Agpu, NULLptr, NULLptr, false, true, false, true);

    // And finally copy them back
    int istart = 0;
    if(Verify ("freeze_occupied", true, kptr->ControlMap)) istart = kptr->highest_occupied + 1;

    for(int idx = istart;idx < num_states * pbasis;idx++) kptr->orbital_storage[idx] = tmp_arrayT[idx];

    delete(RT1);

    // free memory
    delete [] eigs;
#if !GPU_ENABLED
    delete [] Sij;
    delete [] Bij;
    delete [] Aij;
#endif

}

