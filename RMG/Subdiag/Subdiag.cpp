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
#include "RmgParallelFft.h"

#include "common_prototypes.h"
#include "common_prototypes1.h"
#include "transition.h"

#if GPU_ENABLED
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cublas_v2.h>
#endif


template void Kpoint<double>::Subdiag(double *, int);
template void Kpoint<std::complex<double>>::Subdiag(double *, int);

template <class KpointType> void Kpoint<KpointType>::Subdiag (double *vtot_eig, int subdiag_driver)
{
    RmgTimer RT0("4-Diagonalization");

    double vel = L->get_omega() / ((double)(G->get_NX_GRID(1) * G->get_NY_GRID(1) * G->get_NZ_GRID(1)));


    // For MPI routines
    int factor = 1;
    if(!ct.is_gamma) factor = 2;

    BaseThread *T = BaseThread::getBaseThread(0);
    // State array is 4 * the number of states in length but memory above
    // the first set of nstates is unused in this routine so we can use it
    // as temporary space.
    KpointType *tmp_arrayT = Kstates[0].psi;
    tmp_arrayT += nstates * pbasis;

    static KpointType *global_matrix1;

#if GPU_ENABLED

    KpointType *Aij = (KpointType *)GpuMallocManaged(nstates * nstates * sizeof(KpointType));
    KpointType *Bij = (KpointType *)GpuMallocManaged(nstates * nstates * sizeof(KpointType));
    KpointType *Sij = (KpointType *)GpuMallocManaged(nstates * nstates * sizeof(KpointType));
    KpointType *tmp_array2T = (KpointType *)GpuMallocManaged(pbasis * nstates * sizeof(KpointType));     
    if(!global_matrix1) global_matrix1 = (KpointType *)GpuMallocManaged(nstates * nstates * sizeof(KpointType));     
    double *eigs = (double *)GpuMallocManaged(2*nstates * sizeof(double));
    GpuFill((double *)Aij, factor*nstates * nstates, 0.0);
    GpuFill((double *)Sij, factor*nstates * nstates, 0.0);

#else
    KpointType *Aij = new KpointType[nstates * nstates]();
    KpointType *Bij = new KpointType[nstates * nstates];
    KpointType *Sij = new KpointType[nstates * nstates];
    KpointType *tmp_array2T = new KpointType[pbasis * nstates];
    if(!global_matrix1) {
        int retval1 = MPI_Alloc_mem(ct.max_states * ct.max_states * sizeof(KpointType) , MPI_INFO_NULL, &global_matrix1);
        if(retval1 != MPI_SUCCESS) rmg_error_handler (__FILE__, __LINE__, "Memory allocation failure in Subdiag");
    }
    for(int ix=0;ix < nstates*nstates;ix++)global_matrix1[ix] = 0.0;
    double *eigs = new double[2*nstates];
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
    AppNls(this, newsint_local, Kstates[0].psi, nv, ns, Bns, 0, std::min(ct.non_local_block_size, nstates));
    delete RT2;
    int first_nls = 0;

    // Each thread applies the operator to one wavefunction
    KpointType *a_psi = (KpointType *)tmp_arrayT;
    KpointType *b_psi = (KpointType *)tmp_array2T;

#if GPU_ENABLED
    // Until the finite difference operators are being applied on the GPU it's faster
    // to make sure that the result arrays are present on the cpu side.
    int device = -1;
    cudaMemPrefetchAsync ( a_psi, nstates*pbasis*sizeof(KpointType), device, NULL);
    cudaMemPrefetchAsync ( b_psi, nstates*pbasis*sizeof(KpointType), device, NULL);
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
             AppNls(this, newsint_local, Kstates[st1].psi, nv, &ns[st1 * pbasis], Bns,
                    st1, std::min(ct.non_local_block_size, nstates - st1));
#if GPU_ENABLED
             cudaDeviceSynchronize();
#endif
             first_nls = 0;
             delete RT3;
         }

        for(int ist = 0;ist < active_threads;ist++) {
            thread_control.job = HYBRID_SUBDIAG_APP_AB;
            thread_control.sp = &Kstates[st1 + ist];
            thread_control.p1 = (void *)&a_psi[(st1 + ist) * pbasis];
            thread_control.p2 = (void *)&b_psi[(st1 + ist) * pbasis];
            thread_control.p3 = (void *)this;
            thread_control.vtot = vtot_eig;
            thread_control.nv = (void *)&nv[(first_nls + ist) * pbasis];
            thread_control.Bns = (void *)&Bns[(first_nls + ist) * pbasis];
            thread_control.basetag = Kstates[st1 + ist].istate;
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
             AppNls(this, newsint_local, Kstates[st1].psi, nv, &ns[st1 * pbasis], Bns, st1, std::min(ct.non_local_block_size, nstates - st1));
#if GPU_ENABLED
             cudaDeviceSynchronize();
#endif
             first_nls = 0;
             delete RT3;
         }
        ApplyOperators (this, st1, &a_psi[st1 * pbasis], &b_psi[st1 * pbasis], vtot_eig, 
                       &nv[first_nls * pbasis], &Bns[first_nls * pbasis]);
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
    KpointType alpha(1.0);
    KpointType beta(0.0);

    if(ct.is_gamma)
        RmgSyrkx("L", "T", nstates, pbasis, alpha, orbital_storage, pbasis, tmp_arrayT, pbasis, beta, Aij, nstates);
    else
        RmgGemm(trans_a, trans_n, nstates, nstates, pbasis, alpha, orbital_storage, pbasis, tmp_arrayT, pbasis, beta, Aij, nstates);

#if HAVE_ASYNC_ALLREDUCE
    // Asynchronously reduce it
    MPI_Request MPI_reqAij;
    MPI_Request MPI_reqSij;
    MPI_Request MPI_reqBij;
    if(ct.use_async_allreduce)
       MPI_Iallreduce(MPI_IN_PLACE, (double *)Aij, nstates * nstates * factor, MPI_DOUBLE, MPI_SUM, grid_comm, &MPI_reqAij);
    else
       MPI_Allreduce(MPI_IN_PLACE, (double *)Aij, nstates * nstates * factor, MPI_DOUBLE, MPI_SUM, grid_comm);
#else
    MPI_Allreduce(MPI_IN_PLACE, (double *)Aij, nstates * nstates * factor, MPI_DOUBLE, MPI_SUM, grid_comm);
#endif

    // Compute S matrix
    KpointType alpha1(vel);
    if(ct.norm_conserving_pp && (ct.discretization != MEHRSTELLEN_DISCRETIZATION) && ct.is_gamma)
    {
        RmgSyrkx("L", "T", nstates, pbasis, alpha1, orbital_storage, pbasis,  orbital_storage, pbasis, beta, Sij, nstates);
    }
    else
    {
        RmgGemm (trans_a, trans_n, nstates, nstates, pbasis, alpha1, orbital_storage, pbasis, ns, pbasis, beta, Sij, nstates);
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


    // We need B matrix for US pseudopotentials and/or MEHRSTELLEN_DISCRETIZATION
    if(!ct.norm_conserving_pp || (ct.norm_conserving_pp && ct.discretization == MEHRSTELLEN_DISCRETIZATION)) {

        // Compute B matrix
        RmgGemm (trans_a, trans_n, nstates, nstates, pbasis, alpha1, orbital_storage, pbasis, tmp_array2T, pbasis, beta, Bij, nstates);

        // Reduce matrix and store copy in Bij
#if HAVE_ASYNC_ALLREDUCE
        if(ct.use_async_allreduce)
            MPI_Iallreduce(MPI_IN_PLACE, (double *)Bij, nstates * nstates * factor, MPI_DOUBLE, MPI_SUM, grid_comm, &MPI_reqBij);
        else
            MPI_Allreduce(MPI_IN_PLACE, (double *)Bij, nstates * nstates * factor, MPI_DOUBLE, MPI_SUM, grid_comm);
#else
        MPI_Allreduce(MPI_IN_PLACE, (double *)Bij, nstates * nstates * factor, MPI_DOUBLE, MPI_SUM, grid_comm);
#endif

    }

#if HAVE_ASYNC_ALLREDUCE
    // Wait for S request to finish and when done store copy in Sij
    if(ct.use_async_allreduce) MPI_Wait(&MPI_reqAij, MPI_STATUS_IGNORE);
    if(ct.use_async_allreduce) MPI_Wait(&MPI_reqSij, MPI_STATUS_IGNORE);
    if(!ct.norm_conserving_pp || (ct.norm_conserving_pp && ct.discretization == MEHRSTELLEN_DISCRETIZATION))
        if(ct.use_async_allreduce) MPI_Wait(&MPI_reqBij, MPI_STATUS_IGNORE);
#endif

    delete(RT1);

    // Free up tmp_array2T
#if GPU_ENABLED
    GpuFreeManaged(tmp_array2T);
#else
    delete [] tmp_array2T;
#endif

    // Dispatch to correct subroutine, eigs will hold eigenvalues on return and global_matrix1 will hold the eigenvectors.
    // The eigenvectors may be stored in row-major or column-major format depending on the type of diagonaliztion method
    // used. This is handled during the rotation of the orbitals by trans_b which is set by the driver routine.
    RT1 = new RmgTimer("4-Diagonalization: Eigensolver");
    char *trans_b = "n";
    switch(subdiag_driver) {

        case SUBDIAG_LAPACK:
            trans_b = Subdiag_Lapack (this, Aij, Bij, Sij, eigs, global_matrix1);
            break;
        case SUBDIAG_SCALAPACK:
            trans_b = Subdiag_Scalapack (this, Aij, Bij, Sij, eigs, global_matrix1);
            break;
        case SUBDIAG_ELPA:
            trans_b = Subdiag_Elpa (this, Aij, Bij, Sij, eigs, global_matrix1);
            break;
        case SUBDIAG_MAGMA:
        case SUBDIAG_CUSOLVER:
#if GPU_ENABLED
            trans_b = Subdiag_Cusolver (this, Aij, Bij, Sij, eigs, global_matrix1);
#else
            trans_b = Subdiag_Lapack (this, Aij, Bij, Sij, eigs, global_matrix1);
#endif
            break;
        default:
            rmg_error_handler(__FILE__, __LINE__, "Invalid subdiag_driver type");

    } // end switch
    delete(RT1);

    // If subspace diagonalization is used every step, use eigenvalues obtained here 
    // as the correct eigenvalues
    if (ct.diag == 1) {
        for (int st1 = 0; st1 < nstates; st1++) {
            Kstates[st1].eig[0] = eigs[st1];
        }
    }

    // free memory
#if GPU_ENABLED
    GpuFreeManaged(Sij);
    GpuFreeManaged(Bij);
    GpuFreeManaged(Aij);
#else
    delete [] Sij;
    delete [] Bij;
    delete [] Aij;
#endif


    // Update the orbitals
    RT1 = new RmgTimer("4-Diagonalization: Update orbitals");

    RmgGemm(trans_n, trans_b, pbasis, nstates, nstates, alpha, 
            orbital_storage, pbasis, global_matrix1, nstates, beta, tmp_arrayT, pbasis);

    // And finally copy them back
    int istart = 0;
    int tlen = nstates * pbasis * sizeof(KpointType); 
    if(Verify ("freeze_occupied", true, ControlMap))
    {
        for(int istate = 0;istate < nstates;istate++)
        {
            if(Kstates[istate].occupation[0] > 1.0e-10) highest_occupied = istate;
        }
        istart = (highest_occupied + 1)*pbasis;
        tlen = nstates * pbasis - (highest_occupied + 1) * pbasis;
    }

    // And finally make sure they follow the same sign convention when using hybrid XC
    // Optimize this for GPUs!
    if(ct.xc_is_hybrid)
    {
        for(int istate=0;istate < nstates;istate++)
        {
            if(std::real(global_matrix1[istate*nstates + istate]) < 0.0)
            {
                for(int idx=0;idx < pbasis;idx++) Kstates[istate].psi[idx] = -Kstates[istate].psi[idx];
            }
        }
    }

#if GPU_ENABLED
    cudaMemcpy(&orbital_storage[istart], &tmp_arrayT[istart], tlen, cudaMemcpyDefault);
    //cudaMemPrefetchAsync (orbital_storage , nstates*sizeof(double), cudaCpuDeviceId, NULL);
    // Not sure why but the cudaMemcpy behaves strangely here sometimes.
    //for(int idx=0;idx<nstates*pbasis;idx++) orbital_storage[istart+idx] = tmp_arrayT[istart+idx];
#else
    memcpy(&orbital_storage[istart], &tmp_arrayT[istart], tlen);
#endif

    delete(RT1);

#if GPU_ENABLED
    // After the first step this matrix does not need to be as large
    if(ct.scf_steps == 0) {GpuFreeManaged(global_matrix1);global_matrix1 = NULL;}
    GpuFreeManaged(eigs);
#else
    delete [] eigs;
#endif

}

