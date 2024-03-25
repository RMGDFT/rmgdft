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
#include "Solvers.h"
#include "Functional.h"
#include "RmgMatrix.h"

#include "common_prototypes.h"
#include "common_prototypes1.h"
#include "transition.h"

#if CUDA_ENABLED
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cublas_v2.h>
#endif



template void Kpoint<double>::ComputeHcore(double *, double *, double *Hij, double *Hij_kin, double *Hij_localpp);
template void Kpoint<std::complex<double>>::ComputeHcore(double *, double *, std::complex<double> *Hij,
    std::complex<double> *Hij_kin, std::complex<double> *Hij_localpp);

template <class KpointType> void Kpoint<KpointType>::ComputeHcore (double *vtot_eig, double *vxc_psi, 
        KpointType *Hij, KpointType *Hij_kin, KpointType *Hij_localpp)
{
    RmgTimer RT0("4-Hcore");

    double vel = L->get_omega() / ((double)(G->get_NX_GRID(1) * G->get_NY_GRID(1) * G->get_NZ_GRID(1)));

    bool potential_acceleration = false;


    bool is_xc_hybrid = ct.xc_is_hybrid;
    ct.xc_is_hybrid = false;


    // For MPI routines
    int factor = 1;
    if(!ct.is_gamma) factor = 2;

    int pbasis_noncoll = pbasis * ct.noncoll_factor;
    BaseThread *T = BaseThread::getBaseThread(0);
    // State array is 4 * the number of states in length but memory above
    // the first set of nstates is unused in this routine so we can use it
    // as temporary space.
    KpointType *tmp_arrayT = Kstates[0].psi;
    tmp_arrayT += nstates * pbasis_noncoll ;

    char *trans_t = "t";
    char *trans_n = "n";
    char *trans_c = "c";
    char *trans_a = trans_t;
    if(typeid(KpointType) == typeid(std::complex<double>)) trans_a = trans_c;

    // Apply operators on each wavefunction
    RmgTimer *RT1 = new RmgTimer("4-Hcore: apply operators");
    RmgTimer *RT2 = new RmgTimer("4-Hcore: AppNls");

    // Apply Nls
    AppNls(this, newsint_local, Kstates[0].psi, nv, ns, 0, std::min(ct.non_local_block_size, nstates));

    delete RT2;
    int first_nls = 0;

    // Each thread applies the operator to one wavefunction
    KpointType *h_psi = (KpointType *)tmp_arrayT;

#if CUDA_ENABLED || HIP_ENABLED
    // Until the finite difference operators are being applied on the GPU it's faster
    // to make sure that the result arrays are present on the cpu side.
    int device = -1;
    gpuMemPrefetchAsync ( h_psi, nstates*pbasis_noncoll*sizeof(KpointType), device, NULL);
    DeviceSynchronize();
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
            DeviceSynchronize();
            AppNls(this, newsint_local, Kstates[st1].psi, nv, &ns[st1 * pbasis_noncoll],
                    st1, std::min(ct.non_local_block_size, nstates - st1));
            DeviceSynchronize();
            first_nls = 0;
            delete RT3;
        }

        for(int ist = 0;ist < active_threads;ist++) {
            thread_control.job = HYBRID_APPLY_HAMILTONIAN;
            thread_control.vtot = vtot_eig;
            thread_control.vxc_psi = vxc_psi;
            thread_control.extratag1 = potential_acceleration;
            thread_control.istate = st1 + ist;
            thread_control.sp = &Kstates[st1 + ist];
            thread_control.p1 = (void *)Kstates[st1 + ist].psi;
            thread_control.p2 = (void *)&h_psi[(st1 + ist) * pbasis_noncoll];
            thread_control.p3 = (void *)this;
            thread_control.nv = (void *)&nv[(first_nls + ist) * pbasis_noncoll];
            thread_control.ns = (void *)&ns[(st1 + ist) * pbasis_noncoll];  // ns is not blocked!
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
#if CUDA_ENABLED
            DeviceSynchronize();
#endif
            AppNls(this, newsint_local, Kstates[st1].psi, nv, &ns[st1 * pbasis_noncoll], st1, std::min(ct.non_local_block_size, nstates - st1));
#if CUDA_ENABLED
            DeviceSynchronize();
#endif
            first_nls = 0;
            delete RT3;
        }
        ApplyHamiltonian<KpointType, KpointType> (this, st1, Kstates[st1].psi, &h_psi[st1 * pbasis_noncoll], vtot_eig, vxc_psi, &nv[first_nls * pbasis_noncoll], potential_acceleration);

        first_nls++;
    }
    delete(RT1);
    /* Operators applied and we now have
tmp_arrayT:  A|psi> + BV|psi> + B|beta>dnm<beta|psi> */

#if CUDA_ENABLED
    DeviceSynchronize();
#endif

    // Compute A matrix
    RT1 = new RmgTimer("4-Hcore: matrix setup/reduce");
    KpointType alphavel(vel);
    KpointType beta(0.0);

    RmgGemm(trans_a, trans_n, nstates, nstates, pbasis_noncoll, alphavel, orbital_storage, pbasis_noncoll, tmp_arrayT, pbasis_noncoll, beta, Hij, nstates);

    BlockAllreduce((double *)Hij, (size_t)(nstates)*(size_t)nstates * (size_t)factor, grid_comm);
    if(Hij_kin == NULL && Hij_localpp == NULL)
    {
        ct.xc_is_hybrid = is_xc_hybrid;
        return;
    }


    int density = 1;
    int dimx = this->G->get_PX0_GRID(density);
    int dimy = this->G->get_PY0_GRID(density);
    int dimz = this->G->get_PZ0_GRID(density);
    for(int st = 0; st < nstates; st++) {
        double gridhx = this->G->get_hxgrid(density);
        double gridhy = this->G->get_hygrid(density);
        double gridhz = this->G->get_hzgrid(density);
        ApplyAOperator<KpointType>(Kstates[st].psi, &h_psi[st * pbasis_noncoll], 
                dimx, dimy, dimz, gridhx, gridhy, gridhz, ct.kohn_sham_fd_order, this->kp.kvec);
        if(ct.noncoll){
            ApplyAOperator<KpointType>(&Kstates[st].psi[pbasis], &h_psi[st * pbasis_noncoll+pbasis], 
                    dimx, dimy, dimz, gridhx, gridhy, gridhz, ct.kohn_sham_fd_order, this->kp.kvec);
        }

        double tmag(0.5 * this->kp.kmag);
        for(int idx = 0;idx < pbasis_noncoll;idx++){
            h_psi[st * pbasis_noncoll + idx] = -0.5 * h_psi[st * pbasis_noncoll + idx] + tmag*Kstates[st].psi[idx];
        }

    }

    // - 0.5 for laplacian 
    alphavel = vel;
    RmgGemm(trans_a, trans_n, nstates, nstates, pbasis_noncoll, alphavel, orbital_storage, 
            pbasis_noncoll, h_psi, pbasis_noncoll, beta, Hij_kin, nstates);

    BlockAllreduce((double *)Hij_kin, (size_t)(nstates)*(size_t)nstates * (size_t)factor, grid_comm);

    delete(RT1);

#if CUDA_ENABLED
    DeviceSynchronize();
#endif
    for(int st = 0; st < nstates; st++) {
        for (int idx = 0; idx < pbasis; idx++)
        {
            h_psi[st * pbasis_noncoll + idx] = vtot_eig[idx] * Kstates[st].psi[idx];
        }
        if(ct.noncoll)
        {
            for (int idx = 0; idx < pbasis; idx++)
            {
                h_psi[st * pbasis_noncoll + pbasis + idx] = vtot_eig[idx] * Kstates[st].psi[pbasis + idx];
            }
        }
    }

    alphavel = vel;
    RmgGemm(trans_a, trans_n, nstates, nstates, pbasis_noncoll, alphavel, orbital_storage, 
            pbasis_noncoll, h_psi, pbasis_noncoll, beta, Hij_localpp, nstates);
    BlockAllreduce((double *)Hij_localpp, (size_t)(nstates)*(size_t)nstates * (size_t)factor, grid_comm);

    if(pct.gridpe ==0 && ct.verbose)
    {
        rmg_printf("\n Hcore");
        for(int i= 0; i < std::min(8, nstates); i++)
        {
            rmg_printf("\n");
            for(int j = 0; j < std::min(8, nstates); j++)
                rmg_printf(" %10.6f ", std::real(Hij[i*nstates + j]) );
        }
        rmg_printf("\n Hkin");
        for(int i= 0; i < std::min(8, nstates); i++)
        {
            rmg_printf("\n");
            for(int j = 0; j < std::min(8, nstates); j++)
                rmg_printf(" %10.6f ", std::real(Hij_kin[i*nstates + j]) );
        }
        rmg_printf("\n H_localpp");
        for(int i= 0; i < std::min(8, nstates); i++)
        {
            rmg_printf("\n");
            for(int j = 0; j < std::min(8, nstates); j++)
                rmg_printf(" %10.6f ", std::real(Hij_localpp[i*nstates + j]) );
        }

    }


    ct.xc_is_hybrid = is_xc_hybrid;

}

