/*
 *
 * Copyright 2025 The RMG Project Developers. See the COPYRIGHT file 
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



template void Kpoint<double>::KineticEnergyDensity(double *ked);
template void Kpoint<std::complex<double>>::KineticEnergyDensity(double *ked);

template <class KpointType> void Kpoint<KpointType>::KineticEnergyDensity (double *ked)
{
    RmgTimer RT0("4-KEd");

    double vel = L->get_omega() / (double)G->get_GLOBAL_BASIS(1);
    wfobj<KpointType> gx, gy, gz;
    int pbasis_noncoll = pbasis * ct.noncoll_factor;

    BaseThread *T = BaseThread::getBaseThread(0);
    // State array is 4 * the number of states in length but memory above
    // the first set of nstates is unused in this routine so we can use it
    // as temporary space.
    KpointType *tmp_arrayT = Kstates[0].psi;
    tmp_arrayT += nstates * pbasis_noncoll ;

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


#if CUDA_ENABLED
    DeviceSynchronize();
#endif

    int density = 1;
    for(int st = 0; st < nstates; st++) {
        ApplyGradient (Kstates[st].psi, gx.data(), gy.data(), gz.data(), ct.kohn_sham_fd_order, "Coarse");

        //if(ct.noncoll){
        //    ApplyAOperator<KpointType>(&Kstates[st].psi[pbasis], &h_psi[st * pbasis_noncoll+pbasis], 
        //            dimx, dimy, dimz, gridhx, gridhy, gridhz, ct.kohn_sham_fd_order, this->kp.kvec);
        //}

        for(int idx = 0;idx < pbasis_noncoll;idx++){
            double t1 = std::real(gx[idx]*std::conj(gx[idx]) + gy[idx]*std::conj(gy[idx]) + gz[idx]*std::conj(gz[idx]));
            h_psi[st * pbasis_noncoll + idx] = 0.5*t1;
        }

        double weight = this->kp.kweight * Kstates[st].occupation[0];
        for(int idx=0;idx < pbasis_noncoll;idx++)
        {
            ked[idx] += weight * std::real(h_psi[st * pbasis_noncoll + idx]);
        }
    }

#if CUDA_ENABLED
    DeviceSynchronize();
#endif

}

