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

#include "common_prototypes.h"
#include "common_prototypes1.h"
#include "transition.h"
#include "prototypes_tddft.h"

#if HIP_ENABLED
#include <hip/hip_runtime.h>
#endif

#if CUDA_ENABLED
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cublas_v2.h>
#endif


template void HmatrixUpdate<double>(Kpoint<double> *, double *, double *, int tddft_start_state);
template void HmatrixUpdate<std::complex<double> >(Kpoint<std::complex<double>> *, double *, std::complex<double> *, int tddft_start_state);

template <typename KpointType>
void HmatrixUpdate (Kpoint<KpointType> *kptr, double *vtot_eig, KpointType *Aij, int tddft_start_state)
{

    BaseGrid *G = kptr->G;
    Lattice *L = kptr->L;

    int num_states = kptr->nstates - tddft_start_state;
    int pbasis = kptr->pbasis;
    double vel = L->get_omega() / ((double)(G->get_NX_GRID(1) * G->get_NY_GRID(1) * G->get_NZ_GRID(1)));

    static KpointType *tmp_arrayT;
    static KpointType *global_matrix1;
    size_t psi_alloc = (size_t)ct.max_states * (size_t)pbasis * sizeof(KpointType);

    int factor = 1;
    if(!ct.is_gamma) factor = 2;

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

#if CUDA_ENABLED || HIP_ENABLED
    KpointType *psi_dev = (KpointType *)kptr->psi_dev;
    psi_dev = psi_dev + tddft_start_state * pbasis;
    KpointType *work_dev = (KpointType *)kptr->work_dev;
    static double *v_dev;
    static KpointType *mat_dev;
    gpublasStatus_t gstat;

    int block_size = ct.scalapack_block_factor;
     block_size = num_states;
    int nblock = (num_states + block_size -1)/block_size;


    if(!v_dev)
    {
        gpuMalloc((void **)&v_dev, pbasis * sizeof(double));
        gpuMalloc((void **)&mat_dev, num_states * block_size * sizeof(KpointType));
        int retval1 = MPI_Alloc_mem(num_states * block_size * sizeof(KpointType) , MPI_INFO_NULL, &global_matrix1);

        if(retval1 != MPI_SUCCESS) {
            rmg_error_handler (__FILE__, __LINE__, "Memory allocation failure in HmatrixUpdate");
        }
    }

    gpuMemcpy(v_dev, vtot_eig,  pbasis * sizeof(double), gpuMemcpyHostToDevice);
    gstat = gpublasDdgmm(ct.gpublas_handle, GPUBLAS_SIDE_LEFT, pbasis, num_states, 
            (double *)psi_dev, pbasis, (double *)v_dev, 1, (double *)work_dev, pbasis);
    RmgGpuError(__FILE__, __LINE__, gstat, "Error performing gpublasDgmm.");

    // V|psi> is in work_dev now
    // Compute A matrix
    KpointType alpha(vel);
    KpointType beta(0.0);
    for(int j = 0; j < nblock; j++)
    {
        int size_col = std::min(block_size, num_states - j * block_size);
        RmgGemm(trans_a, trans_n, num_states, size_col,  pbasis, alpha, psi_dev, pbasis, work_dev + j * block_size * pbasis, 
                pbasis, beta, mat_dev, num_states);
        gpuMemcpy(global_matrix1, mat_dev,  (size_t)num_states * (size_t)size_col * sizeof(KpointType), gpuMemcpyDeviceToHost);

        BlockAllreduce((double *)global_matrix1, (size_t)num_states * (size_t)size_col * (size_t)factor , pct.grid_comm);

        for(int jst = 0; jst < size_col; jst++)
        {
            for(int ist = 0; ist < num_states; ist++)
            {
                int idx1 = ist + (j * block_size + jst) * num_states;
                int idx2 = ist * num_states + (j * block_size + jst);
                Aij[idx1] = global_matrix1[jst * num_states + ist];
                Aij[idx2] = MyConj(global_matrix1[jst * num_states + ist]);

            }
        }
    }

#else

    // First time through allocate pinned memory for global_matrix1
    if(!global_matrix1) {

        int retval1 = MPI_Alloc_mem(num_states * num_states * sizeof(KpointType) , MPI_INFO_NULL, &global_matrix1);

        if(retval1 != MPI_SUCCESS) {
            rmg_error_handler (__FILE__, __LINE__, "Memory allocation failure in HmatrixUpdate");
        }

    }

    // V|psi> is in tmp_arrayT
    tmp_arrayT = new KpointType[psi_alloc];
    for (int st1 = 0; st1 < num_states; st1++)
    {
        int st2 = st1 + tddft_start_state;
        for(int idx = 0; idx <pbasis; idx++)
        {
            tmp_arrayT[st1 * pbasis + idx] = kptr->Kstates[st2].psi[idx] * vtot_eig[idx];
        } 
    }

    // Compute A matrix
    KpointType alpha(vel);
    KpointType beta(0.0);
    KpointType *psi = &kptr->orbital_storage[tddft_start_state * pbasis];
    RmgGemm(trans_a, trans_n, num_states, num_states, pbasis, alpha, psi, pbasis, tmp_arrayT, 
            pbasis, beta, global_matrix1, num_states);

    delete [] tmp_arrayT;


    //    MPI_Allreduce(MPI_IN_PLACE, (double *)global_matrix1, num_states * num_states * factor, MPI_DOUBLE, MPI_SUM, pct.grid_comm);
    BlockAllreduce((double *)global_matrix1, (size_t)num_states * (size_t)num_states * (size_t)factor , pct.grid_comm);

    // Store reduced Aij back in Aij matrix
    for(int idx = 0;idx < num_states*num_states;idx++) Aij[idx] = global_matrix1[idx];
#endif

}

