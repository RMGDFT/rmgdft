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

#if CUDA_ENABLED || HIP_ENABLED
void Veff_x_psi(double *psi_dev,  double *work_dev, double *vtot_eig, int pbasis, int num_states);
void Veff_x_psi(std::complex<double> *psi_dev,  std::complex<double> *work_dev, double *vtot_eig, int pbasis, int num_states);
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
    KpointType alpha(vel);
    KpointType beta(0.0);

    static KpointType *global_matrix1;

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
    static KpointType *mat_dev;
    gpublasStatus_t gstat;

    int block_size = std::max(1024, ct.scalapack_block_factor);
    //block_size = num_states;
    int nblock = (num_states + block_size -1)/block_size;

    Veff_x_psi(psi_dev, work_dev, vtot_eig, pbasis, num_states);


    if(!mat_dev)
    {
        gpuMalloc((void **)&mat_dev, num_states * block_size * sizeof(KpointType));
        int retval1 = MPI_Alloc_mem(num_states * block_size * sizeof(KpointType) , MPI_INFO_NULL, &global_matrix1);

        if(retval1 != MPI_SUCCESS) {
            rmg_error_handler (__FILE__, __LINE__, "Memory allocation failure in HmatrixUpdate");
        }
    }

    // V|psi> is in work_dev now
    // Compute A matrix
    for(int j = 0; j < nblock; j++)
    {
        int size_col = std::min(block_size, num_states - j * block_size);
        int size_row = num_states - j * block_size;
        RmgGemm(trans_a, trans_n, size_row, size_col,  pbasis, alpha, psi_dev+ j*block_size*pbasis, pbasis, work_dev + j * block_size * pbasis, 
                pbasis, beta, mat_dev, size_row);
        gpuMemcpy(global_matrix1, mat_dev,  (size_t)size_row * (size_t)size_col * sizeof(KpointType), gpuMemcpyDeviceToHost);

        BlockAllreduce((double *)global_matrix1, (size_t)size_row * (size_t)size_col * (size_t)factor , pct.grid_comm);

        for(int jst = 0; jst < size_col; jst++)
        {
            for(int ist = 0; ist < size_row; ist++)
            {
                int idx1 = (ist + j * block_size) + (j * block_size + jst) * num_states;
                int idx2 = (ist + j * block_size) * num_states + (j * block_size + jst);
                Aij[idx1] = global_matrix1[jst * size_row + ist];
                Aij[idx2] = MyConj(global_matrix1[jst * size_row + ist]);

            }
        }
    }

#else

    int block_size = ct.state_block_size;
    //block_size = num_states;
    int nblock = (num_states + block_size -1)/block_size;
    // First time through allocate pinned memory for global_matrix1
    if(!global_matrix1) {

        int retval1 = MPI_Alloc_mem(num_states * block_size * sizeof(KpointType) , MPI_INFO_NULL, &global_matrix1);

        if(retval1 != MPI_SUCCESS) {
            rmg_error_handler (__FILE__, __LINE__, "Memory allocation failure in HmatrixUpdate");
        }

    }

    // V|psi> is in tmp_arrayT
    KpointType *psi = kptr->orbital_storage + tddft_start_state * pbasis;
    KpointType *vpsi = &kptr->orbital_storage[kptr->nstates * pbasis];  // use the memory of psi extra 3* state_block_size.

    for(int j = 0; j < nblock; j++)
    {
        int size_col = std::min(block_size, num_states - j * block_size);
        int size_row = num_states - j * block_size;

        for (int st1 = 0; st1 < size_col; st1++)
        {
            for(int idx = 0; idx <pbasis; idx++)
            {
                vpsi[st1 * pbasis + idx] = psi[(j * block_size +st1) * pbasis + idx] * vtot_eig[idx];
            } 
        }

        RmgGemm(trans_a, trans_n, size_row, size_col,  pbasis, alpha, psi+j*block_size*pbasis, pbasis, vpsi, 
                pbasis, beta, global_matrix1, size_row);
        BlockAllreduce((double *)global_matrix1, (size_t)size_row * (size_t)size_col * (size_t)factor , pct.grid_comm);

        for(int jst = 0; jst < size_col; jst++)
        {
            for(int ist = 0; ist < size_row; ist++)
            {
                int idx1 = (ist + j * block_size) + (j * block_size + jst) * num_states;
                int idx2 = (ist + j * block_size) * num_states + (j * block_size + jst);
                Aij[idx1] = global_matrix1[jst * size_row + ist];
                Aij[idx2] = MyConj(global_matrix1[jst * size_row + ist]);

            }
        }
    }

#endif

}

#if CUDA_ENABLED || HIP_ENABLED
void Veff_x_psi(double *psi_dev,  double *work_dev, double *vtot_eig, int pbasis, int num_states)
{
    gpublasStatus_t gstat;
    static double *v_dev;
    if(!v_dev)
    {
        gpuMalloc((void **)&v_dev, pbasis * sizeof(double));
    }
    gpuMemcpy(v_dev, vtot_eig,  pbasis * sizeof(double), gpuMemcpyHostToDevice);
    gstat = gpublasDdgmm(ct.gpublas_handle, GPUBLAS_SIDE_LEFT, pbasis, num_states, 
            (double *)psi_dev, pbasis, (double *)v_dev, 1, (double *)work_dev, pbasis);
    RmgGpuError(__FILE__, __LINE__, gstat, "Error performing gpublasDgmm.");
}
void Veff_x_psi(std::complex<double> *psi_dev,  std::complex<double> *work_dev, double *vtot_eig, int pbasis, int num_states)
{
    gpublasStatus_t gstat;
    static std::complex<double> *v_dev, *vtot_eig_C;
    if(!v_dev)
    {
        gpuMalloc((void **)&v_dev, pbasis * sizeof(std::complex<double>));
        vtot_eig_C = new std::complex<double>[pbasis];
    }
    for(int i = 0; i < pbasis; i++) vtot_eig_C[i] = vtot_eig[i];
    gpuMemcpy(v_dev, vtot_eig_C,  pbasis * sizeof(std::complex<double>), gpuMemcpyHostToDevice);

#if  HIP_ENABLED
    gstat = hipblasZdgmm(ct.gpublas_handle, GPUBLAS_SIDE_LEFT, pbasis, num_states, 
            (hipblasDoubleComplex *)psi_dev, pbasis, (hipblasDoubleComplex *)v_dev, 1, (hipblasDoubleComplex *)work_dev, pbasis);
#endif
#if  CUDA_ENABLED 
    gstat = cublasZdgmm(ct.gpublas_handle, GPUBLAS_SIDE_LEFT, pbasis, num_states, 
            reinterpret_cast<cuDoubleComplex*>(psi_dev), pbasis,
            reinterpret_cast<cuDoubleComplex*>(v_dev), 1,
            reinterpret_cast<cuDoubleComplex*>(work_dev), pbasis);
#endif
    RmgGpuError(__FILE__, __LINE__, gstat, "Error performing gpublasDgmm.");
}
#endif
