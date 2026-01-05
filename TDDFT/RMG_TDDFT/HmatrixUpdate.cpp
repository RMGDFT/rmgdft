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
#include "rmg_gemm.h"
#include "Subdiag.h"
#include "GpuAlloc.h"
#include "ErrorFuncs.h"
#include "blas.h"
#include "RmgParallelFft.h"

#include "common_prototypes.h"
#include "common_prototypes1.h"
#include "transition.h"
#include "prototypes_tddft.h"
#include "GatherScatter.h"

#if HIP_ENABLED
#include <hip/hip_runtime.h>
#endif

#if CUDA_ENABLED
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cublas_v2.h>
#endif

#if CUDA_ENABLED || HIP_ENABLED
void Veff_x_psi(double *psi_dev,  double *work_dev, double *vtot_psi, int pbasis, int num_states);
void Veff_x_psi(std::complex<double> *psi_dev,  std::complex<double> *work_dev, std::complex<double> *vtot_psi, int pbasis, int num_states);
void Veff_x_psi(float *psi_dev,  float *work_dev, float *vtot_psi, int pbasis, int num_states);
void Veff_x_psi(std::complex<float> *psi_dev,  std::complex<float> *work_dev, std::complex<float> *vtot_psi, int pbasis, int num_states);
#endif

template void HmatrixUpdate<double, double>(Kpoint<double> *, wfobj<double> vtot_psi, wf_spinobj<double> vxc_psi, double *, int tddft_start_state, int num_states);
template void HmatrixUpdate<double, float>(Kpoint<double> *, wfobj<double> vtot_psi, wf_spinobj<double> vxc_psi, double *, int tddft_start_state, int num_states);
template void HmatrixUpdate<std::complex<double>, std::complex<double> >(Kpoint<std::complex<double>> *, wfobj<double> vtot_psi, wf_spinobj<double> vxc_psi, std::complex<double> *, int tddft_start_state, int num_states);
template void HmatrixUpdate<std::complex<double>, std::complex<float> >(Kpoint<std::complex<double>> *, wfobj<double> vtot_psi, wf_spinobj<double> vxc_psi, std::complex<double> *, int tddft_start_state, int num_states);

template <typename KpointType, typename CalType>
void HmatrixUpdate (Kpoint<KpointType> *kptr, wfobj<double> vtot_psi, wf_spinobj<double> vxc_psi, KpointType *Aij, int tddft_start_state, int num_states)
{

    BaseGrid *G = kptr->G;
    Lattice *L = kptr->L;
    int pbasis = kptr->pbasis;
    int pbasis_noncoll = kptr->pbasis * ct.noncoll_factor;
    double vel = L->get_omega() / ((double)(G->get_NX_GRID(1) * G->get_NY_GRID(1) * G->get_NZ_GRID(1)));

    //value type of CalType
    using TypeV = get_scalar_t<CalType>;
    static CalType *global_matrix1;

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
    CalType alpha(vel);
    CalType beta(0.0);
    CalType *psi_dev;
    CalType *work_dev;

    static wfobj<CalType> vtot_psi_caltype;
    static wf_spinobj<CalType> vxc_psi_caltype;
    CopyAndConvert(pbasis, vtot_psi.data(), vtot_psi_caltype.data());
    if(ct.noncoll)
    {
        CopyAndConvert(4*pbasis, vxc_psi.data(), vxc_psi_caltype.data());
    }
    if(typeid(KpointType) == typeid(CalType))
    {
        psi_dev = (CalType *)kptr->psi_dev;
        work_dev = (CalType *)kptr->work_dev;
    }
    else
    {
        psi_dev = (CalType *)kptr->psi_dev_float;
        work_dev = (CalType *)kptr->work_dev_float;
    }

    psi_dev = psi_dev + tddft_start_state * pbasis_noncoll;

    static CalType *v_dev, *vxc_dev;
    static CalType *mat_dev;
    if(!v_dev)
    {
        gpuMalloc((void **)&v_dev, pbasis * sizeof(CalType));
        if(ct.noncoll)
        {
            gpuMalloc((void **)&vxc_dev, 4*pbasis * sizeof(CalType));
        }
    }

    if(ct.noncoll)
    {
        rmg::error(" failure in HmatrixUpdate");
        gpuMemcpy(v_dev, vtot_psi_caltype.data(),  pbasis * sizeof(CalType), gpuMemcpyHostToDevice);
        gpuMemcpy(vxc_dev, vxc_psi_caltype.data(),  4*pbasis * sizeof(CalType), gpuMemcpyHostToDevice);
        Veff_x_psi(psi_dev, work_dev, v_dev, pbasis_noncoll, num_states);
    }
    else
    {
        gpuMemcpy(v_dev, vtot_psi_caltype.data(),  pbasis * sizeof(CalType), gpuMemcpyHostToDevice);
        Veff_x_psi(psi_dev, work_dev, v_dev, pbasis_noncoll, num_states);
    }

    gpublasStatus_t gstat;

    int block_size = std::max(1024, ct.scalapack_block_factor);
    //block_size = num_states;
    int nblock = (num_states + block_size -1)/block_size;


    if(!mat_dev)
    {
        gpuMalloc((void **)&mat_dev, num_states * block_size * sizeof(CalType));
        int retval1 = MPI_Alloc_mem(num_states * block_size * sizeof(CalType) , MPI_INFO_NULL, &global_matrix1);

        if(retval1 != MPI_SUCCESS) {
            rmg::error("Memory allocation failure in HmatrixUpdate");
        }
    }

    // V|psi> is in work_dev now
    // Compute A matrix
    for(int j = 0; j < nblock; j++)
    {
        int size_col = std::min(block_size, num_states - j * block_size);
        int size_row = num_states - j * block_size;
        //if(ct.tddft_floatprecision && 0)
        //{
        //    hipblasStatus_t hipstat;
        //    hipblasOperation_t hip_transA = HIPBLAS_OP_N, hip_transN = HIPBLAS_OP_N;

        //    if(!strcmp(trans_a, "t")) hip_transA = HIPBLAS_OP_T;
        //    if(!strcmp(trans_a, "T")) hip_transA = HIPBLAS_OP_T;
        //    if(!strcmp(trans_a, "c")) hip_transA = HIPBLAS_OP_C;
        //    if(!strcmp(trans_a, "C")) hip_transA = HIPBLAS_OP_C;


        //    if(ct.is_gamma)
        //    {
        //        float alpha_f = std::real(alpha);
        //        float beta_f = 0.0;
        //        hipstat =  hipblasGemmEx(ct.gpublas_handle, hip_transA, hip_transN, size_row, size_col, pbasis, 
        //                &alpha, psi_dev + j*block_size * pbasis, HIP_R_32F, pbasis, 
        //                work_dev + j * block_size * pbasis, HIP_R_32F, pbasis, &beta, 
        //                mat_dev, HIP_R_32F, size_row, HIPBLAS_COMPUTE_64F, HIPBLAS_GEMM_DEFAULT);
        //         // rocblas_gemm_ex(ct.roc_handle, rocblas_operation_transpose, rocblas_operation_none, size_row, size_col, pbasis, 
        //         //       &alpha_f, psi_dev + j*block_size * pbasis, rocblas_datatype_f32_r, pbasis, 
        //         //       work_dev + j * block_size * pbasis,  rocblas_datatype_f32_r, pbasis, &beta_f, 
        //         //       mat_dev1, rocblas_datatype_f32_r, size_row, 
        //         //       mat_dev, rocblas_datatype_f32_r, size_row, 
        //         //       rocblas_datatype_f64_r, rocblas_gemm_algo_standard, 0,0);
        //    }
        //    else
        //    {
        //                  //rocblas_operation_conjugate_
        //        hipstat =  hipblasGemmEx(ct.gpublas_handle, hip_transA, hip_transN, size_row, size_col, pbasis, 
        //                &alpha, psi_dev + j*block_size * pbasis, HIP_C_32F, pbasis, 
        //                work_dev + j * block_size * pbasis, HIP_C_32F, pbasis, &beta, 
        //                mat_dev, HIP_C_64F, size_row, HIPBLAS_COMPUTE_64F, HIPBLAS_GEMM_DEFAULT);
        //    }
        //}
        //else
        {
            rmg::gemm(trans_a, trans_n, size_row, size_col,  pbasis_noncoll, alpha, (psi_dev+ j*block_size*pbasis_noncoll), pbasis_noncoll, (work_dev + j * block_size * pbasis_noncoll), 
                    pbasis_noncoll, beta, mat_dev, size_row);
        }
        gpuMemcpy(global_matrix1, mat_dev,  (size_t)size_row * (size_t)size_col * sizeof(CalType), gpuMemcpyDeviceToHost);

        BlockAllreduce(global_matrix1, (size_t)size_row * (size_t)size_col, pct.grid_comm);

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

    if(typeid(KpointType) != typeid(CalType))
    {
        rmg::error("float precision not programmed with cpu  failure in HmatrixUpdate");
    }
    KpointType alpha(vel);
    KpointType beta(0.0);
    int block_size = ct.state_block_size;
    //block_size = num_states;
    int nblock = (num_states + block_size -1)/block_size;
    // First time through allocate pinned memory for global_matrix1
    if(!global_matrix1) {

        int retval1 = MPI_Alloc_mem(num_states * block_size * sizeof(KpointType) , MPI_INFO_NULL, &global_matrix1);

        if(retval1 != MPI_SUCCESS) {
            rmg::error("Memory allocation failure in HmatrixUpdate");
        }

    }

    // V|psi> is in tmp_arrayT
    KpointType *psi = kptr->orbital_storage + tddft_start_state * pbasis_noncoll;
    KpointType *vpsi = &kptr->orbital_storage[kptr->nstates * pbasis_noncoll];  // use the memory of psi extra 3* state_block_size.

    for(int j = 0; j < nblock; j++)
    {
        int size_col = std::min(block_size, num_states - j * block_size);
        int size_row = num_states - j * block_size;

        for (int st1 = 0; st1 < size_col; st1++)
        {
            for(int idx = 0; idx <pbasis_noncoll; idx++)
            {
                vpsi[st1 * pbasis_noncoll + idx] = psi[(j * block_size +st1) * pbasis_noncoll + idx] * vtot_psi[idx];
            } 
        }
        // >>>>????
        if(ct.noncoll)
        {
            for (int st1 = 0; st1 < size_col; st1++)
            {
                std::complex<double> *vpsi_C = (std::complex<double> *)&vpsi[st1 * pbasis_noncoll];
                std::complex<double> *psi_C = (std::complex<double> *)&psi[(j * block_size +st1) * pbasis_noncoll];

                for(int idx = 0; idx <pbasis; idx++)
                {
                      vpsi_C[idx] += psi_C[idx] * std::complex<double>(vxc_psi.cz[idx], 0.0);
                      vpsi_C[idx] += psi_C[idx+pbasis] * std::complex<double>(vxc_psi.cx[idx], vxc_psi.cy[idx]);
                      vpsi_C[idx + pbasis] += - psi_C[idx + pbasis] * std::complex<double>(vxc_psi.cz[idx], 0.0);
                      vpsi_C[idx + pbasis] += psi_C[idx] * std::complex<double>(vxc_psi.cx[idx], -vxc_psi.cy[idx]);

                } 
            }
        }

        rmg::gemm(trans_a, trans_n, size_row, size_col,  pbasis_noncoll, alpha, psi+j*block_size*pbasis_noncoll, pbasis_noncoll, vpsi, 
                pbasis_noncoll, beta, (KpointType *)global_matrix1, size_row);
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
void Veff_x_psi(double *psi_dev,  double *work_dev, double *v_dev, int pbasis, int num_states)
{
    gpublasStatus_t gstat;
    gstat = gpublasDdgmm(ct.gpublas_handle, GPUBLAS_SIDE_LEFT, pbasis, num_states, 
            (double *)psi_dev, pbasis, (double *)v_dev, 1, (double *)work_dev, pbasis);
    RmgGpuError(__FILE__, __LINE__, gstat, "Error performing gpublasDgmm.");
}
void Veff_x_psi(float *psi_dev,  float *work_dev, float *v_dev, int pbasis, int num_states)
{
    gpublasStatus_t gstat;
    gstat = gpublasSdgmm(ct.gpublas_handle, GPUBLAS_SIDE_LEFT, pbasis, num_states, 
            (float *)psi_dev, pbasis, (float *)v_dev, 1, (float *)work_dev, pbasis);
    RmgGpuError(__FILE__, __LINE__, gstat, "Error performing gpublasDgmm.");
}
void Veff_x_psi(std::complex<double> *psi_dev,  std::complex<double> *work_dev, std::complex<double> *v_dev, int pbasis, int num_states)
{
    gpublasStatus_t gstat;

#if  HIP_ENABLED
    gstat = hipblasZdgmm(ct.gpublas_handle, GPUBLAS_SIDE_LEFT, pbasis, num_states, 
            (hipDoubleComplex *)psi_dev, pbasis, (hipDoubleComplex *)v_dev, 1, (hipDoubleComplex *)work_dev, pbasis);
#endif
#if  CUDA_ENABLED 
    gstat = cublasZdgmm(ct.gpublas_handle, GPUBLAS_SIDE_LEFT, pbasis, num_states, 
            reinterpret_cast<cuDoubleComplex*>(psi_dev), pbasis,
            reinterpret_cast<cuDoubleComplex*>(v_dev), 1,
            reinterpret_cast<cuDoubleComplex*>(work_dev), pbasis);
#endif
    RmgGpuError(__FILE__, __LINE__, gstat, "Error performing gpublasDgmm.");
}
void Veff_x_psi(std::complex<float> *psi_dev,  std::complex<float> *work_dev, std::complex<float> *v_dev, int pbasis, int num_states)
{
    gpublasStatus_t gstat;

#if  HIP_ENABLED
    gstat = hipblasCdgmm(ct.gpublas_handle, GPUBLAS_SIDE_LEFT, pbasis, num_states, 
            (hipFloatComplex *)psi_dev, pbasis, (hipFloatComplex *)v_dev, 1, (hipFloatComplex *)work_dev, pbasis);
#endif
#if  CUDA_ENABLED 
    gstat = cublasCdgmm(ct.gpublas_handle, GPUBLAS_SIDE_LEFT, pbasis, num_states, 
            reinterpret_cast<cuFloatComplex*>(psi_dev), pbasis,
            reinterpret_cast<cuFloatComplex*>(v_dev), 1,
            reinterpret_cast<cuFloatComplex*>(work_dev), pbasis);
#endif
    RmgGpuError(__FILE__, __LINE__, gstat, "Error performing gpublasDgmm.");
}
#endif


