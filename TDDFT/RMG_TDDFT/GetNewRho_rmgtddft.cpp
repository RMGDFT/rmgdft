/************************** SVN Revision Information **************************
 **    $Id: get_new_rho_local.c 3140 2015-08-06 15:48:24Z luw $    **
 ******************************************************************************/

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include "const.h"
#include "RmgTimer.h"
#include "RmgException.h"
#include "rmgtypedefs.h"
#include "params.h"
#include "typedefs.h"
#include "rmg_error.h"
#include "blas.h"
#include "init_var.h"
#include "transition.h"
#include "prototypes_on.h"
#include "Kbpsi.h"
#include "Gpufuncs.h"
#include "Kpoint.h"


#include "../Headers/common_prototypes.h"
#include "../Headers/common_prototypes1.h"
#include "prototypes_tddft.h"
#include "RmgParallelFft.h"
#include "rmg_gemm.h"
#include "blas_driver.h"
#include "Prolong.h"
#include "GatherScatter.h"

#include <type_traits>

// 1. Primary template: handles scalar types (float, double, etc.)
template <typename T>
struct get_scalar {
    using type = T;
};

// 2. Partial specialization: handles std::complex types
template <typename T>
struct get_scalar<std::complex<T>> {
    using type = T;
};

// 3. Helper alias (Optional but recommended for modern C++)
template <typename T>
using get_scalar_t = typename get_scalar<T>::type;

// --- Usage ---
//using TypeA = get_scalar_t<float>;                // Result: float
//using TypeB = get_scalar_t<std::complex<float>>;  // Result: float

template void GetNewRho_rmgtddft<double, double>(Kpoint<double> *,spinobj<double> &rho, double *rho_matrix, int numst, int tddft_start_state, double *rho_matrix_caltype);
template void GetNewRho_rmgtddft<double, float>(Kpoint<double> *,spinobj<double> &rho, double *rho_matrix, int numst, int tddft_start_state, float *rho_matrix_caltype);
template void GetNewRho_rmgtddft<std::complex<double>, std::complex<double> >(Kpoint<std::complex<double>> *, spinobj<double> &rho, std::complex<double> *rho_matrix, int numst, int tddft_start_state, std::complex<double> *rho_matrix_caltype);
template void GetNewRho_rmgtddft<std::complex<double>, std::complex<float> >(Kpoint<std::complex<double>> *, spinobj<double> &rho, std::complex<double> *rho_matrix, int numst, int tddft_start_state, std::complex<float> *rho_matrix_caltype);
template <typename KpointType, typename CalType>
void GetNewRho_rmgtddft (Kpoint<KpointType> *kptr, spinobj<double> &rho_k, KpointType *rho_matrix, int numst, int tddft_start_state, CalType *rho_matrix_caltype)
{
    int idx;

    /* for parallel libraries */

    int st1;

    using TypeV = get_scalar_t<CalType>;                // Result: float
    int pbasis = get_P0_BASIS();
    int pbasis_noncoll = pbasis * ct.noncoll_factor;
    // in the noncollinear case, noncoll_factor = 2, and number of rho components = 4
    int n_rho = ct.noncoll_factor * ct.noncoll_factor;

    if(!ct.norm_conserving_pp) {
        rmg_error_handler (__FILE__, __LINE__, "\n tddft not programed for ultrasoft \n");
    }
    if(ct.noncoll)
    {
        rmg_error_handler (__FILE__, __LINE__, "\n tddft not programed for noncoll \n");
    }

    for (int istate = 0; istate < numst; istate++)
    {
        rho_matrix[istate * numst + istate] -=
            kptr->Kstates[istate + tddft_start_state].occupation[0];

    }

    CopyAndConvert(numst*numst, rho_matrix, rho_matrix_caltype);
#if CUDA_ENABLED || HIP_ENABLED 
    // xpsi is a device buffer in this case and GpuProductBr is a GPU functions to do
    // the reduction over numst.
    
    CalType one = 1.0, zero = 0.0;
    TypeV *rho_temp, *rho_temp_dev;


    rho_temp = (TypeV *)GpuMallocHost(pbasis*n_rho * sizeof(TypeV));
    gpuMalloc((void **)&rho_temp_dev, pbasis*n_rho * sizeof(TypeV));

    CalType *psi_dev;
    CalType *xpsi;
    if(typeid(KpointType) == typeid(CalType))
    {
        psi_dev = (CalType *)kptr->psi_dev;
        xpsi = (CalType *)kptr->work_dev;
    }
    else
    {
        psi_dev = (CalType *)kptr->psi_dev_float;
        xpsi = (CalType *)kptr->work_dev_float;
    }

    psi_dev = psi_dev + tddft_start_state * pbasis_noncoll;

    rmg::gemm ("N", "N", pbasis_noncoll, numst, numst, one, 
            psi_dev, pbasis_noncoll, rho_matrix_caltype, numst, zero, xpsi, pbasis_noncoll);

    GpuProductBr(psi_dev, xpsi, rho_temp_dev, numst, pbasis);
    gpuMemcpy(rho_temp, rho_temp_dev,  n_rho * pbasis * sizeof(TypeV), gpuMemcpyDeviceToHost);
#else
    if(typeid(KpointType) != typeid(CalType))
    {
        rmg_error_handler (__FILE__, __LINE__, "\n float precision not for CPU \n");
    }
    KpointType one = 1.0, zero = 0.0;
    RmgTimer *RT = new RmgTimer("TDDFT: rho: gemm");

    double *rho_temp = new double[n_rho * pbasis]();

    KpointType *psi = &kptr->orbital_storage[tddft_start_state * pbasis_noncoll];
    KpointType *xpsi = kptr->work_cpu;
    rmg::gemm ("N", "N", pbasis_noncoll, numst, numst, one, 
            psi, pbasis_noncoll, rho_matrix, numst, zero, xpsi, pbasis_noncoll);

    delete RT;
    RT = new RmgTimer("TDDFT: rho: dot");
    if(!ct.noncoll)
    {
        for(st1 = 0; st1 < numst; st1++)
        {
            for(idx = 0; idx < pbasis; idx++)
            {
                rho_temp[idx] += std::real(psi[st1 * pbasis_noncoll + idx] * std::conj(xpsi[st1 * pbasis_noncoll + idx]));
            }
        }
    }
    else
    {
        for(st1 = 0; st1 < numst; st1++)
        {

            KpointType *one_psi = &psi[st1 * pbasis_noncoll];
            KpointType *one_xpsi = &xpsi[st1 * pbasis_noncoll];
            for(idx = 0; idx < pbasis; idx++)
            {
                double rho_up = std::real(one_psi[idx] * std::conj(one_xpsi[idx]));
                double rho_dn = std::real(one_psi[idx + pbasis] * std::conj(one_xpsi[idx + pbasis]));
                std::complex<double> psiud = one_psi[idx] * std::conj(one_xpsi[idx + pbasis]);
                std::complex<double> psidu = one_psi[idx+pbasis] * std::conj(one_xpsi[idx]);
                rho_temp[idx] += rho_up + rho_dn;
                rho_temp[idx + pbasis] += std::real(psiud + psidu);
                rho_temp[idx + pbasis *2] += std::imag(psiud + psidu);
                rho_temp[idx + pbasis *3] += rho_up - rho_dn;
            }
        }
    }

    delete RT;
#endif


    /* Interpolate onto fine grid, result will be stored in rho*/
    RmgTimer *RT1 = new RmgTimer("TDDFT: rho: interp");
    int ratio = Rmg_G->default_FG_RATIO;
    int dimx = Rmg_G->get_PX0_GRID(ratio);
    int dimy = Rmg_G->get_PY0_GRID(ratio);
    int dimz = Rmg_G->get_PZ0_GRID(ratio);
    int half_dimx = Rmg_G->get_PX0_GRID(1);
    int half_dimy = Rmg_G->get_PY0_GRID(1);
    int half_dimz = Rmg_G->get_PZ0_GRID(1);



    static Prolong P(ratio, ct.prolong_order, ct.cmix, *Rmg_T,  Rmg_L, *Rmg_G);

    spinobj<double> rho_k_double;
    int fpbasis = dimx * dimy * dimz;
    CopyAndConvert(n_rho * pbasis, rho_temp, rho_k_double.data());
    for(int irho = 0; irho< n_rho; irho++)
    {
        P.prolong(rho_k.data()+irho*fpbasis, rho_k_double.data() + irho*pbasis, dimx, dimy, dimz, half_dimx, half_dimy, half_dimz);
        
    }

    delete RT1;

#if CUDA_ENABLED || HIP_ENABLED 
    gpuFree(rho_temp_dev);
    GpuFreeHost(rho_temp);
#else
    delete [] rho_temp;
#endif
    for (int istate = 0; istate < numst; istate++)
    {
        rho_matrix[istate * numst + istate] +=
            kptr->Kstates[istate + tddft_start_state].occupation[0];

    }
}

