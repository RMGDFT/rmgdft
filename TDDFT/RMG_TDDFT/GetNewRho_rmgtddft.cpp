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
#include "RmgGemm.h"
#include "blas_driver.h"

template void GetNewRho_rmgtddft<double>(Kpoint<double> *,double *rho, double *rho_matrix, int numst, int tddft_start_state, double *rho_ground);
template void GetNewRho_rmgtddft<std::complex<double> >(Kpoint<std::complex<double>> *, double *rho, std::complex<double> *rho_matrix, int numst, int tddft_start_state, double *rho_ground);
template <typename KpointType>
void GetNewRho_rmgtddft (Kpoint<KpointType> *kptr, double *rho, KpointType *rho_matrix, int numst, int tddft_start_state, double *rho_ground)
{
    int idx;

    /* for parallel libraries */

    int st1;

    KpointType one = 1.0, zero = 0.0;
    int pbasis = get_P0_BASIS();
    
    if(!ct.norm_conserving_pp) {
        rmg_error_handler (__FILE__, __LINE__, "\n tddft not programed for ultrasoft \n");
    }
    if(ct.num_kpts > 1)
    {
        rmg_error_handler (__FILE__, __LINE__, "\n tddft not programed for kpoint \n");
    } 
    if(ct.nspin > 1)
    {
        rmg_error_handler (__FILE__, __LINE__, "\n tddft not programed for spin-polarized \n");
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

#if CUDA_ENABLED || HIP_ENABLED 
    double *rho_temp, *rho_temp_dev;
    rho_temp = (double *)GpuMallocHost(pbasis * sizeof(double));
    gpuMalloc((void **)&rho_temp_dev, pbasis * sizeof(double));
#else
    double *rho_temp = new double[pbasis];
    for(idx = 0; idx < pbasis; idx++)rho_temp[idx] = 0.0;
#endif


    if(numst > 0)
    {
#if CUDA_ENABLED || HIP_ENABLED 
        // xpsi is a device buffer in this case and GpuProductBr is a GPU functions to do
        // the reduction over numst.
        KpointType *psi_dev = &kptr->psi_dev[tddft_start_state * pbasis];
        KpointTypedouble *xpsi = kptr->work_dev;
        RmgGemm ("N", "N", pbasis, numst, numst, one, 
                psi_dev, pbasis, rho_matrix, numst, zero, xpsi, pbasis);
        GpuProductBr(psi_dev, xpsi, rho_temp_dev, numst, pbasis);
        gpuMemcpy(rho_temp, rho_temp_dev,  pbasis * sizeof(double), gpuMemcpyDeviceToHost);
#else
        KpointType *psi = &kptr->orbital_storage[tddft_start_state * pbasis];
        KpointType *xpsi = kptr->work_cpu;
        RmgGemm ("N", "N", pbasis, numst, numst, one, 
                psi, pbasis, rho_matrix, numst, zero, xpsi, pbasis);

        for(st1 = 0; st1 < numst; st1++)
            for(idx = 0; idx < pbasis; idx++)
                rho_temp[idx] += std::real(psi[st1 * pbasis + idx] * xpsi[st1 * pbasis + idx]);
#endif
    }


    /* Interpolate onto fine grid, result will be stored in rho*/
    switch (ct.interp_flag)
    {
        case CUBIC_POLYNOMIAL_INTERPOLATION:
            pack_rho_ctof (rho_temp, rho);
            break;
        case PROLONG_INTERPOLATION:
            mg_prolong_MAX10 (rho, rho_temp, get_FPX0_GRID(), get_FPY0_GRID(), get_FPZ0_GRID(), get_PX0_GRID(), get_PY0_GRID(), get_PZ0_GRID(), get_FG_RATIO(), 6);
            break;
        case FFT_INTERPOLATION:
            FftInterpolation (*Rmg_G, rho_temp, rho, Rmg_G->default_FG_RATIO, ct.sqrt_interpolation);
            break;

        default:

            //Dprintf ("charge interpolation is set to %d", ct.interp_flag);
            rmg_error_handler (__FILE__, __LINE__, "ct.interp_flag is set to an invalid value.");


    }



    /* Check total charge. */
    ct.tcharge = ZERO;
    int FP0_BASIS = Rmg_G->get_P0_BASIS(Rmg_G->default_FG_RATIO);
    for (int idx = 0; idx < FP0_BASIS; idx++)
        rho[idx] += rho_ground[idx];
    for (int idx = 0; idx < FP0_BASIS; idx++)
        ct.tcharge += rho[idx];

    /* ct.tcharge = real_sum_all (ct.tcharge); */
    ct.tcharge = real_sum_all (ct.tcharge, pct.img_comm);
    ct.tcharge = ct.tcharge * get_vel_f();

    /* Renormalize charge, there could be some discrpancy because of interpolation */
    double t1 = ct.nel / ct.tcharge;
    if(std::abs(t1-1) > 1.0e-6) rmg_printf ("normalization constant-1 for new charge in tddft is %e\n", t1-1);
    for(int i = 0;i < FP0_BASIS;i++) rho[i] *= t1;

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

