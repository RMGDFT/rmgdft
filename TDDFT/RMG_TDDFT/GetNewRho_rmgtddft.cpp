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
#include "Prolong.h"

template void GetNewRho_rmgtddft<double>(Kpoint<double> *,double *rho, double *rho_matrix, int numst, int tddft_start_state);
template void GetNewRho_rmgtddft<std::complex<double> >(Kpoint<std::complex<double>> *, double *rho, std::complex<double> *rho_matrix, int numst, int tddft_start_state);
template <typename KpointType>
void GetNewRho_rmgtddft (Kpoint<KpointType> *kptr, double *rho_k, KpointType *rho_matrix, int numst, int tddft_start_state)
{
    int idx;

    /* for parallel libraries */

    int st1;

    KpointType one = 1.0, zero = 0.0;
    int pbasis = get_P0_BASIS();

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
        KpointType *xpsi = kptr->work_dev;
        RmgGemm ("N", "N", pbasis, numst, numst, one, 
                psi_dev, pbasis, rho_matrix, numst, zero, xpsi, pbasis);
        GpuProductBr(psi_dev, xpsi, rho_temp_dev, numst, pbasis);
        gpuMemcpy(rho_temp, rho_temp_dev,  pbasis * sizeof(double), gpuMemcpyDeviceToHost);
#else
        RmgTimer *RT = new RmgTimer("TDDFT: rho: gemm");
        KpointType *psi = &kptr->orbital_storage[tddft_start_state * pbasis];
        KpointType *xpsi = kptr->work_cpu;
        RmgGemm ("N", "N", pbasis, numst, numst, one, 
                psi, pbasis, rho_matrix, numst, zero, xpsi, pbasis);

        delete RT;
        RT = new RmgTimer("TDDFT: rho: dot");
        for(st1 = 0; st1 < numst; st1++)
            for(idx = 0; idx < pbasis; idx++)
                rho_temp[idx] += std::real(psi[st1 * pbasis + idx] * std::conj(xpsi[st1 * pbasis + idx]));
        delete RT;
#endif
    }


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

    P.prolong(rho_k, rho_temp, dimx, dimy, dimz, half_dimx, half_dimy, half_dimz);

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

