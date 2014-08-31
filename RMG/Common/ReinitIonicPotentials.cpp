#include "make_conf.h"
#include "const.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "Kpoint.h"
#include "transition.h"

template void ReinitIonicPotentials<double>(Kpoint<double> **, double *, double *, double *);
template void ReinitIonicPotentials<std::complex<double> >(Kpoint<std::complex<double>> **, double *, double *, double *);

template <typename KpointType>
void ReinitIonicPotentials (Kpoint<KpointType> **kptr, double * vnuc, double * rhocore, double * rhoc)
{

    /* Update items that change when the ionic coordinates change */
    init_nuc (vnuc, rhoc, rhocore);
    get_QI ();
    get_nlop ();

    /*Other things that need to be recalculated when ionic positions change */
    get_weight ();
    get_qqq ();

#if GPU_ENABLED
    cublasStatus_t custat;

    // If gpu weight buffer has not been setup yet or size has changed must take care of allocation
    if((ct.gpu_weight == NULL) || (gpu_weight_alloc < (get_P0_BASIS() * pct.num_tot_proj * sizeof(double)))) {
        if(ct.gpu_weight != NULL) {
            cudaFree(ct.gpu_weight);
            ct.gpu_weight = NULL;
        }
        if(pct.num_tot_proj) {
            if( cudaSuccess != cudaMalloc((void **)&ct.gpu_weight , get_P0_BASIS() * pct.num_tot_proj * sizeof(double) ))
                rmg_error_handler(__FILE__, __LINE__, "cudaMalloc failed for: gpu_weight\n");
        }
        gpu_weight_alloc = get_P0_BASIS() * pct.num_tot_proj * sizeof(double);
    }

    // Transfer copy of weights to GPU
    custat = cublasSetVector( get_P0_BASIS() * pct.num_tot_proj, sizeof( double ), pct.weight, 1, ct.gpu_weight, 1 );
    RmgCudaError(__FILE__, __LINE__, custat, "Problem transferring weight mmatrix from GPU to system memory.");

    if((ct.gpu_Bweight == NULL) || (gpu_weight_alloc < (get_P0_BASIS() * pct.num_tot_proj * sizeof(double)))) {
        if(ct.gpu_Bweight != NULL) {
            cudaFree(ct.gpu_Bweight);
            ct.gpu_Bweight = NULL;
        }
        if(pct.num_tot_proj) {
            if( cudaSuccess != cudaMalloc((void **)&ct.gpu_Bweight , get_P0_BASIS() * pct.num_tot_proj * sizeof(double) ))
                rmg_error_handler(__FILE__, __LINE__, "cudaMalloc failed for: gpu_Bweight\n");
        }
    }

    // Transfer copy of weights to GPU
    custat = cublasSetVector( get_P0_BASIS() * pct.num_tot_proj, sizeof( double ), pct.Bweight, 1, ct.gpu_Bweight, 1 );
    RmgCudaError(__FILE__, __LINE__, custat, "Problem transferring Bweight mmatrix from GPU to system memory.");

#endif

#if 0
    if (!verify ("calculation_mode", "Band Structure Only"))
    {
        betaxpsi (states);
        mix_betaxpsi(0);
    }
#endif

}                               /* end reinit_ionic_pp */

