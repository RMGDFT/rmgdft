#include "make_conf.h"
#include "const.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "Kpoint.h"
#include "transition.h"

template void ReinitIonicPotentials<double>(Kpoint<double> **, double *, double *, double *);
template void ReinitIonicPotentials<std::complex<double> >(Kpoint<std::complex<double>> **, double *, double *, double *);

template <typename KpointType>
void ReinitIonicPotentials (Kpoint<KpointType> **Kptr, double * vnuc, double * rhocore, double * rhoc)
{

    KpointType ZERO_t(0.0);
    int pbasis = Kptr[0]->pbasis;

    /* Update items that change when the ionic coordinates change */
    init_nuc (vnuc, rhoc, rhocore);
    get_QI ();
    GetNlop(Kptr);

    // Number of total projectors required is computed in GetNlop so we allocate per
    // k-point storage for the weights here.
    for(int kpt=0;kpt < ct.num_kpts;kpt++) {

        if(ct.is_gamma) {

            // Identical for gamma point
            Kptr[kpt]->nl_weight = (KpointType *)pct.weight;
            Kptr[kpt]->nl_Bweight = (KpointType *)pct.Bweight;

        }
        else {

            // Release old memory storage for weights
            if(Kptr[kpt]->nl_weight != NULL) delete [] Kptr[kpt]->nl_weight;
            if(Kptr[kpt]->nl_Bweight != NULL) delete [] Kptr[kpt]->nl_Bweight;

            // Allocate new storage
            Kptr[kpt]->nl_weight = new KpointType[pct.num_tot_proj * pbasis];
            Kptr[kpt]->nl_Bweight = new KpointType[pct.num_tot_proj * pbasis];
            for(int idx = 0;idx < pct.num_tot_proj * pbasis;idx++) Kptr[kpt]->nl_weight[idx] = ZERO_t;    
            for(int idx = 0;idx < pct.num_tot_proj * pbasis;idx++) Kptr[kpt]->nl_Bweight[idx] = ZERO_t;    
        }

    }


    /*Other things that need to be recalculated when ionic positions change */
    GetWeight (Kptr);
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

