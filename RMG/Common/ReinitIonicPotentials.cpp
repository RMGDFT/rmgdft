#include "make_conf.h"
#include "const.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "Kpoint.h"
#include "ErrorFuncs.h"
#include "transition.h"

template void ReinitIonicPotentials<double>(Kpoint<double> **, double *, double *, double *);
template void ReinitIonicPotentials<std::complex<double> >(Kpoint<std::complex<double>> **, double *, double *, double *);

template <typename KpointType>
void ReinitIonicPotentials (Kpoint<KpointType> **Kptr, double * vnuc, double * rhocore, double * rhoc)
{

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
            if(pct.num_tot_proj) {
                Kptr[kpt]->nl_weight = new KpointType[pct.num_tot_proj * pbasis]();
                Kptr[kpt]->nl_Bweight = new KpointType[pct.num_tot_proj * pbasis]();
            }
        }

#if GPU_ENABLED
        cudaError_t cuerr;

        if(Kptr[kpt]->nl_weight_gpu != NULL) cudaFree(Kptr[kpt]->nl_weight_gpu);
        if(Kptr[kpt]->nl_Bweight_gpu != NULL) cudaFree(Kptr[kpt]->nl_Bweight_gpu);

        // Allocate new storage
        if(pct.num_tot_proj) {

            cuerr = cudaMalloc((void **)&Kptr[kpt]->nl_weight_gpu , pbasis * pct.num_tot_proj * sizeof(KpointType) );
            if(cuerr != cudaSuccess)
                RmgCudaError(__FILE__, __LINE__, cuerr, "GPU memory allocation error");

            cuerr = cudaMalloc((void **)&Kptr[kpt]->nl_Bweight_gpu , pbasis * pct.num_tot_proj * sizeof(KpointType) );
            if(cuerr != cudaSuccess)
                RmgCudaError(__FILE__, __LINE__, cuerr, "GPU memory allocation error");

        }
#endif

    } // end loop over kpts


    /*Other things that need to be recalculated when ionic positions change */
    GetWeight (Kptr);
    get_qqq ();


#if GPU_ENABLED
    cublasStatus_t custat;
    for(int kpt=0;kpt < ct.num_kpts;kpt++) {

        if(pct.num_tot_proj) {

            // Transfer copy of weights to GPU
            custat = cublasSetVector( pbasis * pct.num_tot_proj, sizeof( KpointType ), Kptr[kpt]->nl_weight, 1, Kptr[kpt]->nl_weight_gpu, 1 );
            if(custat != CUBLAS_STATUS_SUCCESS)
                rmg_error_handler(__FILE__, __LINE__, "Problem transferring non-local weight matrix from system memory to GPU.");

            custat = cublasSetVector( pbasis * pct.num_tot_proj, sizeof( KpointType ), Kptr[kpt]->nl_Bweight, 1, Kptr[kpt]->nl_Bweight_gpu, 1 );
            if(custat != CUBLAS_STATUS_SUCCESS)
                rmg_error_handler(__FILE__, __LINE__, "Problem transferring non-local weight matrix from system memory to GPU.");

        }

    }

#endif

#if 0
    if (!verify ("calculation_mode", "Band Structure Only"))
    {
        betaxpsi (states);
        mix_betaxpsi(0);
    }
#endif

}                               /* end reinit_ionic_pp */

