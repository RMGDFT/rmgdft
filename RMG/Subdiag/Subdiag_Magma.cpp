#include <complex>
#include "FiniteDiff.h"
#include "const.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "rmg_error.h"
#include "RmgTimer.h"
#include "GlobalSums.h"
#include "Kpoint.h"
#include "Subdiag.h"
#include "blas.h"

#include "prototypes.h"
#include "common_prototypes.h"
#include "common_prototypes1.h"
#include "transition.h"

#if GPU_ENABLED
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cublas_v2.h>
#endif


template void Subdiag_Magma<double> (Kpoint<double> *kptr, double *Aij, double *Bij, double *Sij, double *eigs, double *eigvectors);
template void Subdiag_Magma<std::complex<double> > (Kpoint<std::complex<double>> *kptr, std::complex<double> *Aij, std::complex<double> *Bij, std::complex<double> *Sij, double *eigs, std::complex<double> *eigvectors);

// eigvectors holds Bij on input and the eigenvectors of the matrix on output
template <typename KpointType>
void Subdiag_Magma (Kpoint<KpointType> *kptr, KpointType *Aij, KpointType *Bij, KpointType *Sij, double *eigs, KpointType *eigvectors)
{

#if !MAGMA_LIBS
    rmg_printf("This version of RMG was not built with MAGMA support. Redirecting to LAPACK.");
    Subdiag_Lapack(kptr, Aij, Bij, Sij, eigs, eigvectors);
    return;
#endif

#if GPU_ENABLED
    KpointType ZERO_t(0.0);
    KpointType ONE_t(1.0);
    int num_states = kptr->nstates;
    int ione = 1;
    int factor = 1;
    if(!ct.is_gamma) factor=2;

    cublasStatus_t custat;
    cublasOperation_t cu_transT = CUBLAS_OP_T, cu_transN = CUBLAS_OP_N;



    KpointType *Cij = new KpointType[num_states * num_states];

    // Create unitary matrix
    for (int idx = 0; idx < num_states * num_states; idx++) {
        Cij[idx] = ZERO_t;
    }


    for (int idx = 0; idx < num_states; idx++) {
        Cij[idx * num_states + idx] = ONE_t;
    }
    custat = cublasSetVector(num_states * num_states , sizeof( rmg_double_t ), distCij, ione, gpuCij, ione );
#endif

     
}

