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
#include "RmgGemm.h"
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


template void Subdiag_Lapack<double> (Kpoint<double> *kptr, double *Aij, double *Bij, double *Sij, double *eigs, double *eigvectors);
template void Subdiag_Lapack<std::complex<double> > (Kpoint<std::complex<double>> *kptr, std::complex<double> *Aij, std::complex<double> *Bij, std::complex<double> *Sij, double *eigs, std::complex<double> *eigvectors);

// eigvectors holds Bij on input and the eigenvectors of the matrix on output
template <typename KpointType>
void Subdiag_Lapack (Kpoint<KpointType> *kptr, KpointType *Aij, KpointType *Bij, KpointType *Sij, double *eigs, KpointType *eigvectors)
{

    int num_states = kptr->nstates;
    int ione = 1;

    KpointType ZERO_t(0.0);
    KpointType ONE_t(1.0);
    KpointType *NULLptr = NULL;
    KpointType *Cij = new KpointType[num_states * num_states];

    // Create unitary matrix
    for (int idx = 0; idx < num_states * num_states; idx++) {
        Cij[idx] = ZERO_t;
    }


    for (int idx = 0; idx < num_states; idx++) {
        Cij[idx * num_states + idx] = ONE_t;
    }


    // Invert Bij

    int *ipiv = new int[2*num_states];
    for(int idx=0;idx < num_states;idx++) ipiv[idx] = 0;
    int info = 0;
    if(ct.is_gamma) {

        // Inverse of B should be in Cij
        RmgTimer *RT1 = new RmgTimer("Diagonalization: Invert Bij");
        dgesv (&num_states, &num_states, (double *)eigvectors, &num_states, ipiv, (double *)Cij, &num_states, &info);
        delete(RT1);

    }
    else {

        // Inverse of B should be in Cij
        RmgTimer *RT1 = new RmgTimer("Diagonalization: Invert Bij");
        zgesv (&num_states, &num_states, (double *)eigvectors, &num_states, ipiv, (double *)Cij, &num_states, &info);
        delete(RT1);

    }
    if (info) {
        rmg_printf ("\n PE %d: p{d,z}gesv failed, info is %d", pct.gridpe, info);
        rmg_error_handler (__FILE__, __LINE__, " p{d,z}gesv failed");
    }
    delete [] ipiv;


    /*Multiply inverse of B and and A */
    /*B^-1*A */
    KpointType alpha(1.0);
    KpointType beta(0.0);;

    RmgTimer *RT1 = new RmgTimer("Diagonalization: matrix setup");
    RmgGemm ("n", "n", num_states, num_states, num_states, alpha,
                    Cij, num_states, Aij, num_states, beta, Bij,
                    num_states, NULLptr, NULLptr, NULLptr,
                    true, true, true, true);

    /*Multiply the result with Sij, result is in Cij */
    RmgGemm ("n", "n", num_states, num_states, num_states, alpha,
                    Sij, num_states, Bij, num_states, beta, Cij,
                    num_states, NULLptr, NULLptr, NULLptr,
                    true, true, true, true);
    delete(RT1);


    RT1 = new RmgTimer("Diagonalization: lapack");
    int *ifail = new int[num_states];
    int lwork = 2 * num_states * num_states + 6 * num_states + 2;
    int liwork = 6*num_states;
    int eigs_found;
    double *work2 = new double[2*lwork];
    int *iwork = new int[liwork];
    double vx = 0.0;
    double tol = 1e-15;

    if(ct.is_gamma) {

        dsygvx (&ione, "v", "A", "l", &num_states, (double *)Cij, &num_states, (double *)Sij, &num_states,
                        &vx, &vx, &ione, &ione,  &tol, &eigs_found, eigs, (double *)eigvectors, &num_states, work2,
                        &lwork, iwork, ifail, &info);

    }
    else {

        double *rwork = new double[8 * num_states];
        zhegvx (&ione, "v", "A", "l", &num_states, (double *)Cij, &num_states, (double *)Sij, &num_states,
                        &vx, &vx, &ione, &ione,  &tol, &eigs_found, eigs, (double *)eigvectors, &num_states, work2,
                        &lwork, rwork, iwork, ifail, &info);
        delete [] rwork;

    }

    delete [] iwork;
    delete [] work2;
    delete [] ifail;

    if (info) {
        rmg_printf ("\n Lapack eigensolver failed, info is %d", info);
        rmg_error_handler (__FILE__, __LINE__, "Lapack eigensolver failed");
    }

    delete(RT1);


    delete [] Cij;
}
