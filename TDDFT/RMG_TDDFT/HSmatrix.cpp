#include <complex>
#include "FiniteDiff.h"
#include "const.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "rmg_error.h"
#include "rmgthreads.h"
#include "RmgTimer.h"
#include "RmgThread.h"
#include "rmg_reduce.h"
#include "Kpoint.h"
#include "rmg_gemm.h"
#include "Gpufuncs.h"
#include "Subdiag.h"
#include "GpuAlloc.h"

#include "blas.h"
#include "Solvers.h"
#include "blas_driver.h"

#include "common_prototypes.h"
#include "common_prototypes1.h"
#include "transition.h"
#include "prototypes_tddft.h"

#if CUDA_ENABLED
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cublas_v2.h>
#endif


template void HSmatrix<double>(Kpoint<double> *, double *, double *, double *, double *);
template void HSmatrix<std::complex<double> >(Kpoint<std::complex<double>> *, double *,double *,  std::complex<double> *, std::complex<double> *); 

    template <typename KpointType>
void HSmatrix (Kpoint<KpointType> *kptr, double *vtot_eig,double *vxc_psi,  KpointType *Hmat, KpointType *Smat) 
{

    double vel = kptr->L->get_omega() / ((double)(kptr->G->get_NX_GRID(1) * kptr->G->get_NY_GRID(1) * kptr->G->get_NZ_GRID(1)));

    int nstates = kptr->nstates;
    MPI_Comm grid_comm = kptr->grid_comm;
    KpointType *orbital_storage = kptr->orbital_storage;
    KpointType *prev_orbital_storage = kptr->prev_orbital_storage;
    int pbasis_noncoll = kptr->pbasis * ct.noncoll_factor;

    // For MPI routines
    int factor = 1;
    if(!ct.is_gamma) factor = 2;
    
    // State array is 4 * the number of states in length but memory above
    // the first set of nstates is unused in this routine so we can use it
    // as temporary space.
    KpointType *tmp_arrayT = kptr->Kstates[0].psi;
    tmp_arrayT += kptr->nstates * pbasis_noncoll;


    char *trans_t = "t";
    char *trans_n = "n";
    char *trans_c = "c";
    char *trans_a = trans_t;
    if(typeid(KpointType) == typeid(std::complex<double>)) trans_a = trans_c;

    KpointType alpha(vel);
    KpointType beta(0.0);
    if(Hmat != NULL)
    {
        // Apply operators on each wavefunction
        RmgTimer *RT1 = new RmgTimer("2-TDDFT: apply operators");
        KpointType *h_psi = (KpointType *)tmp_arrayT;
        ApplyBlockedHamiltonian(kptr, h_psi, vtot_eig, vxc_psi);
        rmg::sync_device();
        delete RT1;

        // Compute A matrix
        RT1 = new RmgTimer("2-TDDFT: Hmatrix setup/reduce");

        rmg::gemm(trans_a, trans_n, nstates, nstates, pbasis_noncoll, alpha, orbital_storage, pbasis_noncoll, tmp_arrayT, pbasis_noncoll, beta, Hmat, nstates);

        //MPI_Allreduce(MPI_IN_PLACE, (double *)Aij, nstates * nstates * factor, MPI_DOUBLE, MPI_SUM, grid_comm);
        rmg::block_allreduce((double *)Hmat, (size_t)nstates * (size_t)nstates * (size_t)factor , grid_comm);
        delete(RT1);
    }


    if(Smat != NULL)
    {
        RmgTimer *RT1 = new RmgTimer("2-TDDFT: Smatrix setup/reduce");
        // Compute S matrix
        rmg::gemm (trans_a, trans_n, nstates, nstates, pbasis_noncoll, alpha, prev_orbital_storage, pbasis_noncoll, orbital_storage, pbasis_noncoll, beta, Smat, nstates);

        //MPI_Allreduce(MPI_IN_PLACE, (double *)Sij, nstates * nstates * factor, MPI_DOUBLE, MPI_SUM, grid_comm);
        rmg::block_allreduce((double *)Smat, (size_t)nstates * (size_t)nstates * (size_t)factor , grid_comm);
        delete(RT1);
    }


}

