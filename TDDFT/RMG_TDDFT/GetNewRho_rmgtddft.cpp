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
#include "rmg_hvector.h"


#include "../Headers/common_prototypes.h"
#include "../Headers/common_prototypes1.h"
#include "prototypes_tddft.h"
#include "RmgParallelFft.h"
#include "rmg_gemm.h"
#include "blas_driver.h"
#include "Prolong.h"
#include "GatherScatter.h"

template void GetNewRho_rmgtddft<double, double, double>(Kpoint<double> *,spinobj<double> &rho, double *rho_matrix, int numst, int tddft_start_state, Scalapack &Sp);
template void GetNewRho_rmgtddft<double, float, double>(Kpoint<double> *,spinobj<double> &rho, double *rho_matrix, int numst, int tddft_start_state, Scalapack &Sp);
template void GetNewRho_rmgtddft<double, double, std::complex<double>>(Kpoint<double> *,spinobj<double> &rho, std::complex<double> *rho_matrix, int numst, int tddft_start_state, Scalapack &Sp);
template void GetNewRho_rmgtddft<double, float, std::complex<double>>(Kpoint<double> *,spinobj<double> &rho, std::complex<double> *rho_matrix, int numst, int tddft_start_state, Scalapack &Sp);
template void GetNewRho_rmgtddft<std::complex<double>, std::complex<double>, std::complex<double> >(Kpoint<std::complex<double>> *, spinobj<double> &rho, std::complex<double> *rho_matrix, int numst, int tddft_start_state, Scalapack &Sp);
template void GetNewRho_rmgtddft<std::complex<double>, std::complex<float>, std::complex<double> >(Kpoint<std::complex<double>> *, spinobj<double> &rho, std::complex<double> *rho_matrix, int numst, int tddft_start_state, Scalapack &Sp);
template <typename KpointType, typename CalType, typename MatrixType>
void GetNewRho_rmgtddft (Kpoint<KpointType> *kptr, spinobj<double> &rho_k, MatrixType *rho_matrix, int numst, int tddft_start_state, Scalapack &Sp)
{

    // rho_matrix are distributed either in tiledMM or scalpack 
    /* for parallel libraries */


    using TypeV = get_scalar_t<CalType>;                // Result: float
    int pbasis = get_P0_BASIS();
    int pbasis_noncoll = pbasis * ct.noncoll_factor;
    // in the noncollinear case, noncoll_factor = 2, and number of rho components = 4
    int n_rho = ct.noncoll_factor * ct.noncoll_factor;

    if(!ct.norm_conserving_pp) {
        rmg::error("\n tddft not programed for ultrasoft \n");
    }

    std::vector<double> occ_ground(numst);
    for (int istate = 0; istate < numst; istate++)
    {
        occ_ground[istate] = 
            kptr->Kstates[istate + tddft_start_state].occupation[0];
    }

#if CUDA_ENABLED || HIP_ENABLED 
    // xpsi is a device buffer in this case and GpuProductBr is a GPU functions to do
    // the reduction over numst.
    
    // rho_matrix is on device, with MatrixType
    CalType one = 1.0, zero = 0.0;
    TypeV *rho_temp_dev;

    CalType *rho_matrix_dev = (CalType *)ct.get_gmatrix_gpu(numst * numst * sizeof(CalType));
    double *occ_dev;


    rmg::hvector<TypeV> rho_temp(pbasis*n_rho);
    gpuMalloc((void **)&rho_temp_dev, pbasis*n_rho * sizeof(TypeV));
    gpuMalloc((void **)&occ_dev, numst * sizeof(double));
    gpuMemcpy(occ_dev, occ_ground.data(),  numst * sizeof(double), gpuMemcpyHostToDevice);

    //rho_matrix_dev[i,i] = rho_matrix[i,i] - occ_dev[i]

    int nprocs = pct.local_comm_npes;
    int myrank = pct.local_rank;
    if(!ct.tddft_tiledMM)
    {
        nprocs = 1;
        myrank = 0;
    }

    GpuRhomatrixConvert(rho_matrix_dev, rho_matrix, occ_dev, numst, myrank, nprocs);

    if(ct.tddft_tiledMM)
    {
#if USE_NCCL
        size_t sendcount = numst * numst/nprocs * sizeof(CalType)/sizeof(TypeV);
        if( typeid(TypeV) == typeid(float) )
        {
            rmg::error(ncclAllGather(&rho_matrix_dev[numst*numst/nprocs*myrank], rho_matrix_dev, sendcount, ncclFloat, ct.nccl_local_comm, 0));
        }
        //else if( typeid(TypeV) == typeid(double) )
        else
        {
            rmg::error(ncclAllGather(&rho_matrix_dev[numst*numst/nprocs*myrank], rho_matrix_dev, sendcount, ncclDouble, ct.nccl_local_comm, 0));
        }
#else
            rmg::error("set tddft_tiledMM=false in the input file, need use nccl for this option");
#endif
    }


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
            psi_dev, pbasis_noncoll, rho_matrix_dev, numst, zero, xpsi, pbasis_noncoll);

    GpuProductBr(psi_dev, xpsi, rho_temp_dev, numst, pbasis);
    gpuMemcpy(rho_temp.data(), rho_temp_dev,  n_rho * pbasis * sizeof(TypeV), gpuMemcpyDeviceToHost);
#else
    if(typeid(KpointType) != typeid(CalType))
    {
        rmg::error("\n float precision not for CPU \n");
    }
    KpointType one = 1.0, zero = 0.0;
    RmgTimer *RT = new RmgTimer("TDDFT: rho: gemm");

    rmg::hvector<double> rho_temp(pbasis*n_rho);
    rho_temp.set(0.0);
    KpointType *rho_matrix_glob  = (KpointType *)ct.get_gmatrix(numst*numst*sizeof(KpointType));


    if(ct.tddft_tiledMM == 1)
    {
        double *rho_R = (double *)rho_matrix_glob;
        std::complex<double> *rho_C = (std::complex<double> *)rho_matrix_glob;
        size_t recvcount = numst * numst/pct.local_comm_npes * sizeof(KpointType)/sizeof(double);
        if(typeid(KpointType) == typeid(double))
        {
            for(int i = 0; i < numst * numst/pct.local_comm_npes; i++)
            {
                rho_R[pct.local_rank * numst * numst/pct.local_comm_npes + i] = std::real(rho_matrix[i]);
            }
        }
        else
        {
            for(int i = 0; i < numst * numst/pct.local_comm_npes; i++)
            {
                rho_C[pct.local_rank * numst * numst/pct.local_comm_npes + i] = rho_matrix[i];
            }
        }
        MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, 
                rho_matrix_glob, recvcount, MPI_DOUBLE, pct.local_comm);       
    }
    else
    {
        std::vector<KpointType> rho_matrix_dist(Sp.GetDistMdim() * Sp.GetDistNdim());
        double *rho_R = (double *)rho_matrix_dist.data();
        std::complex<double> *rho_C = (std::complex<double> *)rho_matrix_dist.data();

        if(typeid(KpointType) == typeid(double))
        {
            for(int i = 0; i < Sp.GetDistMdim() * Sp.GetDistNdim(); i++)
            {
                rho_R[i] = std::real(rho_matrix[i]);
            }
        }
        else
        {
            for(int i = 0; i < Sp.GetDistMdim() * Sp.GetDistNdim(); i++)
            {
                rho_C[i] = rho_matrix[i];
            }
        }

        for(int i = 0; i < numst * numst; i++) rho_matrix_glob[i] = 0.0;
        Sp.GatherEigvectors(rho_matrix_glob, rho_matrix_dist.data());
    }

    for(int i = 0; i< numst; i++) rho_matrix_glob[i*numst+i] -= occ_ground[i];
    KpointType *psi = &kptr->orbital_storage[tddft_start_state * pbasis_noncoll];
    KpointType *xpsi = kptr->work_cpu;
    rmg::gemm ("N", "N", pbasis_noncoll, numst, numst, one, 
            psi, pbasis_noncoll, rho_matrix_glob, numst, zero, xpsi, pbasis_noncoll);

    delete RT;
    RT = new RmgTimer("TDDFT: rho: dot");
    if(!ct.noncoll)
    {
        for(int st1 = 0; st1 < numst; st1++)
        {
            for(int idx = 0; idx < pbasis; idx++)
            {
                rho_temp[idx] += std::real(psi[st1 * pbasis_noncoll + idx] * std::conj(xpsi[st1 * pbasis_noncoll + idx]));
            }
        }
    }
    else
    {
        for(int st1 = 0; st1 < numst; st1++)
        {

            KpointType *one_psi = &psi[st1 * pbasis_noncoll];
            KpointType *one_xpsi = &xpsi[st1 * pbasis_noncoll];
            for(int idx = 0; idx < pbasis; idx++)
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
    CopyAndConvert(n_rho * pbasis, rho_temp.data(), rho_k_double.data());
    for(int irho = 0; irho< n_rho; irho++)
    {
        P.prolong(rho_k.data()+irho*fpbasis, rho_k_double.data() + irho*pbasis, dimx, dimy, dimz, half_dimx, half_dimy, half_dimz);

    }

    delete RT1;

#if CUDA_ENABLED || HIP_ENABLED 
    gpuFree(rho_temp_dev);
    gpuFree(occ_dev);
#endif
}

