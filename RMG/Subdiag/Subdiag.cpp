/*
 *
 * Copyright 2014 The RMG Project Developers. See the COPYRIGHT file 
 * at the top-level directory of this distribution or in the current
 * directory.
 * 
 * This file is part of RMG. 
 * RMG is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * any later version.
 *
 * RMG is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
*/

#include <complex>
#include "FiniteDiff.h"
#include "const.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "rmg_error.h"
#include "rmgthreads.h"
#include "RmgTimer.h"
#include "RmgThread.h"
#include "GlobalSums.h"
#include "Kpoint.h"
#include "RmgGemm.h"
#include "Gpufuncs.h"
#include "Subdiag.h"
#include "GpuAlloc.h"
#include "ErrorFuncs.h"
#include "blas.h"
#include "Solvers.h"
#include "Functional.h"
#include "RmgMatrix.h"

#include "common_prototypes.h"
#include "common_prototypes1.h"
#include "transition.h"

#if CUDA_ENABLED
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cublas_v2.h>
#endif



template void Kpoint<double>::Subdiag(double *, double *, int);
template void Kpoint<std::complex<double>>::Subdiag(double *, double *, int);

template <class KpointType> void Kpoint<KpointType>::Subdiag (double *vtot_eig, double *vxc_psi, int subdiag_driver)
{
    RmgTimer RT0("4-Diagonalization");

    bool potential_acceleration = (ct.potential_acceleration_constant_step > 0.0) && (ct.scf_steps > 0);
    double vel = L->get_omega() / ((double)(G->get_NX_GRID(1) * G->get_NY_GRID(1) * G->get_NZ_GRID(1)));


    // For MPI routines
    int factor = 1;
    if(!ct.is_gamma) factor = 2;

    int pbasis_noncoll = pbasis * ct.noncoll_factor;
    // State array is 4 * the number of states in length but memory above
    // the first set of nstates is unused in this routine so we can use it
    // as temporary space.
    KpointType *tmp_arrayT = Kstates[0].psi;
    tmp_arrayT += nstates * pbasis_noncoll ;

    static KpointType *global_matrix1;

    // We pad Bij since we use it as scratch space for the all reduce ops on Hij and Sij
#if HIP_ENABLED || CUDA_ENABLED
    if(!global_matrix1) gpuMallocHost((void **)&global_matrix1, nstates * nstates * sizeof(KpointType));     
    KpointType *Hij = (KpointType *)GpuMallocHost(nstates * nstates * sizeof(KpointType));
    KpointType *Bij = (KpointType *)GpuMallocHost(nstates * nstates * sizeof(KpointType));
    KpointType *Sij = (KpointType *)GpuMallocHost(nstates * nstates * sizeof(KpointType));
    double *eigs;
    gpuMallocHost((void **)&eigs, 2*nstates * sizeof(double));
#else
    if(!global_matrix1) global_matrix1 = new KpointType[nstates * nstates];
    KpointType *Hij = new KpointType[nstates * nstates];
    KpointType *Bij = new KpointType[nstates * nstates];
    KpointType *Sij = new KpointType[nstates * nstates];
    double *eigs = new double[2*nstates];
#endif

//  For CPU only case and CUDA with managed memory psi_d is the same as orbital_storage but
//  for HIP its a GPU buffer.
    KpointType *psi_d = orbital_storage;
#if HIP_ENABLED || CUDA_ENABLED
    // For HIP which does not yet have managed memory copy wavefunctions into array on GPU
    // and use it repeatedly to compute the matrix elements. This is much faster but puts
    // more pressure on GPU memory. A blas implementation that overlapped communication and
    // computation would make this unnecessary.
    gpuMalloc((void **)&psi_d, nstates * pbasis_noncoll * sizeof(KpointType));
    gpuMemcpy(psi_d, orbital_storage, nstates * pbasis_noncoll * sizeof(KpointType), gpuMemcpyHostToDevice);
#endif

    char *trans_t = "t";
    char *trans_n = "n";
    char *trans_c = "c";
    char *trans_a = trans_t;
    if(typeid(KpointType) == typeid(std::complex<double>)) trans_a = trans_c;

    // Apply operators on each wavefunction
    RmgTimer *RT1 = new RmgTimer("4-Diagonalization: apply operators");

    // Each thread applies the operator to one wavefunction
    KpointType *h_psi = (KpointType *)tmp_arrayT;

    ComputeHpsi(vtot_eig, vxc_psi, h_psi);

    delete(RT1);

    /* Operators applied and we now have
tmp_arrayT:  A|psi> + BV|psi> + B|beta>dnm<beta|psi> */

    DeviceSynchronize();

    // Compute A matrix
    RT1 = new RmgTimer("4-Diagonalization: matrix setup/reduce");
    KpointType alpha(1.0);
    KpointType alphavel(vel);
    KpointType beta(0.0);

    if(ct.is_gamma)
        RmgSyrkx("L", "T", nstates, pbasis_noncoll, alphavel, psi_d, pbasis_noncoll, tmp_arrayT, pbasis_noncoll, beta, Hij, nstates);
    else
        RmgGemm(trans_a, trans_n, nstates, nstates, pbasis_noncoll, alphavel, psi_d, pbasis_noncoll, tmp_arrayT, pbasis_noncoll, beta, Hij, nstates);

    // Hij is symmetric or Hermetian so pack into triangular array for reduction call. Use Bij for scratch space
    PackSqToTr("L", nstates, Hij, nstates, Bij);

#if HAVE_ASYNC_ALLREDUCE
    // Asynchronously reduce it
    MPI_Request MPI_reqHij;
    MPI_Request MPI_reqSij;
    if(ct.use_async_allreduce)
        MPI_Iallreduce(MPI_IN_PLACE, (double *)Bij, (nstates+2) * nstates * factor / 2, MPI_DOUBLE, MPI_SUM, grid_comm, &MPI_reqHij);
    else
        BlockAllreduce((double *)Bij, (size_t)(nstates+2)*(size_t)nstates * (size_t)factor / 2, grid_comm);
#else
    BlockAllreduce((double *)Bij, (size_t)(nstates+2)*(size_t)nstates * (size_t)factor / 2, grid_comm);
#endif

    // Compute S matrix
    if(ct.norm_conserving_pp && ct.is_gamma)
    {
        RmgSyrkx("L", "T", nstates, pbasis_noncoll, alphavel, psi_d, pbasis_noncoll,  psi_d, pbasis_noncoll, beta, Sij, nstates);
    }
    else
    {
        RmgGemm (trans_a, trans_n, nstates, nstates, pbasis_noncoll, alphavel, psi_d, pbasis_noncoll, ns, pbasis_noncoll, beta, Sij, nstates);
    }

    // Sij is symmetric or Hermetian so pack into triangular array for reduction call. Use global_matrix1 for scratch space
    PackSqToTr("L", nstates, Sij, nstates, global_matrix1);

#if HAVE_ASYNC_ALLREDUCE
    // Asynchronously reduce Sij request
    if(ct.use_async_allreduce)
        MPI_Iallreduce(MPI_IN_PLACE, (double *)global_matrix1, (nstates+2) * nstates * factor / 2, MPI_DOUBLE, MPI_SUM, grid_comm, &MPI_reqSij);
    else
        BlockAllreduce((double *)global_matrix1, (size_t)(nstates+2)*(size_t)nstates * (size_t)factor / 2, grid_comm);
#else
    BlockAllreduce((double *)global_matrix1, (size_t)(nstates+2)*(size_t)nstates * (size_t)factor / 2, grid_comm);
#endif


#if HAVE_ASYNC_ALLREDUCE
    // Wait for S request to finish and when done store copy in Sij
    if(ct.use_async_allreduce) MPI_Wait(&MPI_reqHij, MPI_STATUS_IGNORE);
    if(ct.use_async_allreduce) MPI_Wait(&MPI_reqSij, MPI_STATUS_IGNORE);
#endif

    UnPackSqToTr("L", nstates, Hij, nstates, Bij);
    UnPackSqToTr("L", nstates, Sij, nstates, global_matrix1);

    // Fill in upper triangle of S
    Scalapack::FillUpper(Sij, nstates);
    delete(RT1);

    // Dispatch to correct subroutine, eigs will hold eigenvalues on return and global_matrix1 will hold the eigenvectors.
    // The eigenvectors may be stored in row-major or column-major format depending on the type of diagonaliztion method
    // used. This is handled during the rotation of the orbitals by trans_b which is set by the driver routine.
    RT1 = new RmgTimer("4-Diagonalization: Eigensolver");
    char *trans_b = "n";
    switch(subdiag_driver) {

        case SUBDIAG_LAPACK:
            trans_b = Subdiag_Lapack (this, Hij, Bij, Sij, eigs, global_matrix1);
            break;
        // Redirect to scalapack for now
#if MAGMA_LIBS
        case SUBDIAG_MAGMA:
            trans_b = Subdiag_Magma (this, Hij, Bij, Sij, eigs, global_matrix1);
            break;
#endif
        case SUBDIAG_ELPA:
        case SUBDIAG_SCALAPACK:
            trans_b = Subdiag_Scalapack (this, Hij, Bij, Sij, eigs, global_matrix1);
            break;
        case SUBDIAG_CUSOLVER:
#if CUDA_ENABLED
            trans_b = Subdiag_Cusolver (this, Hij, Bij, Sij, eigs, global_matrix1);
#else
            trans_b = Subdiag_Lapack (this, Hij, Bij, Sij, eigs, global_matrix1);
#endif
            break;
        case SUBDIAG_ROCSOLVER:
#if HIP_ENABLED
            trans_b = Subdiag_Rocsolver (this, Hij, Bij, Sij, eigs, global_matrix1);
#else
            trans_b = Subdiag_Lapack (this, Hij, Bij, Sij, eigs, global_matrix1);
#endif
            break;
        default:
            rmg_error_handler(__FILE__, __LINE__, "Invalid subdiag_driver type");

    } // end switch
    delete(RT1);

    // If subspace diagonalization is used every step, use eigenvalues obtained here 
    // as the correct eigenvalues
    if (ct.diag == 1) {
        for (int st1 = 0; st1 < nstates; st1++) {
            Kstates[st1].eig[0] = eigs[st1];
        }
    }

    // Update the orbitals
    RT1 = new RmgTimer("4-Diagonalization: Update orbitals");

    RmgGemm(trans_n, trans_b, pbasis_noncoll, nstates, nstates, alpha, 
            psi_d, pbasis_noncoll, global_matrix1, nstates, beta, tmp_arrayT, pbasis_noncoll);

    // And finally copy them back
    size_t istart = 0;
    size_t tlen = (size_t)nstates * (size_t)pbasis_noncoll * sizeof(KpointType); 
    if(Verify ("freeze_occupied", true, ControlMap))
    {
        for(int istate = 0;istate < nstates;istate++)
        {
            if(Kstates[istate].occupation[0] > 1.0e-10) highest_occupied = istate;
        }
        istart = (size_t)(highest_occupied + 1)*(size_t)pbasis_noncoll;
        tlen = (size_t)nstates * (size_t)pbasis_noncoll - (size_t)(highest_occupied + 1) * (size_t)pbasis_noncoll;
    }

    // And finally make sure they follow the same sign convention when using hybrid XC
    // Optimize this for GPUs!
    if(ct.xc_is_hybrid)
    {
        for(int istate=0;istate < nstates;istate++)
        {
            if(std::real(global_matrix1[istate*nstates + istate]) < 0.0)
            {
                for(int idx=0;idx < pbasis_noncoll;idx++) Kstates[istate].psi[idx] = -Kstates[istate].psi[idx];
            }
        }
    }

    memcpy(&orbital_storage[istart], &tmp_arrayT[istart], tlen);

    // Rotate EXX
    if(ct.xc_is_hybrid && Functional::is_exx_active())
    {
        tlen = nstates * pbasis_noncoll * sizeof(KpointType);
        // vexx is not in managed memory yet so that might create an issue
        RmgGemm(trans_n, trans_b, pbasis_noncoll, nstates, nstates, alpha, 
                this->vexx, pbasis_noncoll, global_matrix1, nstates, beta, tmp_arrayT, pbasis_noncoll);
        memcpy(this->vexx, tmp_arrayT, tlen);
    }

    delete(RT1);

#if HIP_ENABLED || CUDA_ENABLED
    gpuFree(psi_d);
    gpuFreeHost(eigs);
    GpuFreeHost(Sij);
    GpuFreeHost(Bij);
    GpuFreeHost(Hij);
#else
    delete [] eigs;
    delete [] Sij;
    delete [] Bij;
    delete [] Hij;
#endif

#if CUDA_ENABLED || HIP_ENABLED
    // After the first step this matrix does not need to be as large
    if(ct.scf_steps == 0) {gpuFreeHost(global_matrix1);global_matrix1 = NULL;}
#endif

}

