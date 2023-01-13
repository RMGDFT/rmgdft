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
#include "RmgTimer.h"
#include "GlobalSums.h"
#include "GpuAlloc.h"
#include "Kpoint.h"
#include "Subdiag.h"
#include "RmgGemm.h"
#include "blas.h"
#include "RmgMatrix.h"
#include "Functional.h"


#include "common_prototypes.h"
#include "common_prototypes1.h"
#include "transition.h"


#include "Scalapack.h"

Scalapack *MainSp;

template <typename KpointType>
char * Subdiag_Scalapack (Kpoint<KpointType> *kptr, KpointType *hpsi);


template char * Subdiag_Scalapack<double> (Kpoint<double> *kptr, double *hpsi);
template char * Subdiag_Scalapack<std::complex<double> > (Kpoint<std::complex<double>> *kptr, std::complex<double> *hpsi);

extern Scalapack *MainSp;

// eigvectors holds Bij on input and the eigenvectors of the matrix on output
template <typename KpointType>
char * Subdiag_Scalapack (Kpoint<KpointType> *kptr, KpointType *hpsi)
{
    RmgTimer *RT1;
    int my_rank;
    MPI_Comm_rank(pct.grid_comm, &my_rank);

    // Create 1 scalapack instance per grid_comm. We use a static Scalapack here
    // since initialization on large systems is expensive
    if(!MainSp)
    {
        // Need some code here to decide how to set the number of scalapack groups but for now use just 1
        int last = !ct.use_folded_spectrum;
        MainSp = new Scalapack(ct.subdiag_groups, pct.thisimg, ct.images_per_node, kptr->nstates,
                     ct.scalapack_block_factor, last, pct.grid_comm);
    }

    int *desca = MainSp->GetDistDesca();
    bool participates = MainSp->Participates();

    static KpointType *distAij;
    static KpointType *distBij;
    static KpointType *distSij;

    int dist_length=0;
    static int saved_dist_length;
    
    static KpointType *global_matrix1;
    int nstates = kptr->nstates;

    if (participates)
    {
        dist_length = MainSp->GetDistMdim() * MainSp->GetDistNdim();

        if(dist_length != saved_dist_length)
        {
            if(distSij) {MPI_Free_mem(distSij);distSij = NULL;}
            if(distBij) {MPI_Free_mem(distBij);distBij = NULL;}
            if(distAij) {MPI_Free_mem(distAij);distAij = NULL;}
        }
        if(!distAij) {
            int retval1 = MPI_Alloc_mem(dist_length * sizeof(KpointType) , MPI_INFO_NULL, &distAij);
            int retval2 = MPI_Alloc_mem(dist_length * sizeof(KpointType) , MPI_INFO_NULL, &distBij);
            int retval3 = MPI_Alloc_mem(dist_length * sizeof(KpointType) , MPI_INFO_NULL, &distSij);
            if((retval1 != MPI_SUCCESS) || (retval2 != MPI_SUCCESS) || (retval3 != MPI_SUCCESS)) {
                rmg_error_handler (__FILE__, __LINE__, "Memory allocation failure in Subdiag_Scalapack");
            }
            saved_dist_length = dist_length;
        }
    }

    // For Matrix multiplications
    double vel = kptr->L->get_omega() /
                ((double)(kptr->G->get_NX_GRID(1) * kptr->G->get_NY_GRID(1) * kptr->G->get_NZ_GRID(1)));
    KpointType alpha(1.0);
    KpointType alphavel(vel);
    KpointType beta(0.0);

    // For MPI routines
    int factor = 1;
    if(!ct.is_gamma) factor = 2;

    int pbasis_noncoll = kptr->pbasis * ct.noncoll_factor;
    // State array is 2 * the number of states in length but memory above
    // the first set of nstates is unused in this routine so we can use it
    // as temporary space.
    KpointType *tmp_arrayT = kptr->Kstates[0].psi;
    tmp_arrayT += nstates * pbasis_noncoll ;

    char *trans_t = "t";
    char *trans_n = "n";
    char *trans_c = "c";
    char *trans_a = trans_t;
    if(typeid(KpointType) == typeid(std::complex<double>)) trans_a = trans_c;


    float *Tij;
#if HIP_ENABLED || CUDA_ENABLED
    if(!global_matrix1) gpuMallocHost((void **)&global_matrix1, nstates * nstates * sizeof(KpointType));
    Tij = (float *)GpuMallocHost((nstates+2)*nstates*factor/2 * sizeof(float));
    double *eigs;
    eigs = (double *)GpuMallocHost(2*nstates * sizeof(double));
#else
    if(!global_matrix1) global_matrix1 = new KpointType[nstates * nstates];
    Tij = new float[((nstates + 2)*nstates*factor)/2]();
    double *eigs = new double[2*nstates];
#endif

    KpointType *D = new KpointType[nstates];

//  For CPU only case and CUDA with managed memory psi_d is the same as orbital_storage but
//  for HIP its a GPU buffer.
    KpointType *psi_d = kptr->orbital_storage;
#if HIP_ENABLED
    // For HIP which does not yet have managed memory copy wavefunctions into array on GPU
    // and use it repeatedly to compute the matrix elements. This is much faster but puts
    // more pressure on GPU memory. A blas implementation that overlapped communication and
    // computation would make this unnecessary.
    gpuMalloc((void **)&psi_d, nstates * pbasis_noncoll * sizeof(KpointType));
    gpuMemcpy(psi_d, kptr->orbital_storage, nstates * pbasis_noncoll * sizeof(KpointType), gpuMemcpyHostToDevice);
#endif
#if CUDA_ENABLED
    if(ct.gpu_managed_memory == false && ct.use_cublasxt == false)
    {
        gpuMalloc((void **)&psi_d, nstates * pbasis_noncoll * sizeof(KpointType));
        gpuMemcpy(psi_d, kptr->orbital_storage, nstates * pbasis_noncoll * sizeof(KpointType), gpuMemcpyHostToDevice);
    }
#endif

    RT1 = new RmgTimer("4-Diagonalization: matrix");
    RmgTimer *RT1a = new RmgTimer("4-Diagonalization: matrix: Gemm");
    if(ct.is_gamma)
        RmgSyrkx("L", "T", nstates, pbasis_noncoll, alphavel, psi_d, pbasis_noncoll, hpsi, pbasis_noncoll, beta, global_matrix1, nstates);
    else
        RmgGemm(trans_a, trans_n, nstates, nstates, pbasis_noncoll, alphavel, psi_d, pbasis_noncoll, hpsi, pbasis_noncoll, beta, global_matrix1, nstates);
    delete RT1a;

    // Hij is symmetric or Hermetian so pack into triangular array for reduction call. Use Tij for scratch space
    RT1a = new RmgTimer("4-Diagonalization: matrix: SqToTr");
    if(typeid(KpointType) == typeid(std::complex<double>))
        PackSqToTr("L", nstates, (std::complex<double> *)global_matrix1, nstates, (std::complex<float> *)Tij);
    else
        PackSqToTr("L", nstates, (double *)global_matrix1, nstates, (float *)Tij);
    delete RT1a;

    // Save diagonal elements
    RT1a = new RmgTimer("4-Diagonalization: matrix: Allreduce");
    for(int i=0;i < nstates;i++) D[i] = global_matrix1[i*nstates + i];
    BlockAllreduce((float *)Tij, (size_t)(nstates+2)*(size_t)nstates * (size_t)factor / 2, pct.grid_comm);
    delete RT1a;


    RT1a = new RmgTimer("4-Diagonalization: matrix: TrToSq");
    if(typeid(KpointType) == typeid(std::complex<double>))
    {
        UnPackSqToTr("L", nstates, (std::complex<double> *)global_matrix1, nstates, (std::complex<float> *)Tij);
    }
    else
    {
        UnPackSqToTr("L", nstates, (double *)global_matrix1, nstates, (float *)Tij);
    }
    delete RT1a;

    // Reduce diagonal elements in double precision and pack back into global_matrix1
    RT1a = new RmgTimer("4-Diagonalization: matrix: Allreduce");
    MPI_Allreduce(MPI_IN_PLACE, (double *)D, nstates * factor, MPI_DOUBLE, MPI_SUM, pct.grid_comm);
    delete RT1a;

    for(int i=0;i < nstates;i++) global_matrix1[i*nstates + i] = D[i];

#if USE_ELPA
    // Fill in upper triangle of Aij (elpa is only routine that needs it)
    Scalapack::FillUpper(global_matrix1, nstates);
#endif

    // Distribute Aij
    RmgTimer *RT3 = new RmgTimer("4-Diagonalization: matrix: distribute");
    MainSp->CopySquareMatrixToDistArray(global_matrix1, distAij, nstates, desca);
    delete RT3;

    // Compute S matrix
    RT1a = new RmgTimer("4-Diagonalization: matrix: Gemm");
    if(ct.norm_conserving_pp && ct.is_gamma)
    {
        RmgSyrkx("L", "T", nstates, pbasis_noncoll, alphavel, psi_d, pbasis_noncoll,  psi_d, pbasis_noncoll, beta, global_matrix1, nstates);
    }
    else
    {
        RmgGemm (trans_a, trans_n, nstates, nstates, pbasis_noncoll, alphavel, psi_d, pbasis_noncoll, kptr->ns, pbasis_noncoll, beta, global_matrix1, nstates);
    }
    delete RT1a;

    // Save diagonal elements in double precision
    for(int i=0;i < nstates;i++) D[i] = global_matrix1[i*nstates + i];

    // Sij is symmetric or Hermetian so pack into triangular array for reduction call.
    RT1a = new RmgTimer("4-Diagonalization: matrix: SqToTr");
    if(typeid(KpointType) == typeid(std::complex<double>))
        PackSqToTr("L", nstates, (std::complex<double> *)global_matrix1, nstates, (std::complex<float> *)Tij);
    else
        PackSqToTr("L", nstates, (double *)global_matrix1, nstates, (float *)Tij);
    delete RT1a;

    RT1a = new RmgTimer("4-Diagonalization: matrix: Allreduce");
    BlockAllreduce((float *)Tij, (size_t)(nstates+2)*(size_t)nstates * (size_t)factor / 2, pct.grid_comm);
    delete RT1a;

    RT1a = new RmgTimer("4-Diagonalization: matrix: TrToSq");
    if(typeid(KpointType) == typeid(std::complex<double>))
    {
        UnPackSqToTr("L", nstates, (std::complex<double> *)global_matrix1, nstates, (std::complex<float> *)Tij);
    }
    else
    {
        UnPackSqToTr("L", nstates, (double *)global_matrix1, nstates, (float *)Tij);
    }
    delete RT1a;

    // Reduce diagonal elements in double precision and pack back into global_matrix1
    RT1a = new RmgTimer("4-Diagonalization: matrix: Allreduce");
    MPI_Allreduce(MPI_IN_PLACE, (double *)D, nstates * factor, MPI_DOUBLE, MPI_SUM, pct.grid_comm);
    delete RT1a;
    for(int i=0;i < nstates;i++) global_matrix1[i*nstates + i] = D[i];

    // Fill in upper triangle of S
    Scalapack::FillUpper(global_matrix1, nstates);
    delete RT1;

    // Distribute Sij
    RT3 = new RmgTimer("4-Diagonalization: matrix: distribute");
    MainSp->CopySquareMatrixToDistArray(global_matrix1, distSij, nstates, desca);
    delete RT3;

    static int call_count;
    if(ct.subdiag_driver == SUBDIAG_ELPA)
    {
        rmg_printf("\nDiagonalization using elpa for step=%d  count=%d\n\n",ct.scf_steps,call_count);
    }
    else
    {
        rmg_printf("\nDiagonalization using scalapack for step=%d  count=%d\n\n",ct.scf_steps,call_count);
    }
    call_count++;

    if (participates)
    {
        RmgTimer *RT2 = new RmgTimer("4-Diagonalization: ELPA/PDSYGVX/PZHEGVX");

        // Copy Aij into Bij to pass to eigensolver
        for(int ix=0;ix < dist_length;ix++) distBij[ix] = distAij[ix];

        MainSp->generalized_eigenvectors(distAij, distSij, eigs, distBij);
        delete RT2;
        
        RT2 = new RmgTimer("4-Diagonalization: EigenvectorGathering");
        // Gather distributed results from distAij into eigvectors
        MainSp->GatherEigvectors(global_matrix1, distAij);
        delete RT2;
    }

    // Finally send eigenvalues and vectors to everyone 
    RT1 = new RmgTimer("4-Diagonalization: MPI_Bcast");
    MainSp->BcastRoot(global_matrix1, factor * nstates * nstates, MPI_DOUBLE);
    MainSp->BcastRoot(eigs, nstates, MPI_DOUBLE);
    delete RT1;

    // Begin rotation
    RT1 = new RmgTimer("4-Diagonalization: Update orbitals");
    // If subspace diagonalization is used every step, use eigenvalues obtained here 
    // as the correct eigenvalues
    if (ct.diag == 1) {
        for (int st1 = 0; st1 < nstates; st1++) {
            kptr->Kstates[st1].eig[0] = eigs[st1];
        }
    }

    RmgGemm(trans_n, trans_n, pbasis_noncoll, nstates, nstates, alpha, 
            psi_d, pbasis_noncoll, global_matrix1, nstates, beta, tmp_arrayT, pbasis_noncoll);

    // And finally copy them back
    size_t istart = 0;
    size_t tlen = (size_t)nstates * (size_t)pbasis_noncoll * sizeof(KpointType); 
    if(Verify ("freeze_occupied", true, kptr->ControlMap))
    {
        for(int istate = 0;istate < nstates;istate++)
        {
            if(kptr->Kstates[istate].occupation[0] > 1.0e-10) kptr->highest_occupied = istate;
        }
        istart = (size_t)(kptr->highest_occupied + 1)*(size_t)pbasis_noncoll;
        tlen = (size_t)nstates * (size_t)pbasis_noncoll - (size_t)(kptr->highest_occupied + 1) * (size_t)pbasis_noncoll;
    }

    // And finally make sure they follow the same sign convention when using hybrid XC
    // Optimize this for GPUs!
    if(ct.xc_is_hybrid)
    {
        for(int istate=0;istate < nstates;istate++)
        {
            if(std::real(global_matrix1[istate*nstates + istate]) < 0.0)
            {
                for(int idx=0;idx < pbasis_noncoll;idx++) kptr->Kstates[istate].psi[idx] = -kptr->Kstates[istate].psi[idx];
            }
        }
    }

    memcpy(&kptr->orbital_storage[istart], &tmp_arrayT[istart], tlen);

    // Rotate EXX
    if(ct.xc_is_hybrid && Functional::is_exx_active())
    {
        tlen = nstates * pbasis_noncoll * sizeof(KpointType);
        // vexx is not in managed memory yet so that might create an issue
        RmgGemm(trans_n, trans_n, pbasis_noncoll, nstates, nstates, alpha, 
                kptr->vexx, pbasis_noncoll, global_matrix1, nstates, beta, tmp_arrayT, pbasis_noncoll);
        memcpy(kptr->vexx, tmp_arrayT, tlen);
    }
    delete RT1;
// End rotation

#if HIP_ENABLED || CUDA_ENABLED
#if CUDA_ENABLED
    if(ct.gpu_managed_memory == false && ct.use_cublasxt == false)
    {
        gpuFree(psi_d);
    }
#endif
#if HIP_ENABLED
    gpuFree(psi_d);
#endif
    GpuFreeHost(eigs);
    GpuFreeHost(Tij);
#else
    delete [] D;
    delete [] eigs;
    delete [] Tij;
#endif


#if CUDA_ENABLED || HIP_ENABLED
    // After the first step this matrix does not need to be as large
    if(ct.scf_steps == 0) {gpuFreeHost(global_matrix1);global_matrix1 = NULL;}
#endif

//    if(use_folded) return trans_t;   // Currently using pdsyngst in lower level routine. If
//    switch to FOLDED_GSE must uncomment
    return trans_n;

}
