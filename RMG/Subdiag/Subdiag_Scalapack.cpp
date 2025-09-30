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
extern Scalapack *MainSp;

template <typename KpointType>
char * Subdiag_Scalapack (Kpoint<KpointType> *kptr, KpointType *hpsi, int first_state, int last_state, Scalapack &SP);


template char * Subdiag_Scalapack<double> (Kpoint<double> *kptr, double *hpsi, int first_state, int last_state, Scalapack &SP);
template char * Subdiag_Scalapack<std::complex<double> > (Kpoint<std::complex<double>> *kptr, std::complex<double> *hpsi, int first_state, int last_state, Scalapack &SP);

// eigvectors holds Bij on input and the eigenvectors of the matrix on output
template <typename KpointType>
char * Subdiag_Scalapack (Kpoint<KpointType> *kptr, KpointType *hpsi, int first_state, int num_states, Scalapack &SP)
{
    RmgTimer *RT1;
    int my_rank;
    MPI_Comm_rank(pct.grid_comm, &my_rank);

    if (SP.GetN() != num_states)
    {
        std::cout << "Scalalapack dimension not match" << std::endl;
        std::cout << SP.GetN() << "NOT equal " << num_states << std::endl; 
        rmg_error_handler (__FILE__, __LINE__, "scalapack clase not correct");
    }

    int *desca = SP.GetDistDesca();
    bool participates = SP.Participates();

    static KpointType *distAij;
    static KpointType *distBij;
    static KpointType *distSij;

    int dist_length=0;
    static int saved_dist_length;


    if (participates)
    {
        dist_length = SP.GetDistMdim() * SP.GetDistNdim();

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
    int pbasis_noncoll = kptr->pbasis * ct.noncoll_factor;

    // State array is 2 * the number of states in length but memory above
    // the first set of nstates is unused in this routine so we can use it
    // as temporary space.
    KpointType *tmp_arrayT = kptr->Kstates[0].psi;
    tmp_arrayT += kptr->nstates * pbasis_noncoll ;

    char *trans_t = "t";
    char *trans_n = "n";
    char *trans_c = "c";
    char *trans_a = trans_t;
    if(typeid(KpointType) == typeid(std::complex<double>)) trans_a = trans_c;


    KpointType *psi = kptr->orbital_storage + first_state * pbasis_noncoll;
    KpointType *psi_dev;
#if HIP_ENABLED || CUDA_ENABLED || SYCL_ENABLED
    if(ct.gpu_managed_memory == false)
    {
        gpuMalloc((void **)&psi_dev, num_states * pbasis_noncoll * sizeof(KpointType));
        gpuMemcpy(psi_dev, psi, num_states * pbasis_noncoll * sizeof(KpointType), gpuMemcpyHostToDevice);
    }
    else
    {
        psi_dev = psi;
    }
#else
    psi_dev = psi;
#endif

    HS_Scalapack (num_states, pbasis_noncoll, psi_dev, &hpsi[first_state * pbasis_noncoll], &kptr->ns[first_state * pbasis_noncoll], desca, distAij, distSij);


#if HIP_ENABLED || CUDA_ENABLED || SYCL_ENABLED
    double *eigs;
    eigs = (double *)GpuMallocHost(2*num_states * sizeof(double));
#else
    double *eigs = new double[2*num_states];
#endif

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


        if(0 && ct.norm_conserving_pp)
        {
            SP.symherm_eigenvectors(distAij, eigs, distBij);
            for(int ix=0;ix < dist_length;ix++) distAij[ix] = distBij[ix];
        }
        else
        {
            // Copy Aij into Bij to pass to eigensolver
            for(int ix=0;ix < dist_length;ix++) distBij[ix] = distAij[ix];
            SP.generalized_eigenvectors(distAij, distSij, eigs, distBij);
        }
        delete RT2;

    }

    // Finally send eigenvalues and vectors to everyone 
    RT1 = new RmgTimer("4-Diagonalization: MPI_Bcast");
    SP.BcastRoot(eigs, num_states, MPI_DOUBLE);
    delete RT1;

    // If subspace diagonalization is used every step, use eigenvalues obtained here 
    // as the correct eigenvalues
    if (ct.diag == 1) {
        for (int st1 = 0; st1 < num_states; st1++) {
            kptr->Kstates[st1+first_state].eig[0] = eigs[st1];
        }
    }

    // Begin rotation
    RT1 = new RmgTimer("4-Diagonalization: Update orbitals");
    KpointType *matrix_diag = new KpointType[num_states];
    PsiUpdate(num_states, pbasis_noncoll, distAij, desca, psi_dev, &hpsi[first_state * pbasis_noncoll],  matrix_diag);

    // And finally copy them back
    size_t istart = (size_t)first_state * (size_t)pbasis_noncoll;
    size_t tlen = (size_t)num_states * (size_t)pbasis_noncoll * sizeof(KpointType); 

    memcpy(&kptr->orbital_storage[istart], &hpsi[first_state * pbasis_noncoll], tlen);

    // And finally make sure they follow the same sign convention when using hybrid XC
    // Optimize this for GPUs!
    if(ct.xc_is_hybrid && Functional::is_exx_active())
    {
        tlen = (size_t)num_states * (size_t)pbasis_noncoll * sizeof(KpointType);
#if HIP_ENABLED || CUDA_ENABLED || SYCL_ENABLED
        if(ct.gpu_managed_memory == false)
        {
            gpuMemcpy(psi_dev, &kptr->vexx[first_state * pbasis_noncoll], num_states * pbasis_noncoll * sizeof(KpointType), gpuMemcpyHostToDevice);
        }
        else
        {
            psi_dev = &kptr->vexx[first_state * pbasis_noncoll];
        }
#else
        psi_dev = &kptr->vexx[first_state * pbasis_noncoll];
#endif
        PsiUpdate(num_states, pbasis_noncoll, distAij, desca, psi_dev, &hpsi[first_state * pbasis_noncoll],  matrix_diag);
        memcpy(&kptr->vexx[first_state * pbasis_noncoll], &hpsi[first_state * pbasis_noncoll], tlen);
    }

    delete [] matrix_diag;

    delete RT1;
    // End rotation

#if HIP_ENABLED || CUDA_ENABLED || SYCL_ENABLED
    GpuFreeHost(eigs);
    if(ct.gpu_managed_memory == false)
    {
        gpuFree(psi_dev);
    }
#else
    delete [] eigs;
#endif


    //    if(use_folded) return trans_t;   // Currently using pdsyngst in lower level routine. If
    //    switch to FOLDED_GSE must uncomment
    return trans_n;

}



