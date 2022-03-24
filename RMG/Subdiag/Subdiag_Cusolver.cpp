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
#include "Kpoint.h"
#include "Subdiag.h"
#include "RmgGemm.h"
#include "GpuAlloc.h"
#include "Gpufuncs.h"
#include "ErrorFuncs.h"
#include "blas.h"

#include "common_prototypes.h"
#include "common_prototypes1.h"
#include "transition.h"
#include "RmgMatrix.h"

#if CUDA_ENABLED
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cublas_v2.h>
#include <thrust/fill.h>
#include <thrust/device_vector.h>
#endif


template char * Subdiag_Cusolver<double> (Kpoint<double> *kptr, double *Aij, double *Bij, double *Sij, double *eigs, double *eigvectors);
template char * Subdiag_Cusolver<std::complex<double> > (Kpoint<std::complex<double>> *kptr, std::complex<double> *Aij, std::complex<double> *Bij, std::complex<double> *Sij, double *eigs, std::complex<double> *eigvectors);

// eigvectors holds Bij on input and the eigenvectors of the matrix on output
    template <typename KpointType>
char * Subdiag_Cusolver (Kpoint<KpointType> *kptr, KpointType *Aij, KpointType *Bij, KpointType *Sij, double *eigs, KpointType *eigvectors)
{

    static char *trans_t = "t";
    static char *trans_n = "n";
    static int call_count, folded_call_count;
    int num_states = kptr->nstates;
    bool use_folded = ((ct.use_folded_spectrum && (ct.scf_steps > 6)) || (ct.use_folded_spectrum && (ct.runflag == RESTART)));
    RmgTimer *DiagTimer;
#if !CUDA_ENABLED
        rmg_printf("This version of RMG was not built with GPU support so Cusolver cannot be used. Redirecting to LAPACK.");
        return Subdiag_Lapack(kptr, Aij, Bij, Sij, eigs, eigvectors);
#endif

#if CUDA_ENABLED

    if(ct.num_usable_gpu_devices == 1 || pct.local_rank == 0)
    {


#if SCALAPACK_LIBS
        // For folded spectrum start with scalapack if available since cusolver is slow on really large problems
        if(ct.use_folded_spectrum && (ct.scf_steps <= 6)  && (ct.runflag != RESTART) && (num_states > 10000))
            return Subdiag_Scalapack (kptr, Aij, Bij, Sij, eigs, eigvectors);
#endif

        if(use_folded)
        {
            DiagTimer = new RmgTimer("4-Diagonalization: Eigensolver: cusolver folded");
            folded_call_count++;
            rmg_printf("\nDiagonalization using folded cusolver for step=%d  count=%d\n\n",ct.scf_steps, folded_call_count); 
        }
        else
        {
            DiagTimer = new RmgTimer("4-Diagonalization: Eigensolver: cusolver");
            call_count++;
            rmg_printf("\nDiagonalization using cusolver for step=%d  count=%d\n\n",ct.scf_steps, call_count); 
        }

        // Magma is not parallel across MPI procs so only have the local master proc on a node perform
        // the diagonalization. Then broadcast the eigenvalues and vectors to the remaining local nodes.
        // If folded spectrum is selected we only want the local master to participate on each node as
        // long as there are at least 12 nodes.


        // Copy A into eigvectors
        //    memcpy(eigvectors, Aij, (size_t)num_states * (size_t)num_states * sizeof(KpointType));
        KpointType *eigvectors_gpu, *Sij_gpu;
        double *eigs_gpu;
        cudaMallocManaged((void **)&eigvectors_gpu, (size_t)num_states * (size_t)num_states * sizeof(KpointType));
        cudaMallocManaged((void **)&Sij_gpu, (size_t)num_states * (size_t)num_states * sizeof(KpointType));
        cudaMallocManaged((void **)&eigs_gpu, (size_t)num_states * sizeof(KpointType));
        cudaMemcpy(eigvectors_gpu, Aij, (size_t)num_states * (size_t)num_states * sizeof(KpointType), cudaMemcpyDefault);
        cudaMemcpy(Sij_gpu, Sij, (size_t)num_states * (size_t)num_states * sizeof(KpointType), cudaMemcpyDefault);



        int *ifail = new int[num_states];
        int liwork = 6 * num_states + 4;
        int *iwork = new int[2*liwork];

        if(ct.is_gamma) {

            if(use_folded) {

                RmgTimer RT1("4-Diagonalization: Eigensolver: cusolver: folded");
                int lwork = num_states * num_states / 3 + num_states;
                lwork = std::max(lwork, 128000);
                double *work, *Aij_gpu, *Bij_gpu;
                cudaMalloc((void **)&work, lwork * sizeof(KpointType));
                cudaMallocManaged((void **)&Aij_gpu, (size_t)num_states * (size_t)num_states * sizeof(double));
                cudaMallocManaged((void **)&Bij_gpu, (size_t)num_states * (size_t)num_states * sizeof(double));
                FoldedSpectrum<double> (kptr->G, num_states, (double *)eigvectors_gpu, num_states, (double *)Sij_gpu, num_states, (double *)Aij_gpu, (double *)Bij_gpu, eigs_gpu, work, lwork, iwork, liwork, SUBDIAG_CUSOLVER);
                cudaFree(Bij_gpu);
                cudaFree(Aij_gpu);
                cudaFree(work);

            }
            else {

                int lwork = 3 * num_states * num_states + 8 * num_states;
                lwork = std::max(lwork, 128000);
                double *work;
                cudaMalloc((void **)&work, lwork * sizeof(KpointType));
                if(ct.cuda_version >= 9020 )
                {
                    RmgTimer RT1("4-Diagonalization: Eigensolver: cusolver: Dsygvj");
                    DsygvjDriver((double *)eigvectors_gpu, (double *)Sij_gpu, eigs_gpu, work, lwork, num_states, num_states);
                }
                else
                {
                    RmgTimer RT1("4-Diagonalization: Eigensolver: cusolver: Dsygvd");
                    DsygvdDriver((double *)eigvectors_gpu, (double *)Sij_gpu, eigs_gpu, work, lwork, num_states, num_states);
                }
                cudaFree(work);

            }

        }
        else {

            RmgTimer RT1("4-Diagonalization: Eigensolver: cusolver: Zhegvd");
            int lwork = 3 * num_states * num_states + 8 * num_states;
            lwork = std::max(lwork, 128000);
            ZhegvdDriver((std::complex<double> *)eigvectors_gpu, (std::complex<double> *)Sij_gpu, eigs_gpu, NULL, lwork, num_states, num_states);

        }

        delete [] iwork;
        delete [] ifail;

        cudaMemcpy(eigvectors, eigvectors_gpu, (size_t)num_states * (size_t)num_states * sizeof(KpointType), cudaMemcpyDefault);
        cudaMemcpy(eigs, eigs_gpu, (size_t)num_states * sizeof(double), cudaMemcpyDefault);
        cudaFree(eigs_gpu);
        cudaFree(Sij_gpu);
        cudaFree(eigvectors_gpu);
        delete DiagTimer;
    }

    if(ct.num_usable_gpu_devices > 1)
    {
        // end if is_local_master

        // If only one proc on this host participated broadcast results to the rest
        //    if((pct.procs_per_host > 1) && !(use_folded && (nodes < 12))) 
        RmgTimer RT1("4-Diagonalization: Eigensolver: cusolver: Bcast");
        int factor = 2;
        if(ct.is_gamma) factor = 1;
        MPI_Bcast(eigvectors, factor * num_states*num_states, MPI_DOUBLE, 0, pct.local_comm);
        MPI_Bcast(eigs, num_states, MPI_DOUBLE, 0, pct.local_comm);

    } 

    if(use_folded) return trans_t;
    return trans_n;
#endif
}



