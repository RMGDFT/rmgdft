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

#if HIP_ENABLED
#include <rocsolver.h>
#endif


template char * Subdiag_Rocsolver<double> (Kpoint<double> *kptr, double *Aij, double *Bij, double *Sij, double *eigs, double *eigvectors);
template char * Subdiag_Rocsolver<std::complex<double> > (Kpoint<std::complex<double>> *kptr, std::complex<double> *Aij, std::complex<double> *Bij, std::complex<double> *Sij, double *eigs, std::complex<double> *eigvectors);

// eigvectors holds Bij on input and the eigenvectors of the matrix on output
template <typename KpointType>
char * Subdiag_Rocsolver (Kpoint<KpointType> *kptr, KpointType *Aij, KpointType *Bij, KpointType *Sij, double *eigs, KpointType *eigvectors)
{

#if !HIP_ENABLED
        rmg_printf("This version of RMG was not built with GPU support so Rocsolver cannot be used. Redirecting to LAPACK.");
        return Subdiag_Scalapack (kptr, Aij, Bij, Sij, eigs, eigvectors);
#endif

#if HIP_ENABLED
    static char *trans_t = "t";
    static char *trans_n = "n";
    static int call_count, folded_call_count;
    int num_states = kptr->nstates;
    bool use_folded = ((ct.use_folded_spectrum && (ct.scf_steps >= 6)) || (ct.use_folded_spectrum && (ct.runflag == RESTART)));
    RmgTimer *DiagTimer;

    // For folded spectrum start with scalapack if available since rocsolver is slow on large problems
    if(ct.use_folded_spectrum && (ct.scf_steps <= 0)  && (ct.runflag != RESTART) && (num_states > 10000))
        return Subdiag_Scalapack (kptr, Aij, Bij, Sij, eigs, eigvectors);
    if(ct.scf_steps < 2 && ct.runflag != RESTART)
        return Subdiag_Scalapack (kptr, Aij, Bij, Sij, eigs, eigvectors);

    if(pct.is_local_master || ct.num_usable_gpu_devices == 1)
    {

        gpuSetDevice(ct.hip_dev);
        if(use_folded)
        {
            DiagTimer = new RmgTimer("4-Diagonalization: Eigensolver: rocsolver folded");
            folded_call_count++;
            rmg_printf("\nDiagonalization using folded rocsolver for step=%d  count=%d\n\n",ct.scf_steps, folded_call_count); 
        }
        else
        {
            DiagTimer = new RmgTimer("4-Diagonalization: Eigensolver: rocsolver");
            call_count++;
            rmg_printf("\nDiagonalization using rocsolver for step=%d  count=%d\n\n",ct.scf_steps, call_count); 
        }

        // Copy A into eigvectors
        //    memcpy(eigvectors, Aij, (size_t)num_states * (size_t)num_states * sizeof(KpointType));
        KpointType *eigvectors_gpu, *Sij_gpu;
        double *eigs_gpu;
        gpuMalloc((void **)&eigvectors_gpu, (size_t)num_states * (size_t)num_states * sizeof(KpointType));
        gpuMalloc((void **)&Sij_gpu, (size_t)num_states * (size_t)num_states * sizeof(KpointType));
        gpuMalloc((void **)&eigs_gpu, (size_t)num_states * sizeof(KpointType));
        gpuMemcpy(eigvectors_gpu, Aij, (size_t)num_states * (size_t)num_states * sizeof(KpointType), gpuMemcpyDefault);
        gpuMemcpy(Sij_gpu, Sij, (size_t)num_states * (size_t)num_states * sizeof(KpointType), gpuMemcpyDefault);



        int *ifail = new int[num_states];
        int liwork = 6 * num_states + 4;
        int *iwork = new int[2*liwork];

        if(ct.is_gamma) {

            if(use_folded) {

                RmgTimer RT1("4-Diagonalization: Eigensolver: rocsolver: folded");
                int lwork = num_states * num_states / 3 + num_states;
                lwork = std::max(lwork, 128000);
                double *work, *Aij_gpu, *Bij_gpu;
                gpuMalloc((void **)&work, lwork * sizeof(KpointType));
                gpuMallocManaged((void **)&Aij_gpu, (size_t)num_states * (size_t)num_states * sizeof(double));
                gpuMallocManaged((void **)&Bij_gpu, (size_t)num_states * (size_t)num_states * sizeof(double));
                FoldedSpectrum<double> (kptr->G, num_states, (double *)eigvectors_gpu, num_states, (double *)Sij_gpu, num_states, (double *)Aij_gpu, (double *)Bij_gpu, eigs_gpu, work, lwork, iwork, liwork, SUBDIAG_ROCSOLVER);
                gpuFree(Bij_gpu);
                gpuFree(Aij_gpu);
                gpuFree(work);

            }
            else {
                int lwork = 3 * num_states * num_states + 8 * num_states;
                lwork = std::max(lwork, 128000);
                double *work=NULL;
                // rocsolver does not need the work array
                //gpuMalloc((void **)&work, lwork * sizeof(KpointType));
                RmgTimer RT1("4-Diagonalization: Eigensolver: rocsolver: Dsygvj");
                DsygvjDriver((double *)eigvectors_gpu, (double *)Sij_gpu, eigs_gpu, work, lwork, num_states, num_states);
                //gpuFree(work);

            }

        }
        else {

            RmgTimer RT1("4-Diagonalization: Eigensolver: rocsolver: Zhegvd");
            int lwork = 3 * num_states * num_states + 8 * num_states;
            lwork = std::max(lwork, 128000);
            ZhegvdDriver((std::complex<double> *)eigvectors_gpu, (std::complex<double> *)Sij_gpu, eigs_gpu, NULL, lwork, num_states, num_states);

        }

        delete [] iwork;
        delete [] ifail;

        gpuMemcpy(eigvectors, eigvectors_gpu, (size_t)num_states * (size_t)num_states * sizeof(KpointType), gpuMemcpyDefault);
        gpuMemcpy(eigs, eigs_gpu, (size_t)num_states * sizeof(double), gpuMemcpyDefault);
        gpuFree(eigs_gpu);
        gpuFree(Sij_gpu);
        gpuFree(eigvectors_gpu);
        delete DiagTimer;
    }

    if(ct.num_usable_gpu_devices > 1)
    {
        // If only one proc on this host participated broadcast results to the rest
        //    if((pct.procs_per_host > 1) && !(use_folded && (nodes < 12))) 
        RmgTimer RT1("4-Diagonalization: Eigensolver: rocsolver: Bcast");
        int factor = 2;
        if(ct.is_gamma) factor = 1;
        MPI_Bcast(eigvectors, factor * num_states*num_states, MPI_DOUBLE, 0, pct.local_comm);
        MPI_Bcast(eigs, num_states, MPI_DOUBLE, 0, pct.local_comm);

    } 

    if(use_folded) return trans_t;
    return trans_n;
#endif
}



