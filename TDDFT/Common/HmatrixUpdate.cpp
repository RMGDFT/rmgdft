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
#include "rmg_alloc.h"
#include "rmg_error.h"
#include "rmgthreads.h"
#include "RmgTimer.h"
#include "RmgThread.h"
#include "GlobalSums.h"
#include "Kpoint.h"
#include "RmgGemm.h"
#include "Subdiag.h"
#include "GpuAlloc.h"
#include "ErrorFuncs.h"
#include "blas.h"
#include "RmgParallelFft.h"

#include "../../RMG/Headers/prototypes.h"
#include "common_prototypes.h"
#include "common_prototypes1.h"
#include "transition.h"
#include "prototypes_tddft.h"

#if GPU_ENABLED
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cublas_v2.h>
#endif


template void HmatrixUpdate<double>(Kpoint<double> *, double *, double *);
template void HmatrixUpdate<std::complex<double> >(Kpoint<std::complex<double>> *, double *, std::complex<double> *);

template <typename KpointType>
void HmatrixUpdate (Kpoint<KpointType> *kptr, double *vtot_eig, KpointType *Aij)
{

    BaseGrid *G = kptr->G;
    Lattice *L = kptr->L;

    int num_states = kptr->nstates;
    int pbasis = kptr->pbasis;
    double vel = L->get_omega() / ((double)(G->get_NX_GRID(1) * G->get_NY_GRID(1) * G->get_NZ_GRID(1)));

    static KpointType *tmp_arrayT;
    static KpointType *global_matrix1;

    int factor = 1;
    if(!ct.is_gamma) factor = 2;

    char *trans_t = "t";
    char *trans_n = "n";
    char *trans_c = "c";
    char *trans_a;
    if(typeid(KpointType) == typeid(std::complex<double>)) {
         trans_a = trans_c;
    }
    else {
        trans_a = trans_t;
    }   


    // First time through allocate pinned memory for buffers
    if(!tmp_arrayT) {

        int retval1 = MPI_Alloc_mem(pbasis * ct.max_states * sizeof(KpointType) , MPI_INFO_NULL, &tmp_arrayT);
        int retval2 = MPI_Alloc_mem(ct.max_states * ct.max_states * sizeof(KpointType) , MPI_INFO_NULL, &global_matrix1);

        if((retval1 != MPI_SUCCESS) || (retval2 != MPI_SUCCESS) ) {
            rmg_error_handler (__FILE__, __LINE__, "Memory allocation failure in HmatrixUpdate");
        }

        #if GPU_ENABLED
            RmgCudaError(__FILE__, __LINE__, cudaHostRegister( tmp_arrayT, pbasis * ct.max_states * sizeof(KpointType), cudaHostRegisterPortable), "Error registering memory.\n");
            RmgCudaError(__FILE__, __LINE__, cudaHostRegister( global_matrix1, ct.max_states * ct.max_states * sizeof(KpointType), cudaHostRegisterPortable), "Error registering memory.\n");
        #endif

    }

    

    for (int st1 = 0; st1 < num_states; st1++)
        for(int idx = 0; idx <pbasis; idx++)
        {
            tmp_arrayT[st1 * pbasis + idx] = kptr->Kstates[st1].psi[idx] * vtot_eig[idx];
        } 

    /* tmp_arrayT:   V|psi> */

    // Compute A matrix
    KpointType alpha(vel);
    KpointType beta(0.0);
    RmgGemm(trans_a, trans_n, num_states, num_states, pbasis, alpha, kptr->orbital_storage, pbasis, tmp_arrayT, 
            pbasis, beta, global_matrix1, num_states);

    MPI_Allreduce(MPI_IN_PLACE, (double *)global_matrix1, num_states * num_states * factor, MPI_DOUBLE, MPI_SUM, pct.grid_comm);

    // Store reduced Aij back in Aij matrix
    for(int idx = 0;idx < num_states*num_states;idx++) Aij[idx] = global_matrix1[idx];


}

