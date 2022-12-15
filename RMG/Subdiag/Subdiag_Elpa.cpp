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
#include "blas.h"

#include "common_prototypes.h"
#include "common_prototypes1.h"
#include "transition.h"

#if CUDA_ENABLED
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cublas_v2.h>
#endif

#include "Elpa.h"


template char * Subdiag_Elpa<double> (Kpoint<double> *kptr, double *Aij, double *Bij, double *Sij, double *eigs, double *eigvectors);
template char * Subdiag_Elpa<std::complex<double> > (Kpoint<std::complex<double>> *kptr, std::complex<double> *Aij, std::complex<double> *Bij, std::complex<double> *Sij, double *eigs, std::complex<double> *eigvectors);

Elpa *MainElpa;

// eigvectors holds Bij on input and the eigenvectors of the matrix on output
template <typename KpointType>
char * Subdiag_Elpa (Kpoint<KpointType> *kptr, KpointType *Aij, KpointType *Bij, KpointType *Sij, double *eigs, KpointType *eigvectors)
{

#if !SCALAPACK_LIBS || !USE_ELPA
    rmg_printf("This version of RMG was not built with Scalapack/Elpa support. Redirecting to LAPACK.");
    return Subdiag_Lapack(kptr, Aij, Bij, Sij, eigs, eigvectors);
#else

// Fill in upper triangle of Aij (elpa is only routine that needs it
    int blocksize = 16;
#pragma omp parallel for
    for (int i = 0; i < kptr->nstates; i += blocksize) {
        for (int j = i; j < kptr->nstates; j += blocksize) {
            for (int row = i; row < i + blocksize && row < kptr->nstates; row++) {
                for (int col = j; col < j + blocksize && col < kptr->nstates; col++) {
                    Aij[col*kptr->nstates+row]=Aij[row*kptr->nstates+col];
                }
            }
        }
    }

    KpointType ZERO_t(0.0);
    KpointType ONE_t(1.0);
    int info;

    static char *trans_n = "n";
    static char *trans_t = "t";

//  folded spectrum with scalapack is experimental. Uncomment the second line if you want to try it.
    bool use_folded = ((ct.use_folded_spectrum && (ct.scf_steps > 6)) || (ct.use_folded_spectrum && (ct.runflag == RESTART)));
    //use_folded = false;

    RmgTimer *DiagTimer;
    static int call_count, folded_call_count;
    if(use_folded)
    {
        DiagTimer = new RmgTimer("4-Diagonalization: elpa folded");
        folded_call_count++;
        rmg_printf("\nDiagonalization using folded elpa for step=%d  count=%d\n\n",ct.scf_steps, folded_call_count);
    }
    else
    {
        DiagTimer = new RmgTimer("4-Diagonalization: elpa");
        call_count++;
        rmg_printf("\nDiagonalization using elpa for step=%d  count=%d\n\n",ct.scf_steps, call_count);
    }

    int ione=1;
    int num_states = kptr->nstates;
    int factor = 1;
    if(!ct.is_gamma) factor=2;

    // Create 1 scalapack instance per grid_comm. We use a static Scalapack here since initialization on large systems is expensive
    if(!MainElpa) {
        // Need some code here to decide how to set the number of scalapack groups but for now use just 1
        int scalapack_groups = 1;
        int last = !ct.use_folded_spectrum;
        MainElpa = new Elpa(scalapack_groups, pct.thisimg, ct.images_per_node, num_states,
                     ct.scalapack_block_factor, last, pct.grid_comm);
        MainElpa->Init();
    }

    bool participates = MainElpa->Participates();
    int scalapack_nprow = MainElpa->GetRows();
    int scalapack_npcol = MainElpa->GetCols();
    int scalapack_npes = scalapack_nprow * scalapack_npcol;
    int root_npes = MainElpa->GetRootNpes();

    // Allocate and clear distributed matrices */
    static KpointType *distAij;
    static KpointType *distBij;
    static KpointType *distSij;
    static KpointType *distCij;
    KpointType *Cij = new KpointType[num_states * num_states]();


    RmgTimer *RT1;

    int *desca;
    int dist_length;
    static int saved_dist_length;
    
    if (participates) {

        dist_length = MainElpa->GetDistMdim() * MainElpa->GetDistNdim();
        desca = MainElpa->GetDistDesca();

        if(dist_length != saved_dist_length)
        {
            if(distSij) {MPI_Free_mem(distSij);distSij = NULL;}
            if(distBij) {MPI_Free_mem(distBij);distBij = NULL;}
            if(distAij) {MPI_Free_mem(distAij);distAij = NULL;}
            if(distCij) {MPI_Free_mem(distCij);distCij = NULL;}
        }
        if(!distAij) {
            int retval1 = MPI_Alloc_mem(dist_length * sizeof(KpointType) , MPI_INFO_NULL, &distAij);
            int retval2 = MPI_Alloc_mem(dist_length * sizeof(KpointType) , MPI_INFO_NULL, &distBij);
            int retval3 = MPI_Alloc_mem(dist_length * sizeof(KpointType) , MPI_INFO_NULL, &distSij);
            int retval4 = MPI_Alloc_mem(dist_length * sizeof(KpointType) , MPI_INFO_NULL, &distCij);
            if((retval1 != MPI_SUCCESS) || (retval2 != MPI_SUCCESS) || (retval3 != MPI_SUCCESS) || (retval4 != MPI_SUCCESS)) {
                rmg_error_handler (__FILE__, __LINE__, "Memory allocation failure in Subdiag_Scalapack");
            }
            saved_dist_length = dist_length;
        }

        // Copy matrices to dist arrays
        RT1 = new RmgTimer("4-Diagonalization: distribute matrices.");
        MainElpa->CopySquareMatrixToDistArray(Aij, distAij, num_states, desca);
        MainElpa->CopySquareMatrixToDistArray(Sij, distSij, num_states, desca);
        MainElpa->CopySquareMatrixToDistArray(eigvectors, distBij, num_states, desca);
        delete(RT1);

        // Copy A into Bij to pass to eigensolver
        for(int ix=0;ix < dist_length;ix++) distBij[ix] = distAij[ix];

    }


    if(use_folded) {
#if 0
        // We have to gather distBij back to Bij and then broadcast it to all nodes in the root
        // Sij is still present on all nodes in original form
        //MainElpa->GatherMatrix(Bij, distBij);
        //MainElpa->BcastRoot(Bij, factor * num_states * num_states, MPI_DOUBLE);

        FoldedSpectrumScalapack<double> ((Kpoint<double> *)kptr, num_states, (double *)Bij, (double *)distBij, num_states, (double *)Sij, num_states, eigs, (double *)Aij, MainElpa, SUBDIAG_LAPACK, ct.scalapack_block_factor);

        for(int idx=0;idx< num_states * num_states;idx++)eigvectors[idx] = Bij[idx];
        // Broadcast results if required
        if(root_npes != scalapack_npes) { 
            RT1 = new RmgTimer("4-Diagonalization: MPI_Bcast");
            MainElpa->BcastRoot(eigvectors, factor * num_states * num_states, MPI_DOUBLE);
            MainElpa->BcastRoot(eigs, num_states, MPI_DOUBLE);
            delete(RT1);
        }
#endif
    }
    else {

        if(participates) {

            /****************** Find Matrix of Eigenvectors *****************************/
            /* Using lwork=-1, PDSYGVX should return minimum required size for the work array */
            RT1 = new RmgTimer("4-Diagonalization: elpa_general");
            MainElpa->generalized_eigenvectors(distAij, distSij, eigs, distCij, false, &info);
            delete(RT1);
            
            // Clear eigvectors
            RT1 = new RmgTimer("4-Diagonalization: distribute ev");
            for (int idx = 0; idx < num_states*num_states; idx++) {
                eigvectors[idx] = ZERO_t;
            }
            // Gather distributed results from distCij into eigvectors
            //MainElpa->GatherMatrix(eigvectors, distAij);
            MainElpa->CopyDistArrayToSquareMatrix(eigvectors, distCij, num_states, desca);
            MainElpa->Allreduce(MPI_IN_PLACE, eigvectors, factor *num_states*num_states, MPI_DOUBLE, MPI_SUM);
            delete(RT1);

        }

        // Finally send eigenvalues and vectors to everyone 
        RT1 = new RmgTimer("4-Diagonalization: MPI_Bcast");
        MainElpa->BcastRoot(eigvectors, factor * num_states * num_states, MPI_DOUBLE);
        MainElpa->BcastRoot(eigs, num_states, MPI_DOUBLE);
        delete(RT1);

    }

    delete [] Cij;

    delete DiagTimer;

//    if(use_folded) return trans_t;   // Currently using pdsyngst in lower level routine. If
//    switch to FOLDED_GSE must uncomment
    return trans_n;
#endif

}
