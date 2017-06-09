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
#include <omp.h>
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
#include "ErrorFuncs.h"
#include "Scalapack.h"
#include "blas.h"

#include "prototypes.h"
#include "common_prototypes.h"
#include "common_prototypes1.h"
#include "transition.h"

#if GPU_ENABLED
    #include <cuda.h>
    #include <cuda_runtime_api.h>
    #include <cublas_v2.h>
    #if MAGMA_LIBS
        #include <magma.h>
    #endif
#endif

// Used to convert the generalized eigenvalue problem A*x = lambda*B*x into standard form
// B^(-1)*A*x = lambda*x
//
// Inputs are A and B matrices. Output is Z = B^(-1)*A. Parallelized over nodes, istart
// and istop delineate the rows of the matrix Z that this nodes is responsible for.
// fs_eigcounts and fs_eigstart are used to collect the data in the MPI_Allgatherv call
// at the end. A and B are both overwritten.
//
template void FoldedSpectrumScalapackGSE<double> (double *, double *, double *, int, int, int, int *, int *, int, Scalapack* MainSp);
template <typename DataType>
void FoldedSpectrumScalapackGSE(DataType *A, DataType *B, DataType *Z, int n, int istart, int istop, int *fs_eigcounts, int *fs_eigstart, int iterations, Scalapack *MainSp)
{
    RmgTimer RT0("Diagonalization: fs: GSE");
    RmgTimer *RT1 = new RmgTimer("Diagonalization: fs: GSE-setup");
    DataType ZERO_t(0.0);
    DataType ONE_t(1.0);

    DataType *NULLptr = NULL;
    int istep = istop - istart;
    int ione = 1;
    double rone=1.0, rzero=0.0;

    // For mpi routines. Transfer twice as much data for complex orbitals
    int factor = 2;
    if(ct.is_gamma) factor = 1;

    char *trans_n = "n";

    Scalapack *FSp = MainSp->GetNextScalapack();
    // A PE is actually a scalapack instance

    // descriptor for full matrix in main scalapack instance
    //int m_f_dist_length = MainSp->ComputeMdim(n) *  MainSp->ComputeNdim(n);

    // Get dist length and desca for the submatrices
    int s_s_desca[DLEN];
    int s_s_dist_length = FSp->ComputeDesca(istep, n, s_s_desca);
    int s_f_desca[DLEN];
    //int s_f_dist_length = FSp->ComputeDesca(n, n, s_f_desca);


    DataType *D = new DataType[n];

    DataType *T1 = new DataType[n*n]();
    DataType *Zdist = new DataType[s_s_dist_length];
    DataType *T1dist = new DataType[s_s_dist_length];

    // Set up D^(-1)
    for(int ix = 0;ix < n;ix++) D[ix] = 1.0 / B[ix*n + ix];

    // Initial starting guess is just the identity
    pdlaset_( "A", &istep, &n, &rone, &rzero, (double *)Zdist, &ione, &ione, s_s_desca );

    // (I - D-1 * B)
    for(int st1 = 0;st1 < n;st1++){
        for(int st2 = 0;st2 < n;st2++){
            T1[st1*n + st2] -= D[st2] * B[st1*n + st2];
        }
    }
    FSp->CopySquareMatrixToDistArray(T1, T1dist, n, s_f_desca);

    delete(RT1);


    RT1 = new RmgTimer("Diagonalization: fs: GSE-Second term");
    // Compute D^(-1) * B * X and store in B
    for(int st1 = istart;st1 < istop;st1++){
        for(int st2 = 0;st2 < n;st2++){
            B[st1*n + st2] = D[st2] * A[st1*n + st2];
        }
    }
//    FSp->CopySquareMatrixToDistArray(B, Bdist, n, s_f_desca);

    delete(RT1);

    RT1 = new RmgTimer("Diagonalization: fs: GSE-First term");

    // outer loop over steps
    for(int step = 0;step < iterations;step++) {

        // Compute (I - D-1 * B) * Z(step) and store in A
        RmgGemm(trans_n, trans_n, n, istep, n, ONE_t, T1, n, &Z[istart*n], n, ZERO_t, &A[istart*n], n, NULLptr, NULLptr, NULLptr, false, false, false, true);

        // Finally generate Z(step+1) = (I - D-1 * B) * Z(step) + D^(-1) * B * X 
        //for(int ix=0;ix < n*n;ix++) Z[ix] = A[ix] + B[ix];
        for(int st1 = istart;st1 < istop;st1++){
            for(int st2 = 0;st2 < n;st2++){
                Z[st1*n + st2] =  A[st1*n + st2] +  B[st1*n + st2];
            }
        }

    }    

    delete(RT1);
    delete [] T1;
    delete [] D;


    // Make sure everybody has a copy
    RT1 = new RmgTimer("Diagonalization: fs: GSE-Allgatherv");
    MPI_Allgatherv(MPI_IN_PLACE, istep * n * factor, MPI_DOUBLE, Z, fs_eigcounts, fs_eigstart, MPI_DOUBLE, pct.grid_comm);
    delete(RT1);

}
