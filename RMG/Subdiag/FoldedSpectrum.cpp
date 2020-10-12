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
#include "RmgMatrix.h"
#include "GpuAlloc.h"
#include "Gpufuncs.h"
#include "ErrorFuncs.h"
#include "blas.h"

#include "common_prototypes.h"
#include "common_prototypes1.h"
#include "transition.h"
#if CUDA_ENABLED
    #include <cuda.h>
    #include <cuda_runtime_api.h>
    #include <cublas_v2.h>
#endif



// Communicator for PE's participating in folded spectrum operations
MPI_Comm fs_comm;

static int FS_NPES;
static int FS_RANK;

// I have not finished updating this to work with complex orbitals yet. Given that the folded spectrum method is only
// useful for large systems which are almost always run at gamma with real orbitals it's not a high priority but should
// be straightforward enough to finish.
template int FoldedSpectrum<double> (BaseGrid *, int, double *, int, double *, int, double *, double *, double *, double *, int, int *, int, int);

template <typename KpointType>
int FoldedSpectrum(BaseGrid *Grid, int n, KpointType *A, int lda, KpointType *B, int ldb, KpointType *Asave, KpointType *Bsave,
		double *eigs, double *work, int lwork, int *iwork, int liwork, int driver)
{

    RmgTimer RT0("4-Diagonalization: fs");
    RmgTimer *RT1;

    // For mpi routines. Transfer twice as much data for complex orbitals
    int factor = 2;
    if(ct.is_gamma) factor = 1;

    // Initialize some bookkeeping stuff
    if(!FS_NPES) {

        FS_NPES = Grid->get_PE_X() * Grid->get_PE_Y() * Grid->get_PE_Z();
        fs_comm = pct.grid_comm;
        MPI_Comm_rank(fs_comm, &FS_RANK);

    }

    // Array storage for folded spectrum diagonalization communications
    int *fs_eigstart = new int[FS_NPES];
    int *fs_eigstop = new int[FS_NPES];
    int *fs_eigcounts = new int[FS_NPES];

    // Set up partition indices and bookeeping arrays
    int eig_start, eig_stop, eig_step;
    int n_start, n_win;
    FoldedSpectrumSetup(n, FS_NPES, FS_RANK, eig_start, eig_stop, eig_step,
                        n_start, n_win, fs_eigstart, fs_eigstop, fs_eigcounts, 1);


#if CUDA_ENABLED
    DeviceSynchronize();
    RT1 = new RmgTimer("4-Diagonalization: fs: folded");
    double *Vdiag = (double *)GpuMallocManaged(n * sizeof(double));
    double *tarr = (double *)GpuMallocManaged(n * sizeof(double));
    cudaMemcpy(Bsave, B, n*n*sizeof(double), cudaMemcpyDefault);

    //  Transform problem to standard eigenvalue problem
    RmgTimer *RT2 = new RmgTimer("4-Diagonalization: fs: transform");
    int its=7;
    cudaMemcpy(Asave, A, n*n*sizeof(double), cudaMemcpyDefault);
    DeviceSynchronize();
    FoldedSpectrumGSE<double> (Asave, Bsave, A, n, eig_start, eig_stop, fs_eigcounts, fs_eigstart, its, driver, fs_comm);
    delete(RT2);

    // Zero out matrix of eigenvectors (V) and eigenvalues n. G is submatrix storage
    KpointType *V = (KpointType *)GpuMallocManaged(n * n * sizeof(KpointType));
    KpointType *G = (KpointType *)GpuMallocManaged(n_win * n_win * sizeof(KpointType));
    GpuFill((double *)V, n*n, 0.0);
    double *n_eigs = new double[n]();

    DeviceSynchronize();

    // Do the submatrix along the diagonal to get starting values for folded spectrum
    //--------------------------------------------------------------------
    RT2 = new RmgTimer("4-Diagonalization: fs: submatrix");
    cudaMemcpy2D ( G, n_win*sizeof(double), &A[n_start*n + n_start], n*sizeof(double), n_win*sizeof(double), n_win, cudaMemcpyDefault); 

#else
    double *Vdiag = new double[n];
    double *tarr = new double[n];
    memcpy(Bsave, B, (size_t)n*(size_t)n*sizeof(double));
    RT1 = new RmgTimer("4-Diagonalization: fs: folded");

    //  Transform problem to standard eigenvalue problem
    RmgTimer *RT2 = new RmgTimer("4-Diagonalization: fs: transform");
    int its=7;
    memcpy(Asave, A, (size_t)n*(size_t)n*sizeof(double));
    FoldedSpectrumGSE<double> (Asave, Bsave, A, n, eig_start, eig_stop, fs_eigcounts, fs_eigstart, its, driver, fs_comm);
    delete(RT2);

    // Zero out matrix of eigenvectors (V) and eigenvalues n. G is submatrix storage
    KpointType *V = new KpointType[n*n]();
    KpointType *G = new KpointType[n_win*n_win]();
    double *n_eigs = new double[n]();

    // Do the submatrix along the diagonal to get starting values for folded spectrum
    //--------------------------------------------------------------------
    RT2 = new RmgTimer("4-Diagonalization: fs: submatrix");
    for(int ix = 0;ix < n_win;ix++){
        for(int iy = 0;iy < n_win;iy++){
            G[ix*n_win + iy] = A[(n_start+ix)*n + n_start + iy];
        }
    }

#endif


#if CUDA_ENABLED
    //int device = -1;
    //cudaGetDevice(&device);
    //cudaMemPrefetchAsync ( A, n_win*n_win*sizeof(double), device, NULL);
    DeviceSynchronize();
#endif

#if CUDA_ENABLED
    if(ct.cuda_version >= 9020)
        DsyevjDriver(G, &eigs[n_start], work, lwork, n_win, n_win);
    else
#endif
        DsyevdDriver(G, &eigs[n_start], work, lwork, n_win, n_win);

    // Store the eigen vector from the submatrix
#if CUDA_ENABLED
    DeviceSynchronize();
    GpuFill(Vdiag, n, 1.0);
    GpuNegate(G, n_win, Vdiag, 1, n_win);
    int i1=0, ione = 1;
    for(int i = 0;i < n_win;i++) if((i + n_start) == eig_start) i1 = i;
    cublasStatus_t custat;
    custat = cublasDdgmm(ct.cublas_handle, CUBLAS_SIDE_RIGHT, n_win, eig_step, &G[i1*n_win], n_win, &Vdiag[i1], ione, &V[(i1 + n_start)*n + n_start], n);
    DeviceSynchronize();
    RmgCudaError(__FILE__, __LINE__, custat, "Problem executing cublasDdgmm.");
#else
    // Make sure same sign convention is followed by all eigenvectors
    for(int ix = 0;ix < n_win;ix++) {
        Vdiag[ix] = 1.0;
        if(G[ix*n_win + ix] < 0.0) Vdiag[ix] = -1.0;
    }
    for(int ix = 0;ix < n_win;ix++) {
        if(((n_start+ix) >= eig_start) && ((n_start+ix) < eig_stop)) 
        {
            for(int iy = 0;iy < n_win;iy++) {
                  V[(ix + n_start)*n + n_start + iy] = Vdiag[ix] * G[ix * n_win + iy];
            }
        }
    }
#endif
    delete(RT2);

    // And the eigenvalues
    for(int eig_index = eig_start;eig_index < eig_stop;eig_index++) {
        n_eigs[eig_index] = eigs[eig_index];
    }

#if HAVE_ASYNC_ALLREDUCE
    // Use async MPI to get a copy of the eigs to everyone. Will overlap with the iterator
    MPI_Request MPI_reqeigs;
    MPI_Iallreduce(MPI_IN_PLACE, n_eigs, n, MPI_DOUBLE, MPI_SUM, fs_comm, &MPI_reqeigs);
#else
    MPI_Allreduce(MPI_IN_PLACE, n_eigs, n, MPI_DOUBLE, MPI_SUM, fs_comm);
#endif


    // Apply folded spectrum to this PE's range of eigenvectors
#if CUDA_ENABLED
    DeviceSynchronize();
#endif
    RT2 = new RmgTimer("4-Diagonalization: fs: iteration");
    if(ct.folded_spectrum_iterations)
        FoldedSpectrumIterator(A, n, &eigs[eig_start], eig_stop - eig_start, &V[eig_start*n], -0.5, ct.folded_spectrum_iterations, driver);
#if CUDA_ENABLED
    DeviceSynchronize();
#endif
    delete(RT2);

#if HAVE_ASYNC_ALLREDUCE
    // Wait for eig request to finish and copy summed eigs from n_eigs back to eigs
    MPI_Wait(&MPI_reqeigs, MPI_STATUS_IGNORE);
#endif
    for(int idx = 0;idx < n;idx++) eigs[idx] = n_eigs[idx];

    // Make sure all PE's have all eigenvectors.
    RT2 = new RmgTimer("4-Diagonalization: fs: allreduce1");
    MPI_Allgatherv(MPI_IN_PLACE, eig_step * n * factor, MPI_DOUBLE, V, fs_eigcounts, fs_eigstart, MPI_DOUBLE, fs_comm);
    delete(RT2);


    // Gram-Schmidt ortho for eigenvectors.
    RT2 = new RmgTimer("4-Diagonalization: fs: Gram-Schmidt");

#if CUDA_ENABLED
    DeviceSynchronize();
#endif

    FoldedSpectrumOrtho(n, eig_start, eig_stop, fs_eigcounts, fs_eigstart, V, B, Asave, Bsave, driver, fs_comm);

#if CUDA_ENABLED
    DeviceSynchronize();
    cudaMemcpy(A, V, (size_t)n*(size_t)n*sizeof(double), cudaMemcpyDefault);
    //memcpy(A, V, (size_t)n*(size_t)n*sizeof(double));
    DeviceSynchronize();
#else
    memcpy(A, V, (size_t)n*(size_t)n*sizeof(double));
#endif
    delete(RT2);


    delete(RT1);
#if CUDA_ENABLED
    GpuFreeManaged(G);
    GpuFreeManaged(V);
    GpuFreeManaged(tarr);
    GpuFreeManaged(Vdiag);
#else
    delete [] G;
    delete [] V;
    delete [] tarr;
    delete [] Vdiag;
#endif
    delete [] n_eigs;

    delete [] fs_eigcounts;
    delete [] fs_eigstop;
    delete [] fs_eigstart;

    return 0;
} 



