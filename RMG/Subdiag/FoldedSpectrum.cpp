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


#define FOLDED_GSE 1



// Array storage for folded spectrum diagonalization communications
static int *fs_eigstart = NULL;
static int *fs_eigstop = NULL;
static int *fs_eigcounts = NULL;


// I have not finished updating this to work with complex orbitals yet. Given that the folded spectrum method is only
// useful for large systems which are almost always run at gamma with real orbitals it's not a high priority but should
// be straightforward enough to finish.
template int FoldedSpectrum<double> (Kpoint<double> *, int, double *, int, double *, int, double *, double *, int, int *, int, double *, int);

template <typename KpointType>
int FoldedSpectrum(Kpoint<KpointType> *kptr, int n, KpointType *A, int lda, KpointType *B, int ldb, 
		double *eigs, double *work, int lwork, int *iwork, int liwork, KpointType *C, int driver)
{

    RmgTimer RT0("Diagonalization: fs:");
    RmgTimer *RT1;

    KpointType ZERO_t(0.0);
    KpointType ONE_t(1.0);

    BaseGrid *Grid = kptr->G;
    Lattice *L = kptr->L;

    char *trans_t="t", *trans_n="n";
    char *cuplo = "l", *side="l", *diag="n", *jobz="V";

    int ione=1, itype=1, info=0;
    double rone = 1.0;

    // For mpi routines. Transfer twice as much data for complex orbitals
    int factor = 2;
    if(ct.is_gamma) factor = 1;


    int NPES = Grid->get_PE_X() * Grid->get_PE_Y() * Grid->get_PE_Z();


    // Allocate some memory for our communication book keeping arrays
    if(!fs_eigstart) {
        fs_eigstart = new int[NPES];
        fs_eigstop = new int[NPES];
        fs_eigcounts = new int[NPES];
    }

    // Set up partition indices and bookeeping arrays
    int eig_start, eig_stop, eig_step;
    int n_start, n_win;
//if(ct.scf_steps > 15)ct.folded_spectrum_width=0.25;
    FoldedSpectrumSetup(n, NPES, pct.gridpe, &eig_start, &eig_stop, &eig_step,
                        &n_start, &n_win, fs_eigstart, fs_eigstop, fs_eigcounts);


    double *Vdiag = new double[n];
    double *tarr = new double[n];
#if GPU_ENABLED
    double *Asave = (double *)GpuMallocHost(n * n * sizeof(double));
    double *Bsave = (double *)GpuMallocHost(n * n * sizeof(double));
#else
    double *Asave = new double[n*n];
    double *Bsave = new double[n*n];
#endif
    for(int ix = 0;ix < n*n;ix++) Bsave[ix] = B[ix];


#if !FOLDED_GSE
    RT1 = new RmgTimer("Diagonalization: fs: cholesky");
    //  Form a Cholesky factorization of B
    dpotrf(cuplo, &n, B, &ldb, &info);
    if( info != 0 )
        rmg_error_handler(__FILE__, __LINE__, "dpotrf failure");
    delete(RT1);
#endif


    RT1 = new RmgTimer("Diagonalization: fs: folded");
    KpointType *NULLptr = NULL;

    //  Transform problem to standard eigenvalue problem
    RmgTimer *RT2 = new RmgTimer("Diagonalization: fs: tridiagonal");

#if !FOLDED_GSE
    dsygst(&itype, cuplo, &n, A, &lda, B, &ldb, &info);
    if( info != 0 )
        rmg_error_handler(__FILE__, __LINE__, "dsygst failure");
#else

    int its=7;
#if GPU_ENABLED
    double *T = (double *)GpuMallocHost(n * n * sizeof(double));
#else
    double *T = new double[n*n];
#endif
    for(int idx = 0;idx < n*n;idx++) Asave[idx] = A[idx];
    FoldedSpectrumGSE<double> (Asave, Bsave, T, n, eig_start, eig_stop, fs_eigcounts, fs_eigstart, its, driver);

    // Copy back to A
    for(int ix=0;ix < n*n;ix++) A[ix] = T[ix];
#if GPU_ENABLED
    GpuFreeHost(T);
#else
    delete [] T;
#endif

#endif
    delete(RT2);
    // Zero out matrix of eigenvectors (V) and eigenvalues n. G is submatrix storage
#if GPU_ENABLED
    KpointType *V = (KpointType *)GpuMallocHost(n * n * sizeof(KpointType));
    KpointType *G = (KpointType *)GpuMallocHost(n * n * sizeof(KpointType));
    for(int ix = 0;ix < n * n;ix++) V[ix] = ZERO_t;
    for(int ix = 0;ix < n * n;ix++) G[ix] = ZERO_t;
#else
    KpointType *V = new KpointType[n*n]();
    KpointType *G = new KpointType[n*n]();
#endif
    double *n_eigs = new double[n]();

    // AX=lambdaX  store a copy of A in Asave
    for(int idx = 0;idx < n*n;idx++) Asave[idx] = A[idx];
    //QMD_dcopy (n*n, a, 1, Asave, 1);

 
    // Do the submatrix along the diagonal to get starting values for folded spectrum
    //--------------------------------------------------------------------
    RT2 = new RmgTimer("Diagonalization: fs: submatrix");
    for(int ix = 0;ix < n_win;ix++){
        for(int iy = 0;iy < n_win;iy++){
            G[ix*n_win + iy] = Asave[(n_start+ix)*n + n_start + iy];
        }
    }

    for(int idx = 0;idx < n_win * n_win;idx++) A[idx] = G[idx];
    //QMD_dcopy (n_win * n_win, G, 1, a, 1);

#if (GPU_ENABLED && MAGMA_LIBS)
    magma_dsyevd(MagmaVec, MagmaLower, n_win, A, n_win, &eigs[n_start],
                    work, lwork,
                    iwork, liwork,
                    &info);
#else
    dsyevd(jobz, cuplo, &n_win, A, &n_win, &eigs[n_start],
                    work, &lwork,
                    iwork, &liwork,
                    &info);
#endif
    if( info != 0 ) 
            rmg_error_handler(__FILE__, __LINE__, "dsyevd failure");

    //--------------------------------------------------------------------

    for(int idx = 0;idx < n_win * n_win;idx++) G[idx] = A[idx];
    //QMD_dcopy (n_win * n_win, a, 1, G, 1);

    for(int ix = 0;ix < n_win;ix++) {
        Vdiag[ix] = ONE_t;
        if(G[ix*n_win + ix] < 0.0) Vdiag[ix] = -ONE_t;
    }

    // Store the eigen vector from the submatrix
    for(int ix = 0;ix < n_win;ix++) {

        if(((n_start+ix) >= eig_start) && ((n_start+ix) < eig_stop)) {

            for(int iy = 0;iy < n_win;iy++) {
                  V[(ix + n_start)*n + n_start + iy] = Vdiag[ix] * G[ix * n_win + iy];
            }

        }

    }
    delete(RT2);

    // And the eigenvalues
    for(int eig_index = eig_start;eig_index < eig_stop;eig_index++) {
        n_eigs[eig_index] = eigs[eig_index];
    }

    // Use async MPI to get a copy of the eigs to everyone. Will overlap with the iterator
    MPI_Request MPI_reqeigs;
    MPI_Iallreduce(MPI_IN_PLACE, n_eigs, n, MPI_DOUBLE, MPI_SUM, pct.grid_comm, &MPI_reqeigs);


    // Apply folded spectrum to this PE's range of eigenvectors
    RT2 = new RmgTimer("Diagonalization: fs: iteration");
    FoldedSpectrumIterator(Asave, n, &eigs[eig_start], eig_stop - eig_start, &V[eig_start*n], -0.5, 10, driver);
    delete(RT2);

    // Wait for eig request to finish and copy summed eigs from n_eigs back to eigs
    MPI_Wait(&MPI_reqeigs, MPI_STATUS_IGNORE);
    for(int idx = 0;idx < n;idx++) eigs[idx] = n_eigs[idx];

    // Make sure all PE's have all eigenvectors.
    RT2 = new RmgTimer("Diagonalization: fs: allreduce1");
    MPI_Allgatherv(MPI_IN_PLACE, eig_step * n * factor, MPI_DOUBLE, V, fs_eigcounts, fs_eigstart, MPI_DOUBLE, pct.grid_comm);
    delete(RT2);


    // Gram-Schmidt ortho for eigenvectors.
    RT2 = new RmgTimer("Diagonalization: fs: Gram-Schmidt");

#if !FOLDED_GSE
    FoldedSpectrumOrtho(n, eig_start, eig_stop, fs_eigcounts, fs_eigstart, V, NULLptr, driver);
#else
    FoldedSpectrumOrtho(n, eig_start, eig_stop, fs_eigcounts, fs_eigstart, V, B, driver);
#endif
    for(int idx = 0;idx < n*n;idx++) A[idx] = V[idx];
    delete(RT2);


#if !FOLDED_GSE
    RT2 = new RmgTimer("Diagonalization: fs: dtrsm");
    dtrsm (side, cuplo, trans_t, diag, &n, &n, &rone, B, &ldb, A, &lda);
    delete(RT2);
#endif

    delete(RT1);
#if GPU_ENABLED
    GpuFreeHost(G);
    GpuFreeHost(V);
    GpuFreeHost(Bsave);
    GpuFreeHost(Asave);
#else
    delete [] G;
    delete [] V;
    delete [] Bsave;
    delete [] Asave;
#endif
    delete [] n_eigs;

    delete [] tarr;
    delete [] Vdiag;

    return 0;
} 
