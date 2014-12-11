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
#include "RmgException.h"
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


#define FOLDED_GSE 0



// Array storage for folded spectrum diagonalization communications
static int *fs_eigstart = NULL;
static int *fs_eigstop = NULL;
static int *fs_eigcounts = NULL;


// I have not finished updating this to work with complex orbitals yet. Given that the folded spectrum method is only
// useful for large systems which are almost always run at gamma with real orbitals it's not a high priority but should
// be straightforward enough to finish.
template int FoldedSpectrumScalapack<double> (Kpoint<double> *, int, double *, int, double *, int, double *, double *, Scalapack*, int);

// Just to note here the inputs rdistA,rdistB and rdistC are distributed matrices in the root communicator
template <typename KpointType>
int FoldedSpectrumScalapack(Kpoint<KpointType> *kptr, int n, KpointType *rA, int lda, KpointType *rB, int ldb, 
		double *eigs, KpointType *rC, Scalapack* MainSp, int driver)
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
    double rone = 1.0, scale = 1.0;;

    // For mpi routines. Transfer twice as much data for complex orbitals
    int factor = 2;
    if(ct.is_gamma) factor = 1;


    int NPES = Grid->get_PE_X() * Grid->get_PE_Y() * Grid->get_PE_Z();
    int root_pes = MainSp->GetNpes();
    if(root_pes < 8) {
        // Not much point in trying to use folded spectrum in this case
        throw RmgFatalException() << "Insufficient PE's to use folded spectrum in " << __FILE__ << " at line " << __LINE__ << ". Terminating.\n";
    }
    int FS_NPES = 16;               // Need to find a better way of setting this at some point.
    if(FS_NPES > root_pes) FS_NPES = root_pes;

    // Allocate some memory for our communication book keeping arrays
    if(!fs_eigstart) {
        fs_eigstart = new int[FS_NPES];
        fs_eigstop = new int[FS_NPES];
        fs_eigcounts = new int[FS_NPES];
    }

    // Set up partition indices and bookeeping arrays
    int eig_start, eig_stop, eig_step;
    int n_start, n_win;
    FoldedSpectrumSetup(n, FS_NPES, pct.gridpe, &eig_start, &eig_stop, &eig_step,
                        &n_start, &n_win, fs_eigstart, fs_eigstop, fs_eigcounts);

    // Folded spectrum scalapack instances use the Main scalapack instance as their root communicator.
    Scalapack *FSp = MainSp->GetNextScalapack();

    bool participates = FSp->Participates();
    int THIS_PE = FSp->GetGroupIndex();

    // Get dist_length and desca for the submatrices
    int dist_length = FSp->GetDistMdim() * FSp->GetDistNdim();
    int *p_desca = FSp->GetDistDesca();

    // Get dist length and desca for the full matrices
    int *f_desca = new int[FS_NPES*DLEN];
    if(participates) FSp->ComputeDesca(n, n, f_desca); // DistDesca is the base desca for a matrix dimenioned (n_win,n_win) but we also need a desca for the full (n,n) matrix

    int f_distmdim = FSp->ComputeMdim(n);
    int f_distndim = FSp->ComputeNdim(n);
    int f_dist_length = f_distmdim * f_distndim;

    MPI_Comm scalapack_comm = FSp->GetComm();
    int scalapack_nprow = FSp->GetRows();
    int scalapack_npcol = FSp->GetCols();
    if(dist_length == 0) dist_length = 1;   // Just to keep allocations from complaining

    static double *distAsave;
    static double *distBsave;
    static double *distA;
    static double *distB;
    static double *distV;
    static double *Asave;
    static double *Bsave;
    if(!distAsave) {
         int retval1 = MPI_Alloc_mem(f_dist_length * sizeof(double) , MPI_INFO_NULL, &distAsave);
         int retval2 = MPI_Alloc_mem(f_dist_length * sizeof(double) , MPI_INFO_NULL, &distBsave);
         int retval3 = MPI_Alloc_mem(f_dist_length * sizeof(double) , MPI_INFO_NULL, &distA);
         int retval4 = MPI_Alloc_mem(f_dist_length * sizeof(double) , MPI_INFO_NULL, &distB);
         int retval5 = MPI_Alloc_mem(f_dist_length * sizeof(double) , MPI_INFO_NULL, &distV);
         int retval6 = MPI_Alloc_mem(n * n * sizeof(double) , MPI_INFO_NULL, &Asave);
         int retval7 = MPI_Alloc_mem(n * n * sizeof(double) , MPI_INFO_NULL, &Bsave);
         if((retval1 != MPI_SUCCESS) || (retval2 != MPI_SUCCESS) || (retval3 != MPI_SUCCESS) || (retval4 != MPI_SUCCESS) || 
            (retval5 != MPI_SUCCESS) || (retval6 != MPI_SUCCESS) || (retval7 != MPI_SUCCESS)) {
            rmg_error_handler (__FILE__, __LINE__, "Memory allocation failure in FoldedSpectrum_Scalapack");
         }
    }

    // Gather distributed results back into compact array located on MainSp node with group_index=0;
    MainSp->GatherMatrix(Asave, rA);
    MainSp->GatherMatrix(Bsave, rB);

    // broadcast data to rest of nodes in main communicator that participate
    MPI_Bcast(Asave, n*n, MPI_DOUBLE, 0, MainSp->GetComm()); 
    MPI_Bcast(Bsave, n*n, MPI_DOUBLE, 0, MainSp->GetComm()); 

    // Copy A and B into distA, distB and distAsave and distBsave
    FSp->CopySquareMatrixToDistArray(Asave, distA, n, p_desca);
    FSp->CopySquareMatrixToDistArray(Bsave, distB, n, p_desca);
    for(int i = 0;i < f_dist_length;i++) distBsave[i] = distB[i];

    double *Vdiag = new double[n];
    double *tarr = new double[n];
    KpointType *V = new KpointType[n*n]();
    KpointType *G = new KpointType[n*n]();
    //for(int ix = 0;ix < n*n;ix++) Bsave[ix] = B[ix];

    RmgTimer *RT2;
    if(participates) {

#if !FOLDED_GSE
        RT1 = new RmgTimer("Diagonalization: fs: cholesky");
        //  Form a Cholesky factorization of B
        pdpotrf_ (cuplo, &n, distB, &ione, &ione, p_desca, &info);
        if( info != 0 )
            rmg_error_handler(__FILE__, __LINE__, "pdpotrf failure");
        delete(RT1);
#endif


        RT1 = new RmgTimer("Diagonalization: fs: folded");
        KpointType *NULLptr = NULL;

        //  Transform problem to standard eigenvalue form
        RT2 = new RmgTimer("Diagonalization: fs: transform");

#if !FOLDED_GSE
        {
            // Get pdsyngst_ workspace
            int lwork = -1;
            double lwork_tmp;
            pdsyngst_(&itype, cuplo, &n, distA, &ione, &ione, p_desca, distB, &ione, &ione, p_desca, &scale, &lwork_tmp, &lwork, &info);
            lwork = 2*(int)lwork_tmp;
            double *work = new double[lwork];

            pdsyngst_(&itype, cuplo, &n, distA, &ione, &ione, p_desca, distB, &ione, &ione, p_desca, &scale, work, &lwork, &info);
            delete [] work;
        }
        if( info != 0 )
            rmg_error_handler(__FILE__, __LINE__, "pdsyngst failure");
#else

        int its=7;
        double *T = new double[n*n];

        for(int idx = 0;idx < n*n;idx++) Asave[idx] = A[idx];
        FoldedSpectrumGSE<double> (Asave, Bsave, T, n, eig_start, eig_stop, fs_eigcounts, fs_eigstart, its, driver);

        // Copy back to A
        for(int ix=0;ix < n*n;ix++) A[ix] = T[ix];
        delete [] T;

#endif
        delete(RT2);

        // Zero out matrix of eigenvectors (V)
        for(int ix = 0;ix < f_dist_length;ix++) distV[ix] = ZERO_t;

        // AX=lambdaX  store a copy of distA in distAsave
        for(int ix = 0;ix < f_dist_length;ix++) distAsave[ix] = distA[ix];
        //QMD_dcopy (n*n, a, 1, Asave, 1);

     
        // Do the submatrix along the diagonal to get starting values for folded spectrum
        //--------------------------------------------------------------------
        RT2 = new RmgTimer("Diagonalization: fs: submatrix");
        {
            int lwork=-1, liwork=-1, liwork_tmp;
            double lwork_tmp;
            lwork = -1;
            int f_n_start = 1;
            pdsyevd_(jobz, cuplo, &n_win, distA, &f_n_start, &f_n_start, p_desca, &eigs[n_start],
                    distV, &f_n_start, &f_n_start, p_desca, &lwork_tmp, &lwork, &liwork_tmp, &liwork, &info);
            lwork = 2*(int)lwork_tmp;
            liwork = liwork_tmp;
            double *work = new double[lwork];
            int *iwork = new int[liwork];

            f_n_start = n_start + 1;
            pdsyevd_(jobz, cuplo, &n_win, distA, &f_n_start, &f_n_start, p_desca, &eigs[n_start],
                     distV, &f_n_start, &f_n_start, p_desca, work, &lwork, iwork, &liwork, &info);

            delete [] iwork;
            delete [] iwork;
        }
        if( info != 0 ) 
                rmg_error_handler(__FILE__, __LINE__, "pdsyevd failure");

        //--------------------------------------------------------------------

        // distV contains the eigenvectors of the submatrix. Copy back to full array
        FSp->CopyDistArrayToSquareMatrix(G, distV, n, f_desca);
        for(int ix = n_start;ix < n_start + n_win;ix++) {
            Vdiag[ix] = ONE_t;
            if(G[ix*n + ix] < 0.0) Vdiag[ix-n_start] = -ONE_t;
        }

        // Store the eigen vector from the submatrix
        for(int ix = 0;ix < n_win;ix++) {

            if(((n_start+ix) >= eig_start) && ((n_start+ix) < eig_stop)) {

                for(int iy = 0;iy < n_win;iy++) {
                      V[(ix + n_start)*n + n_start + iy] = Vdiag[ix] * G[(ix + n_start)*n + n_start + iy];
                }

            }

        }
        // Copy back to distributed form
        FSp->CopySquareMatrixToDistArray(V, distV, n, f_desca);
        delete(RT2);


        // Apply folded spectrum to this PE's range of eigenvectors
        RT2 = new RmgTimer("Diagonalization: fs: iteration");
        FoldedSpectrumIterator(Asave, n, &eigs[eig_start], eig_stop - eig_start, &V[eig_start*n], -0.5, 10, driver);
        delete(RT2);

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
        //for(int idx = 0;idx < n*n;idx++) A[idx] = V[idx];
        delete(RT2);


#if !FOLDED_GSE
        RT2 = new RmgTimer("Diagonalization: fs: dtrsm");
        //dtrsm (side, cuplo, trans_t, diag, &n, &n, &rone, B, &ldb, A, &lda);
        delete(RT2);
#endif

        delete(RT1);
    } // end if participates

    delete [] G;
    delete [] V;

    delete [] tarr;
    delete [] Vdiag;
    delete [] f_desca;

    return 0;
} 
