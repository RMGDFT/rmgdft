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
#include <cmath>
#include <float.h>
#include "FiniteDiff.h"
#include "const.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "rmgthreads.h"
#include "RmgTimer.h"
#include "RmgThread.h"
#include "GlobalSums.h"
#include "Kpoint.h"
#include "RmgGemm.h"
#include "Mgrid.h"
#include "RmgException.h"
#include "Subdiag.h"
#include "Solvers.h"
#include "GpuAlloc.h"
#include "ErrorFuncs.h"
#include "RmgParallelFft.h"
#include "TradeImages.h"
#include "packfuncs.h"

#include "transition.h"
#include "blas.h"

#if GPU_ENABLED
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cublas_v2.h>
#if MAGMA_LIBS
    #include <magma.h>
#endif
#endif


// Kludge for now
extern Scalapack *MainSp;
extern void *GMatrix1;
extern void *GMatrix2;

// Solves for the eigenvalues and eigenvectors of a generalized Hermitian matrix.
//
//   A*v = lambda*B*v
//
//   INPUT:  A  =  N by N matrix
//   INPUT:  B  =  N by N matrix
//   INPUT:  V  =  N by N matrix
//   INPUT: ld  =  leading dimension of A,B,V
//   INPUT:  M  =  number of eigenvectors to solve for
//
// 

template int GeneralDiag<double>(double *A, double *B, double *eigs, double *V, int N, int M, int ld, int subdiag_driver);
template int GeneralDiag<std::complex<double>>(std::complex<double> *A, std::complex<double> *B, double *eigs, std::complex<double> *V, int N, int M, int ld, int subdiag_driver);

template <typename KpointType>
int GeneralDiag(KpointType *A, KpointType *B, double *eigs, KpointType *V, int N, int M, int ld, int subdiag_driver)
{

    int info = 0;
    switch(subdiag_driver) {
        case SUBDIAG_LAPACK:
            info = GeneralDiagLapack(A, B, eigs, V, N, M, ld);
            break;
        case SUBDIAG_SCALAPACK:
            if(!ct.is_gamma) {
                info = GeneralDiagLapack(A, B, eigs, V, N, M, ld);
            }

            else {
                info = GeneralDiagScaLapack((double *)A, (double *)B, eigs, (double *)V, N, M, ld);
            }
            break;
#if MAGMA_LIBS
        case SUBDIAG_MAGMA:
                info = GeneralDiagMagma(A, B, eigs, V, N, M, ld);
            break;
#endif
        default:
            throw RmgFatalException() << "Invalid diagonalization driver type in " << __FILE__ << " at line " << __LINE__ << "\n";

    }

    return info;
}

template int GeneralDiagLapack<double>(double *A, double *B, double *eigs, double *V, int N, int M, int ld);
template int GeneralDiagLapack<std::complex<double>>(std::complex<double> *A, std::complex<double> *B, double *eigs, std::complex<double> *V,
int N, int M, int ld);

template <typename KpointType>
int GeneralDiagLapack(KpointType *A, KpointType *B, double *eigs, KpointType *V, int N, int M, int ld)
{
    int info = 0;
    //bool use_folded = ((ct.use_folded_spectrum && (ct.scf_steps > 6)) || (ct.use_folded_spectrum && (ct.runflag == RESTART)));

    if(pct.is_local_master) {

        // Increase the resources available to this proc since the others on the local node
        // will be idle
        int nthreads = ct.THREADS_PER_NODE;
        if(pct.procs_per_host > 1) nthreads = pct.ncpus;
        omp_set_num_threads(nthreads);


        {
            int *ifail = new int[N];
            int lwork = 6 * N * N + 6 * N + 2;
            int liwork = 6*N;
            int eigs_found;
            double *work2 = new double[2*lwork];
            int *iwork = new int[liwork];
            double vx = 0.0;
            double tol = 1.0e-15;
            int itype = 1, ione = 1;

            KpointType *Asave = new KpointType[N*ld];
            KpointType *Bsave = new KpointType[N*ld];
            for(int i = 0;i < N*ld;i++) Asave[i] = A[i];
            for(int i = 0;i < N*ld;i++) Bsave[i] = B[i];

            if(N < M) throw RmgFatalException() << "M must be >= N in " << __FILE__ << " at line " << __LINE__ << "\n";

            if(ct.is_gamma)
            {
                if(M == N) {
                    dsygvd_(&itype, "V", "L", &N, (double *)A, &ld, (double *)B, &ld, eigs, work2, &lwork, iwork, &liwork, &info);
                    for(int ix=0;ix < N*ld;ix++) V[ix] = A[ix];
                }
                else if(N > M) {
                    dsygvx (&itype, "V", "I", "L", &N, (double *)A, &ld, (double *)B, &ld,
                            &vx, &vx, &ione, &M,  &tol, &eigs_found, eigs, (double *)V, &ld, work2,
                            &lwork, iwork, ifail, &info);
                }
            }
            else
            {
                if(M == N) {
                    zhegvd_(&itype, "V", "U", &N, (double *)A, &ld, (double *)B, &ld, eigs, work2, &lwork, &work2[lwork], &lwork, iwork, &liwork, &info);
                    for(int ix=0;ix < N*ld;ix++) V[ix] = A[ix];
                }
                else if(N > M) {
                    zhegvx (&itype, "V", "I", "U", &N, (double *)A, &ld, (double *)B, &ld,
                            &vx, &vx, &ione, &M,  &tol, &eigs_found, eigs, (double *)V, &ld, work2,
                            &lwork, &work2[lwork], iwork, ifail, &info);

                    //printf("\n %d infooo ", info);
                }

            }


            for(int i=0;i < N*ld;i++) A[i] = Asave[i];
            for(int i=0;i < N*ld;i++) B[i] = Bsave[i];
            delete [] Bsave;
            delete [] Asave;
            delete [] iwork;
            delete [] work2;
            delete [] ifail;

        }

        // Reset omp_num_threads
        omp_set_num_threads(ct.THREADS_PER_NODE);

    } // end if pct.is_local_master


    // If only one proc on this host participated broadcast results to the rest
    if(pct.procs_per_host > 1) {
        int factor = 2;
        if(ct.is_gamma) factor = 1;
        MPI_Bcast(V, factor * M*ld, MPI_DOUBLE, 0, pct.local_comm);
        MPI_Bcast(eigs, M, MPI_DOUBLE, 0, pct.local_comm);
        MPI_Bcast(&info, 1, MPI_INT, 0, pct.local_comm);
    }
    return info;
}


int GeneralDiagScaLapack(double *A, double *B, double *eigs, double *V, int N, int M, int ld)
{

#if !SCALAPACK_LIBS
    rmg_printf("This version of RMG was not built with Scalapack support. Redirecting to LAPACK.");
    return GeneralDiagLapack(A, B, eigs, V, N, M, ld);
#else

    double *global_matrix1 = (double *)GMatrix1;
    bool participates = MainSp->Participates();

    if (participates) {


        // Get dist length and desca
        int desca[DLEN];
        int dist_length = MainSp->ComputeDesca(N, N, desca);

        //desca = MainSp->GetDistDesca();
        MainSp->ComputeDesca(N, N, desca);

        // Allocate distributed matrices 
        double *distA = new double[dist_length * sizeof(double)];
        double *distB = new double[dist_length * sizeof(double)];
        double *distV = new double[dist_length * sizeof(double)];


        // Distribute matrices
        // This is a kludge for now but the Scalapack object has not yet been extended to 
        // handle the case when (ld != M) so we have to pack and repack the arrays locally
        // Also regular Subdiag_Scalapack must be called before this routine since we use
        // the Scalapack object setup there
        for(int i = 0;i < N;i++) {
            for(int j = 0;j < N;j++) {
                global_matrix1[j + i*N] = A[j + i*ld];
            }
        }
        MainSp->CopySquareMatrixToDistArray(global_matrix1, distA, N, desca);

        for(int i = 0;i < N;i++) {
            for(int j = 0;j < N;j++) {
                global_matrix1[j + i*N] = B[j + i*ld];
            }
        }
        MainSp->CopySquareMatrixToDistArray(global_matrix1, distB, N, desca);

        int ibtype = 1;
        double scale=1.0, rone = 1.0;
        int ione = 1;
        int info = 0;

        pdpotrf_("L", &N, (double *)distB,  &ione, &ione, desca,  &info);

        // Get pdsyngst_ workspace
        int lwork = -1;
        double lwork_tmp;
        pdsyngst_(&ibtype, "L", &N, (double *)distA, &ione, &ione, desca,
                (double *)distB, &ione, &ione, desca, &scale, &lwork_tmp, &lwork, &info);
        lwork = 2*(int)lwork_tmp;
        double *work2 = new double[lwork];

        pdsyngst_(&ibtype, "L", &N, (double *)distA, &ione, &ione, desca,
                (double *)distB, &ione, &ione, desca, &scale, work2, &lwork, &info);

        // Get workspace required
        lwork = -1;
        int liwork=-1;
        int liwork_tmp;
        pdsyevd_("V", "L", &N, (double *)distA, &ione, &ione, desca,
                eigs, (double *)distV, &ione, &ione, desca, &lwork_tmp, &lwork, &liwork_tmp, &liwork, &info);
        lwork = 16*(int)lwork_tmp;
        liwork = 16*N;
        double *nwork = new double[lwork];
        int *iwork = new int[liwork];

        // and now solve it 
        pdsyevd_("V", "L", &N, (double *)distA, &ione, &ione, desca,
                eigs, (double *)distV, &ione, &ione, desca, nwork, &lwork, iwork, &liwork, &info);

        pdtrsm_("Left", "L", "T", "N", &N, &N, &rone, (double *)distB, &ione, &ione, desca,
                (double *)distV, &ione, &ione, desca);
        delete [] iwork;
        delete [] nwork;
        delete [] work2;


        // Gather distributed results from distA into global_matrix1
        MainSp->CopyDistArrayToSquareMatrix(global_matrix1, distV, N, desca);
        MainSp->Allreduce(MPI_IN_PLACE, global_matrix1, N*M, MPI_DOUBLE, MPI_SUM);

        delete [] distV;
        delete [] distB;
        delete [] distA;
    }

    // Finally send eigenvalues and vectors to everyone 
    MainSp->BcastRoot(global_matrix1, N * M, MPI_DOUBLE);
    MainSp->BcastRoot(eigs, M, MPI_DOUBLE);

    // Repack V1 into V. Only need the first M eigenvectors
    for(int i = 0;i < M;i++) {
        for(int j = 0;j < N;j++) {
            V[j + i*ld] = global_matrix1[j + i*N];
        }
    }

    return 0;

#endif

}


#if MAGMA_LIBS
template int GeneralDiagMagma<double>(double *A, double *B, double *eigs, double *V, int N, int M, int ld);
template int GeneralDiagMagma<std::complex<double>>(std::complex<double> *A, std::complex<double> *B, double *eigs, std::complex<double> *V,
int N, int M, int ld);


template <typename KpointType>
int GeneralDiagMagma(KpointType *A, KpointType *B, double *eigs, KpointType *V, int N, int M, int ld)
{

    int info = 0;
    //bool use_folded = ((ct.use_folded_spectrum && (ct.scf_steps > 6)) || (ct.use_folded_spectrum && (ct.runflag == RESTART)));

    if(N < M) throw RmgFatalException() << "M must be >= N in " << __FILE__ << " at line " << __LINE__ << "\n";
    if(pct.is_local_master) {

        int *ifail = new int[N];
        int liwork = 6*N;
        int *iwork = new int[liwork];
        int eigs_found;
        double vx = 0.0;
        double tol = 1.0e-15;

        KpointType *Asave = new KpointType[N*ld];
        KpointType *Bsave = new KpointType[N*ld];
        for(int i = 0;i < N*ld;i++) Asave[i] = A[i];
        for(int i = 0;i < N*ld;i++) Bsave[i] = B[i];

        int itype = 1, ione = 1;
        int lwork = 3 * N * N + 6 * N;
        //int lwork = 6 * N * N + 6 * N + 2;
        double *work = (double *)GpuMallocHost(lwork * sizeof(KpointType));

        if(M == N) {
            magma_dsygvd(itype, MagmaVec, MagmaLower, N, (double *)A, ld, (double *)B, ld, eigs, work, lwork, iwork, liwork, &info);
            for(int ix=0;ix < N*ld;ix++) V[ix] = A[ix];
        }
        else if(N > M) {
            magma_dsygvdx (itype, MagmaVec, MagmaRangeI, MagmaLower, N, (double *)A, ld, (double *)B, ld,
                    vx, vx, ione, M,  &eigs_found, eigs, work, lwork, iwork, liwork, &info);
            for(int ix=0;ix < N*ld;ix++) V[ix] = A[ix];
        }


        GpuFreeHost(work);

        for(int i=0;i < N*ld;i++) A[i] = Asave[i];
        for(int i=0;i < N*ld;i++) B[i] = Bsave[i];
        delete [] Bsave;
        delete [] Asave;
        delete [] iwork;
        delete [] ifail;

    } // end if pct.is_local_master

    // If only one proc on this host participated broadcast results to the rest
    if(pct.procs_per_host > 1) {
        int factor = 2;
        if(ct.is_gamma) factor = 1;
        MPI_Bcast(V, factor * M*ld, MPI_DOUBLE, 0, pct.local_comm);
        MPI_Bcast(eigs, M, MPI_DOUBLE, 0, pct.local_comm);
        MPI_Bcast(&info, 1, MPI_INT, 0, pct.local_comm);
    }
    return info;
}
#endif
