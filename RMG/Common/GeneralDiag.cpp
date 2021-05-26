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

#if CUDA_ENABLED
#if MAGMA_LIBS
    #include <magma.h>
#endif
#endif

#if MKLBLAS_SET_NUM_THREADS
    #include <mkl_service.h>
#endif


// Kludge for now
extern Scalapack *MainSp;

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

// dsygvdx routines are not available in cusolver library until version 10.1 or later
// so redirect to scalapack
#if CUDA_ENABLED
  #if (CUDA_VERSION_MAJOR < 10) || ((CUDA_VERSION_MAJOR == 10) && (CUDA_VERSION_MINOR < 1))
        case SUBDIAG_CUSOLVER:
  #endif
#endif
        case SUBDIAG_SCALAPACK:
            info = GeneralDiagScaLapack(A, B, eigs, V, N, M, ld);
            break;
#if MAGMA_LIBS && CUDA_ENABLED
        case SUBDIAG_MAGMA:
                info = GeneralDiagMagma(A, B, eigs, V, N, M, ld);
#endif
#if CUDA_ENABLED
  #if (CUDA_VERSION_MAJOR > 10) || ((CUDA_VERSION_MAJOR == 10) && (CUDA_VERSION_MINOR > 0))
        case SUBDIAG_CUSOLVER:
            info = GeneralDiagCusolver(A, B, eigs, V, N, M, ld);
            break;
  #endif
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

    if(pct.is_local_master) {

        // Increase the resources available to this proc since the others on the local node
        // will be idle
        int nthreads = ct.OMP_THREADS_PER_NODE;
        if(pct.procs_per_host > 1) nthreads = pct.ncpus;
        nthreads = std::min(nthreads, 4);
        omp_set_num_threads(nthreads);
#if OPENBLAS_SET_NUM_THREADS
        openblas_set_num_threads(nthreads);
#endif
#if MKLBLAS_SET_NUM_THREADS
        mkl_set_num_threads_local(ct.OMP_THREADS_PER_NODE);
#endif


        {
            if(N < M) throw RmgFatalException() << "M must be >= N in " << __FILE__ << " at line " << __LINE__ << "\n";
            
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



            if(ct.is_gamma)
            {
                if(M == N) {
                    dsygvd (&itype, "V", "L", &N, (double *)A, &ld, (double *)B, &ld, eigs, work2, &lwork, iwork, &liwork, &info);
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
                    zhegvd (&itype, "V", "U", &N, (double *)A, &ld, (double *)B, &ld, eigs, work2, &lwork, &work2[lwork], &lwork, iwork, &liwork, &info);
                    for(int ix=0;ix < N*ld;ix++) V[ix] = A[ix];
                }
                else if(N > M) {
                    zhegvx (&itype, "V", "I", "U", &N, (double *)A, &ld, (double *)B, &ld,
                            &vx, &vx, &ione, &M,  &tol, &eigs_found, eigs, (double *)V, &ld, work2,
                            &lwork, &work2[lwork], iwork, ifail, &info);

                    //printf("\n %d infooo ", info);
                }

            }

            if(info) 
            {
                for(int i=0;i < N*ld;i++) B[i] = Bsave[i];
                int info1;
                double *rwork = new double[3*N];
                zheev("V", "L", &N, (double *)B, &ld,(double *)A, work2,&lwork, rwork, &info1);
                delete [] rwork;

                if(pct.gridpe == 0)
                {
                    printf("\n WARNING: GeneralDiag.cpp, failed, the eig of overlap matrix as follow");
                    for(int i = 0; i < 5; i++) printf("\n eig %d %d %d %e ", info, pct.spinpe, i, std::real(A[i]));
                }
                fflush(NULL);
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
        omp_set_num_threads(ct.OMP_THREADS_PER_NODE);
#if OPENBLAS_SET_NUM_THREADS
        openblas_set_num_threads(ct.OMP_THREADS_PER_NODE);
#endif
#if MKLBLAS_SET_NUM_THREADS
        mkl_set_num_threads_local(ct.OMP_THREADS_PER_NODE);
#endif


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


template int GeneralDiagScaLapack<double>(double *A, double *B, double *eigs, double *V, int N, int M, int ld);
template int GeneralDiagScaLapack<std::complex<double>>(std::complex<double> *A, std::complex<double> *B, double *eigs, std::complex<double> *V,
int N, int M, int ld);

template <typename KpointType>
int GeneralDiagScaLapack(KpointType *A, KpointType *B, double *eigs, KpointType *V, int N, int M, int ld)
{

#if !SCALAPACK_LIBS
    rmg_printf("This version of RMG was not built with Scalapack support. Redirecting to LAPACK.");
    return GeneralDiagLapack(A, B, eigs, V, N, M, ld);
#else
    KpointType *global_matrix1 = (KpointType *)RmgMallocHost(ct.max_states * ct.max_states * sizeof(KpointType));

    int info = 0;
    // Create 1 scalapack instance per grid_comm. We use a static Scalapack here since initialization on large systems is expensive
    if(!MainSp) {
        // Need some code here to decide how to set the number of scalapack groups but for now use just 1
        int scalapack_groups = 1;
        int last = !ct.use_folded_spectrum;
        MainSp = new Scalapack(scalapack_groups, pct.thisimg, ct.images_per_node, ct.max_states,
                ct.scalapack_block_factor, last, pct.grid_comm);

    }

    bool participates = MainSp->Participates();

    if (participates) {


        // Get dist length and desca
        int desca[DLEN];
        int dist_length = MainSp->ComputeDesca(N, N, desca);

        //desca = MainSp->GetDistDesca();
        MainSp->ComputeDesca(N, N, desca);

        // Allocate distributed matrices 
        KpointType *distA = new KpointType[dist_length * sizeof(KpointType)];
        KpointType *distB = new KpointType[dist_length * sizeof(KpointType)];
        KpointType *distV = new KpointType[dist_length * sizeof(KpointType)];


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
        int ione = 1;

        if(typeid(KpointType) == typeid(double))
        {
            double scale=1.0, rone = 1.0;
            pdpotrf("L", &N, (double *)distB,  &ione, &ione, desca,  &info);

            // Get pdsyngst_ workspace
            int lwork = -1;
            double lwork_tmp;
            pdsyngst(&ibtype, "L", &N, (double *)distA, &ione, &ione, desca,
                    (double *)distB, &ione, &ione, desca, &scale, &lwork_tmp, &lwork, &info);
            lwork = 2*(int)lwork_tmp;
            double *work2 = new double[lwork];

            pdsyngst(&ibtype, "L", &N, (double *)distA, &ione, &ione, desca,
                    (double *)distB, &ione, &ione, desca, &scale, work2, &lwork, &info);

            // Get workspace required
            lwork = -1;
            int liwork=-1;
            int liwork_tmp;
            pdsyevd("V", "L", &N, (double *)distA, &ione, &ione, desca,
                    eigs, (double *)distV, &ione, &ione, desca, &lwork_tmp, &lwork, &liwork_tmp, &liwork, &info);
            lwork = 16*(int)lwork_tmp;
            liwork = 16*N;
            double *nwork = new double[lwork];
            int *iwork = new int[liwork];

            // and now solve it 
            pdsyevd("V", "L", &N, (double *)distA, &ione, &ione, desca,
                    eigs, (double *)distV, &ione, &ione, desca, nwork, &lwork, iwork, &liwork, &info);

            pdtrsm("Left", "L", "T", "N", &N, &N, &rone, (double *)distB, &ione, &ione, desca,
                    (double *)distV, &ione, &ione, desca);
            delete [] iwork;
            delete [] nwork;
            delete [] work2;

        }
        else if(typeid(KpointType) == typeid(std::complex<double>))
        {
            double scale=1.0, rone[2] = {1.0, 0.0};
            pzpotrf("L", &N, (double *)distB,  &ione, &ione, desca,  &info);
            if(info) return info;

            pzhegst(&ibtype, "L", &N, (double *)distA, &ione, &ione, desca,
                    (double *)distB, &ione, &ione, desca, &scale, &info);

            if(info) return info;
            // Get workspace required
            int lwork = -1, liwork=-1, lrwork=-1;
            double lwork_tmp[2], lrwork_tmp;
            int liwork_tmp;
            pzheevd("V", "L", &N, (double *)distA, &ione, &ione, desca,
                    eigs, (double *)distV, &ione, &ione, desca, lwork_tmp, &lwork, &lrwork_tmp, &lrwork, &liwork_tmp, &liwork, &info);
            lwork = (int)lwork_tmp[0]+1;
            liwork = liwork_tmp + 1;
            lrwork = 2*(int)lrwork_tmp;
            double *rwork = new double[lrwork];
            double *nwork = new double[lwork*2];
            int *iwork = new int[liwork];

            // and now solve it
            pzheevd("V", "L", &N, (double *)distA, &ione, &ione, desca,
                    eigs, (double *)distV, &ione, &ione, desca, nwork, &lwork, (double *)rwork, &lrwork, iwork, &liwork, &info);
            if(info) return info;

            pztrsm("Left", "L", "C", "N", &N, &N, rone, (double *)distB, &ione, &ione, desca,
                    (double *)distV, &ione, &ione, desca);

            delete [] iwork;
            delete [] nwork;
            delete [] rwork;


        }

        // Gather distributed results from distA into global_matrix1
        int fac = sizeof(KpointType)/sizeof(double);
        MainSp->CopyDistArrayToSquareMatrix(global_matrix1, distV, N, desca);
        MainSp->Allreduce(MPI_IN_PLACE, global_matrix1, N*M*fac, MPI_DOUBLE, MPI_SUM);

        delete [] distV;
        delete [] distB;
        delete [] distA;
    }

    // Finally send eigenvalues and vectors to everyone 
    int fac = sizeof(KpointType)/sizeof(double);
    MainSp->BcastRoot(global_matrix1, N * M * fac, MPI_DOUBLE);
    MainSp->BcastRoot(eigs, M, MPI_DOUBLE);

    // Repack V1 into V. Only need the first M eigenvectors
    for(int i = 0;i < M;i++) {
        for(int j = 0;j < N;j++) {
            V[j + i*ld] = global_matrix1[j + i*N];
        }
    }

    RmgFreeHost(global_matrix1);
    return info;

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

    if(N < M) throw RmgFatalException() << "M must be >= N in " << __FILE__ << " at line " << __LINE__ << "\n";
    if(pct.is_local_master) {

        int *ifail = new int[N];
        int liwork = 6*N;
        int *iwork = new int[liwork];
        int eigs_found;
        double vx = 0.0;

        KpointType *Asave = new KpointType[N*ld];
        KpointType *Bsave = new KpointType[N*ld];
        for(int i = 0;i < N*ld;i++) Asave[i] = A[i];
        for(int i = 0;i < N*ld;i++) Bsave[i] = B[i];

        int itype = 1, ione = 1;
        int lwork = 3 * N * N + 6 * N;
        //int lwork = 6 * N * N + 6 * N + 2;
        double *work = (double *)RmgMallocHost(lwork * sizeof(KpointType));

        if(M == N) {
            magma_dsygvd(itype, MagmaVec, MagmaLower, N, (double *)A, ld, (double *)B, ld, eigs, work, lwork, iwork, liwork, &info);
            for(int ix=0;ix < N*ld;ix++) V[ix] = A[ix];
        }
        else if(N > M) {
            magma_dsygvdx (itype, MagmaVec, MagmaRangeI, MagmaLower, N, (double *)A, ld, (double *)B, ld,
                    vx, vx, ione, M,  &eigs_found, eigs, work, lwork, iwork, liwork, &info);
            for(int ix=0;ix < N*ld;ix++) V[ix] = A[ix];
        }


        RmgFreeHost(work);

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

#if CUDA_ENABLED
#if (CUDA_VERSION_MAJOR > 10) || ((CUDA_VERSION_MAJOR == 10) && (CUDA_VERSION_MINOR > 0))
template int GeneralDiagCusolver<double>(double *A, double *B, double *eigs, double *V, int N, int M, int ld);
template int GeneralDiagCusolver<std::complex<double>>(std::complex<double> *A, std::complex<double> *B, double *eigs, std::complex<double> *V,
        int N, int M, int ld);

    template <typename KpointType>
int GeneralDiagCusolver(KpointType *A, KpointType *B, double *eigs, KpointType *V, int N, int M, int ld)
{

    // Redirect to CPU routines until we get this coded
    int info;
    if(!ct.is_gamma) {
        info = GeneralDiagLapack(A, B, eigs, V, N, M, ld);
    }

    else {
        info = GeneralDiagScaLapack((double *)A, (double *)B, eigs, (double *)V, N, M, ld);
    }
    return info;
}
#endif
#endif

