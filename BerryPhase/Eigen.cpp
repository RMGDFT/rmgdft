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

#if MKLBLAS_SET_NUM_THREADS
#include <mkl_service.h>
#endif


// Solves for the eigenvalues and eigenvectors of a Hermitian matrix.
//
//   A*v = lambda*v
//
//   INPUT:  A  =  N by N matrix, distributed as Scalapack class Sp
//   INPUT:  V  =  N by N matrix, distributed as Scalapack class Sp
//   INPUT: ld  =  leading dimension of A,V
//   INPUT:  M  =  number of eigenvectors to solve for
//
// 

template void Eigen<double>(double *A, double *eigs, double *V, int N, int M, Scalapack &Sp);
template void Eigen<std::complex<double>>(std::complex<double> *A, double *eigs, std::complex<double> *V, int N, int M, Scalapack &Sp);

    template <typename KpointType>
void Eigen(KpointType *distA, double *eigs, KpointType *distV, int N, int M, Scalapack &Sp)
{

#if !SCALAPACK_LIBS
    rmg_printf("This version of RMG was not built with Scalapack support. Redirecting to LAPACK.");
#endif

    int ibtype = 1;
    int ione = 1;
    int info = 0;

    if(typeid(KpointType) == typeid(double))
    {

        // Get workspace required
        int izero = 0;
        int nprow = Sp.GetRows();
        int npcol = Sp.GetCols();
        int NB = Sp.GetNB();
        int *desca = Sp.GetDistDesca();
        int NN = std::max( N, NB);
        int NP00 = numroc( &NN, &NB, &izero, &izero, &nprow );
        int NNBmax = std::max(std::max(N, NB), 2);
        int MQ00 = numroc( &NNBmax, &NB, &izero, &izero, &npcol);
        int NNP = std::max(N, std::max(nprow+npcol + 1, 4));
        int lwork = 2 + 5*N + std::max(18*NN, NP00*MQ00 + 2*NB*NB) 
            + (2 + std::ceil((double)N)/(double)(nprow*npcol))*NN;
        lwork *= 2;
        int liwork = 12*NNP + 2*N;
        liwork *= 2;
        // and now solve it 
        double vl = 0.0, vu = 0.0;
        int il = 1, iu = M, eigs_found, eigsv_found;
        double *nwork = new double[lwork];
        int *iwork = new int[liwork];

        pdsyevr("V", "A", "L", &N, (double *)distA, &ione, &ione, desca, &vl, &vu, &il, &iu,
                &eigs_found, &eigsv_found, eigs, (double *)distV, 
                &ione, &ione, desca, nwork, &lwork, iwork, &liwork, &info);

        delete [] iwork;
        delete [] nwork;

    }
    else if(typeid(KpointType) == typeid(std::complex<double>))
    {

        double vl = 0.0, vu = 0.0;
        int il = 1, iu = M, izero = 0, eigs_found, eigsv_found;
        int nprow = Sp.GetRows();
        int npcol = Sp.GetCols();
        int NB = Sp.GetNB();
        int *desca = Sp.GetDistDesca();
        int NN = std::max( N, NB);
        int NP00 = numroc( &NN, &NB, &izero, &izero, &nprow );
        int NNBmax = std::max(std::max(N, NB), 2);
        int MQ00 = numroc( &NNBmax, &NB, &izero, &izero, &npcol);
        int NNP = std::max(N, std::max(nprow+npcol + 1, 4));
        int lwork = N + ( NP00 + MQ00 + NB ) * NB;
        lwork *= 2;
        int lrwork = 2 + 5 * N + std::max(18*N, NP00*MQ00 + 2*NB*NB) +
            (2 + std::ceil( ((double)N)/(double)(nprow*npcol)))*N;
        lrwork *= 2;
        int liwork = 12*NNP + 2*N;
        std::complex<double> *nwork = new std::complex<double>[lwork];
        double *rwork = new double[lrwork];
        int *iwork = new int[liwork];

        pzheevr("V", "A", "L", &N, (std::complex<double> *)distA, &ione, &ione, desca, &vl, &vu, &il, &iu,
                &eigs_found, &eigsv_found, eigs, (std::complex<double> *)distV, &ione, &ione,
                desca, (std::complex<double> *)nwork, &lwork, rwork, &lrwork, iwork, &liwork, &info);

        delete [] iwork;
        delete [] nwork;
        delete [] rwork;


    }

}
