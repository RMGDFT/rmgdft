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
#include "const.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "rmg_error.h"
#include "RmgTimer.h"
#include "Subdiag.h"
#include "RmgGemm.h"
#include "GpuAlloc.h"
#include "ErrorFuncs.h"
#include "Gpufuncs.h"
#include "blas.h"
#include "GlobalSums.h"
#include "RmgException.h"



#include "transition.h"


// Gram-Schmidt ortho for extra eigenvectors in Davidson solver.
// nbase: number of wavefunctions alreadt orthogonalized
// notcon: number of extra wavefunctions 
// psi: first nbase*pbasis_noncoll will not be changed.
//      next notcon * pbasis_noncoll will be orthogonalized. 
// only work for norm-conserving
// mat: matrix to hold <Psi|psi>, need to be allocated nbase * notconv at least
//

template void MgridOrtho (int, int, int pbasis_noncoll, double *);
template void MgridOrtho(int, int, int pbasis_noncoll, std::complex<double> *);

template <typename KpointType>
void MgridOrtho(int nbase, int notcon, int pbasis_noncoll, KpointType *psi)
{

    if(!ct.norm_conserving_pp)
    {
        return;
    }

#if HIP_ENABLED || CUDA_ENABLED || SYCL_ENABLED
    if(!ct.gmatrix) gpuMallocHost((void **)&ct.gmatrix, notcon * notcon * sizeof(KpointType));
#else
    if(!ct.gmatrix) ct.gmatrix = new KpointType[notcon * notcon];
#endif

    char *trans_t = "t";
    char *trans_n = "n";
    char *trans_c = "c";
    char *trans_a = trans_t;
    int factor = 1;
    if(typeid(KpointType) == typeid(std::complex<double>)) {
        trans_a = trans_c;
        factor = 2;
    }

    double vel = Rmg_L.get_omega() /
        ((double)((size_t)Rmg_G->get_NX_GRID(1) * (size_t)Rmg_G->get_NY_GRID(1) * (size_t)Rmg_G->get_NZ_GRID(1)));
    KpointType alphavel(vel);

    KpointType zero(0.0);
    KpointType one(1.0);
    KpointType mone(-1.0);
    KpointType *psi_extra = &psi[nbase * pbasis_noncoll];
    KpointType *mat = (KpointType *)ct.gmatrix;
    
    int st, st1, length, idx, omp_tid;
    KpointType *sarr;
    char *transt = "t";
    char *uplo = "u";
    char *diag = "N";
    char *side = "R";

    RmgTimer *RT1 = new RmgTimer("MgridOrtho: overlaps");
    KpointType *tarr = new KpointType[notcon];

    if (typeid(KpointType) == typeid(double))
    {
        RmgSyrk( uplo, transt, notcon, pbasis_noncoll, one, psi_extra, pbasis_noncoll,
            zero, mat, notcon);
    }
    else
    {
        KpointType cone(1.0), czero(0.0);
        RmgGemm(trans_a, trans_n, notcon, notcon, pbasis_noncoll, cone, psi_extra,
                pbasis_noncoll, psi_extra, pbasis_noncoll, czero, mat, notcon);
    }
    delete RT1;

    /* get the global part */
    length = factor * notcon * notcon;
    RT1 = new RmgTimer("MgridOrtho: allreduce");
    MPI_Allreduce(MPI_IN_PLACE, mat, length, MPI_DOUBLE, MPI_SUM, pct.grid_comm);
    delete RT1;

    /* compute the cholesky factor of the overlap matrix */
    int info=0, info_trtri=0;
    if (typeid(KpointType) == typeid(double))
    {
        double rone = 1.0/sqrt(vel);
        RT1 = new RmgTimer("MgridOrtho: potrf");
        rmg_potrf(uplo, notcon, mat, notcon, &info);
        delete RT1;
        RT1 = new RmgTimer("MgridOrtho: update");
#if HIP_ENABLED
        double *invmat;
        gpuMalloc((void **)&invmat, notcon * notcon * sizeof(double));

        hipblasFillMode_t hip_fill = HIPBLAS_FILL_MODE_UPPER;
        hipblasSideMode_t hip_side = HIPBLAS_SIDE_RIGHT;
        hipblasDiagType_t hip_diag = HIPBLAS_DIAG_NON_UNIT;
        hipblasOperation_t hip_trans = HIPBLAS_OP_N;
        hipblasDtrtri(ct.hipblas_handle, hip_fill, hip_diag, notcon, (double *)mat, notcon,
        invmat, notcon);
        hipblasDtrmm(ct.hipblas_handle, hip_side, hip_fill, hip_trans, 
                     hip_diag, pbasis_noncoll, notcon, &rone, 
                     invmat, notcon, (double *)psi_extra, pbasis_noncoll,
                     (double *)psi_extra, pbasis_noncoll); 
        gpuFree(invmat);

#else
        dtrtri(uplo, diag, &notcon, (double *)mat, &notcon, &info_trtri);
        dtrmm(side, uplo, "N", diag, &pbasis_noncoll, &notcon, &rone, 
             (double *)mat, &notcon, (double *)psi_extra, &pbasis_noncoll); 
#endif
        delete RT1;
    }
    else
    {
        RT1 = new RmgTimer("MgridOrtho: potrf");
        rmg_potrf(uplo, notcon, mat, notcon, &info);
        delete RT1;
        RT1 = new RmgTimer("MgridOrtho: update");
#if HIP_ENABLED
        hipblasDoubleComplex cone = 1.0/sqrt(vel);
        std::complex<double> *invmat;
        gpuMalloc((void **)&invmat, notcon * notcon * sizeof(std::complex<double>));
        hipblasFillMode_t hip_fill = HIPBLAS_FILL_MODE_UPPER;
        hipblasSideMode_t hip_side = HIPBLAS_SIDE_RIGHT;
        hipblasDiagType_t hip_diag = HIPBLAS_DIAG_NON_UNIT;
        hipblasOperation_t hip_trans = HIPBLAS_OP_N;
        hipblasZtrtri(ct.hipblas_handle, hip_fill, hip_diag, notcon, (hipblasDoubleComplex *)mat,
                      notcon, (hipblasDoubleComplex *)invmat, notcon);
        hipblasZtrmm(ct.hipblas_handle, hip_side, hip_fill, hip_trans, 
                     hip_diag, pbasis_noncoll, notcon, &cone, 
                     (hipblasDoubleComplex *)invmat, notcon,
                     (hipblasDoubleComplex *)psi_extra, pbasis_noncoll,
                     (hipblasDoubleComplex *)psi_extra, pbasis_noncoll); 
        gpuFree(invmat);
#else
        std::complex<double> cone = 1.0/sqrt(vel);
        ztrtri(uplo, diag, &notcon, (std::complex<double> *)mat, &notcon, &info_trtri);
        ztrmm(side, uplo, "N", diag, &pbasis_noncoll, &notcon, &cone, (std::complex<double> *)mat, &notcon, (std::complex<double> *)psi_extra, &pbasis_noncoll); 
#endif
        delete RT1;
    }
    if (info != 0)
        throw RmgFatalException() << "Error in " << __FILE__ << " at line " << __LINE__ << ". Matrix not positive definite or argument error. Terminating";

    if (info_trtri != 0)
    throw RmgFatalException() << "Error in " << __FILE__ << " at line " << __LINE__ << "info = " <<info << ". tritri problem";

}
