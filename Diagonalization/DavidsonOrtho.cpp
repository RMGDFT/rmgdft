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
#include "Solvers.h"
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

template void DavidsonOrtho (int, int, int pbasis_noncoll, double *);
template void DavidsonOrtho(int, int, int pbasis_noncoll, std::complex<double> *);

template <typename KpointType>
void DavidsonOrtho(int nbase, int notcon, int pbasis_noncoll, KpointType *psi)
{

    if(!ct.norm_conserving_pp)
    {
        return;
    }
    size_t alloc = (notcon+nbase)*(notcon+nbase)*sizeof(KpointType);
    KpointType *mat = (KpointType *)ct.get_gmatrix(alloc);

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

    if(nbase > 0)
    {
        // ortho to the first nbase states
        RmgGemm(trans_a, trans_n, nbase, notcon, pbasis_noncoll, alphavel, psi, pbasis_noncoll, psi_extra, pbasis_noncoll, zero, mat, nbase);
        BlockAllreduce((double *)mat, (size_t)notcon*(size_t)nbase * (size_t)factor, pct.grid_comm);
        RmgGemm(trans_n, trans_n, pbasis_noncoll, notcon, nbase, mone, psi, pbasis_noncoll, mat, nbase, one, psi_extra, pbasis_noncoll);
    }

    if(!ct.davidson_2stage_ortho) return;

    int st, st1, length, idx, omp_tid;
    KpointType *sarr;
    char *transt = "t";
    char *uplo = "u";
    char *diag = "N";
    char *side = "R";

    KpointType *tarr = new KpointType[notcon];

    if (typeid(KpointType) == typeid(double))
    {
        double rone = 1.0, rzero = 0.0;
        dsyrk( uplo, transt, &notcon, &pbasis_noncoll, &rone, (double *)psi_extra, &pbasis_noncoll,
            &rzero, (double *)mat, &notcon);
    }
    else
    {
        KpointType cone(1.0), czero(0.0);
        RmgGemm(trans_a, trans_n, notcon, notcon, pbasis_noncoll, cone, psi_extra,
                pbasis_noncoll, psi_extra, pbasis_noncoll, czero, mat, notcon);
    }

    /* get the global part */
    length = factor * notcon * notcon;
    MPI_Allreduce(MPI_IN_PLACE, mat, length, MPI_DOUBLE, MPI_SUM, pct.grid_comm);


    /* compute the cholesky factor of the overlap matrix then subtract off projections */
    int info, info_trtri;
    KpointType inv_vel(1.0/sqrt(vel));
    rmg_potrf(uplo, notcon, mat, notcon, &info);
    rmg_trtri(uplo, diag, notcon, mat, notcon, &info_trtri);
    rmg_trmm(side, uplo, "N", diag, pbasis_noncoll, notcon, inv_vel, mat, notcon, psi_extra, pbasis_noncoll);

    if (info != 0)
        throw RmgFatalException() << "Error in " << __FILE__ << " at line " << __LINE__ << ". Matrix not positive definite or argument error. Terminating";

    if (info_trtri != 0)
    throw RmgFatalException() << "Error in " << __FILE__ << " at line " << __LINE__ << "info = " <<info << ". tritri problem";

}
