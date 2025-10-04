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
#include "RmgMatrix.h"
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

template void DavidsonOrtho (int, int, int pbasis_noncoll, double *, bool);
template void DavidsonOrtho(int, int, int pbasis_noncoll, std::complex<double> *, bool);

template <typename KpointType>
void DavidsonOrtho(int nbase, int notcon, int pbasis_noncoll, KpointType *psi, bool dostage2)
{

    if(!ct.norm_conserving_pp)
    {
        return;
    }

    char *trans_t = "t";
    char *trans_n = "n";
    char *trans_c = "c";
    char *trans_a = trans_t;
    int factor = 1;
    if(typeid(KpointType) == typeid(std::complex<double>)) {
        trans_a = trans_c;
        factor = 2;
    }

    size_t tlength = ((notcon + 2) * notcon / 2);
    size_t alloc = (notcon+nbase)*(notcon+nbase);
    KpointType *mat = (KpointType *)ct.get_gmatrix((alloc + tlength + 8192)*sizeof(KpointType));
    size_t offset = 4096 * (alloc / 4096 + 1);
    KpointType *tmat = mat + offset;

    double vel = Rmg_L.get_omega() /
        ((double)((size_t)Rmg_G->get_NX_GRID(1) * (size_t)Rmg_G->get_NY_GRID(1) * (size_t)Rmg_G->get_NZ_GRID(1)));
    KpointType alphavel(vel);

    KpointType zero(0.0);
    KpointType one(1.0);
    KpointType mone(-1.0);
    KpointType *psi_extra = &psi[nbase * pbasis_noncoll];

    if(nbase > 0)
    {
        RmgTimer RT("MgridOrtho: 1st stage");
        // ortho to the first nbase states
        RmgGemm(trans_a, trans_n, nbase, notcon, pbasis_noncoll, alphavel, psi, pbasis_noncoll, psi_extra, pbasis_noncoll, zero, mat, nbase);
        BlockAllreduce((double *)mat, (size_t)notcon*(size_t)nbase * (size_t)factor, pct.grid_comm);
        RmgGemm(trans_n, trans_n, pbasis_noncoll, notcon, nbase, mone, psi, pbasis_noncoll, mat, nbase, one, psi_extra, pbasis_noncoll);
    }

    if(!dostage2) return;

    RmgTimer *RT1 = new RmgTimer("MgridOrtho: 1st stage update");
    char *transt = "t";
    char *uplo = "u";
    char *diag = "N";
    char *side = "R";

    if constexpr (std::is_same_v<KpointType, double>)
    {
        RmgSyrk( uplo, transt, notcon, pbasis_noncoll, one, psi_extra, pbasis_noncoll,
            zero, mat, notcon);
    }
    if constexpr (std::is_same_v<KpointType, std::complex<double>>)
    {
        RmgGemm(trans_a, trans_n, notcon, notcon, pbasis_noncoll, one, psi_extra,
                pbasis_noncoll, psi_extra, pbasis_noncoll, zero, mat, notcon);
    }
    delete RT1;

    /* get the global part */
    RT1 = new RmgTimer("MgridOrtho: allreduce");
    int length = factor * (notcon + 2) * notcon / 2;
    PackSqToTr("U", notcon, mat, notcon, tmat);
    BlockAllreduce((double *)tmat, length, pct.grid_comm);
    UnPackSqToTr("U", notcon, mat, notcon, tmat);
    delete RT1;

    /* compute the cholesky factor of the overlap matrix then subtract off projections */
    RT1 = new RmgTimer("MgridOrtho: cholesky");
    int info, info_trtri;
    KpointType inv_vel(1.0/sqrt(vel));
    rmg_potrf(uplo, notcon, mat, notcon, &info);
    delete RT1;
    RT1 = new RmgTimer("MgridOrtho: inverse");
    rmg_trtri(uplo, diag, notcon, mat, notcon, &info_trtri);
    delete RT1;
    RT1 = new RmgTimer("MgridOrtho: 2nd stage update");
    rmg_trmm(side, uplo, "N", diag, pbasis_noncoll, notcon, inv_vel, mat, notcon, psi_extra, pbasis_noncoll);
    delete RT1;

    if (info != 0)
        throw RmgFatalException() << "Error in " << __FILE__ << " at line " << __LINE__ << ". Matrix not positive definite or argument error. Terminating";

    if (info_trtri != 0)
    throw RmgFatalException() << "Error in " << __FILE__ << " at line " << __LINE__ << "info = " <<info << ". tritri problem";

}
