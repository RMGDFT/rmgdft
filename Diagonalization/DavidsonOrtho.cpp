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



#include "transition.h"


// Gram-Schmidt ortho for extra eigenvectors in Davidson solver.
// nbase: number of wavefunctions alreadt orthogonalized
// notcon: number of extra wavefunctions 
// psi: first nbase*pbasis_noncoll will not be changed.
//      next notcon * pbasis_noncoll will be orthogonalized. 
// only work for norm-conserving
// mat: matrix to hold <Psi|psi>, need to be allocated nbase * notconv at least
//

template void DavidsonOrtho (int, int, int pbasis_noncoll, double *, double *);
template void DavidsonOrtho(int, int, int pbasis_noncoll, std::complex<double> *, std::complex<double> *);

template <typename KpointType>
void DavidsonOrtho(int nbase, int notcon, int pbasis_noncoll, KpointType *psi, KpointType *mat)
{

    if(!ct.norm_conserving_pp)
    {
        rmg_error_handler(__FILE__, __LINE__, "only support norm-conserving pp in DavidsonOrtho now\n");
    }
    if(nbase == 0)
    {
        rmg_error_handler(__FILE__, __LINE__, "nbase cannot be zero\n");
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

    double vel = Rmg_L.get_omega() /
        ((double)((size_t)Rmg_G->get_NX_GRID(1) * (size_t)Rmg_G->get_NY_GRID(1) * (size_t)Rmg_G->get_NZ_GRID(1)));
    KpointType alphavel(vel);

    KpointType zero(0.0);
    KpointType one(1.0);
    KpointType mone(-1.0);
    KpointType *psi_extra = &psi[nbase * pbasis_noncoll];

    // ortho to the first nbase states
    RmgGemm(trans_a, trans_n, nbase, notcon, pbasis_noncoll, alphavel, psi, pbasis_noncoll, psi_extra, pbasis_noncoll, zero, mat, nbase);
    BlockAllreduce((double *)mat, (size_t)notcon*(size_t)nbase * (size_t)factor, pct.grid_comm);
    RmgGemm(trans_n, trans_n, pbasis_noncoll, notcon, nbase, mone, psi, pbasis_noncoll, mat, nbase, one, psi_extra, pbasis_noncoll);

    // ortho for the remainl noncon states

    double norm;
    int pbasis_c = pbasis_noncoll * factor;
    int ione = 1;
    norm =  ddot(&pbasis_c, (double *)psi_extra, &ione, (double *)psi_extra, &ione);
    MPI_Allreduce(MPI_IN_PLACE, &norm, ione, MPI_DOUBLE, MPI_SUM, pct.grid_comm);
    norm = 1.0/sqrt(norm * vel);
    dscal(&pbasis_c, &norm, (double *)psi_extra, &ione);

    for (int st = 1; st < notcon; st++)
    {
        // mat = <psi_extra[0:st-1] |psi_extr[st]>
        RmgGemm(trans_a, trans_n, st, ione, pbasis_noncoll, alphavel, psi_extra, pbasis_noncoll, &psi_extra[st*pbasis_noncoll], pbasis_noncoll, zero, mat, st);
        BlockAllreduce((double *)mat, (size_t)st * (size_t)factor, pct.grid_comm);
        RmgGemm(trans_n, trans_n, pbasis_noncoll, ione, st, mone, psi_extra, pbasis_noncoll, mat, st, one, &psi_extra[st*pbasis_noncoll], pbasis_noncoll);
        norm =  ddot(&pbasis_c, (double *)&psi_extra[st*pbasis_noncoll], &ione, (double *)&psi_extra[st*pbasis_noncoll], &ione);
        MPI_Allreduce(MPI_IN_PLACE, &norm, ione, MPI_DOUBLE, MPI_SUM, pct.grid_comm);
        norm = 1.0/sqrt(norm * vel);
        dscal(&pbasis_c, &norm, (double *)&psi_extra[st*pbasis_noncoll], &ione);
    }


    /*
    RmgGemm(trans_a, trans_n, nbase+notcon, nbase+notcon, pbasis_noncoll, alphavel, psi, pbasis_noncoll, psi, pbasis_noncoll, zero, mat, nbase+notcon);
    BlockAllreduce((double *)mat, (size_t)(nbase+notcon)*(size_t)(nbase+notcon) * (size_t)factor, pct.grid_comm);
    for(int i = 0; i < nbase + notcon; i++) 
        for(int j = 0; j < nbase + notcon; j++) 
        {
            if(std::abs(mat[i * (nbase + notcon) + j]) > 1.0e-5)
            rmg_printf("\n ortho? %d %d  %e %e", i,j, mat[i *(nbase+notcon) + j]);
        }

    rmg_error_handler(__FILE__, __LINE__, "only support norm-conserving pp in DavidsonOrtho now\n");
    */

}
