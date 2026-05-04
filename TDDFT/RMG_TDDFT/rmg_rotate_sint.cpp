/*
 *
 * Copyright 2026 The RMG Project Developers. See the COPYRIGHT file 
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

#include "rmg_tddft.h"
#include "../Headers/prototypes_tddft.h"
#include "rmg_dev_allocate.h"
#include "rmg_hvector.h"
#include "Kpoint.h"
#include "blas_driver.h"
#include "rmg_gemm.h"


template void rmg::rotate_sint<double, double>(Kpoint<double> *, double *, double *);
template void rmg::rotate_sint<double, std::complex<double>>(Kpoint<double> *, double *, std::complex<double> *);
template void rmg::rotate_sint<std::complex<double>, std::complex<double>>(Kpoint<std::complex<double>> *, std::complex<double> *, std::complex<double> *);
template <typename OrbitalType, typename MatrixType>

void rmg::rotate_sint(
            Kpoint<OrbitalType> *Kptr, OrbitalType *sint, MatrixType *rho_matrix)
{
    char *transt = "t", *transn = "n", *transc = "c";
    char *transa;

    transa = transc;
    if(ct.is_gamma) transa = transt;
    int num_nonloc_ions = Kptr->BetaProjector->get_num_nonloc_ions();
    int pstride = Kptr->BetaProjector->get_pstride();

    // Vector potential at gamma case
    if constexpr (std::is_same_v<OrbitalType, double> && std::is_same_v<MatrixType, std::complex<double>>)
    {
        double rzero(0.0);
        double alpha(1.0);
        OrbitalType *oldsint = Kptr->oldsint_local;
        rmg::hvector<double> real_rho_matrix(ct.num_states*ct.num_states); 
        for(int i=0;i < ct.num_states*ct.num_states;i++) real_rho_matrix[i] = std::real(rho_matrix[i]);
        rmg::gemm (transn, transt, ct.num_states, num_nonloc_ions*pstride, ct.num_states, alpha,
            oldsint, num_nonloc_ions*pstride, real_rho_matrix.data(), ct.num_states, rzero, sint, num_nonloc_ions*pstride);
        //for(int i=0;i < ct.num_states*num_nonloc_ions*pstride;i++) sint[i] *= 0.5;
    }

}
