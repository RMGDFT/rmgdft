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


template void rmg::rotate_sint<double>(Kpoint<double> *, double *, double *);
template void rmg::rotate_sint<std::complex<double>>(Kpoint<std::complex<double>> *, std::complex<double> *, std::complex<double> *);
template <typename OrbitalType>

void rmg::rotate_sint(
        Kpoint<OrbitalType> *Kptr, OrbitalType *sint, OrbitalType *rho_matrix)
{
    // No sint to rotate in this case
    if(ct.internal_pseudo_type == ALL_ELECTRON) return;

    char *transt = "t", *transn = "n", *transc = "c";
    char *transa;

    transa = transc;
    if(ct.is_gamma) transa = transt;
    int num_nonloc_ions = Kptr->BetaProjector->get_num_nonloc_ions();
    int pstride = Kptr->BetaProjector->get_pstride();

    OrbitalType rzero(0.0);
    OrbitalType alpha(1.0);

    OrbitalType *oldsint = Kptr->oldsint_local;
    rmg::gemm (transn, transa, num_nonloc_ions*pstride, ct.num_states, ct.num_states, alpha,
            oldsint, num_nonloc_ions*pstride, rho_matrix, ct.num_states, rzero, sint, num_nonloc_ions*pstride);

}
