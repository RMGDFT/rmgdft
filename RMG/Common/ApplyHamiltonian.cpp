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
#include "FiniteDiff.h"
#include "const.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "State.h"
#include "GlobalSums.h"
#include "Subdiag.h"
#include "Solvers.h"

#include "transition.h"


template double ApplyHamiltonian<double>(Kpoint<double> *, double *, double *, double *, double *);
template double ApplyHamiltonian<std::complex<double> >(Kpoint<std::complex<double>> *, std::complex<double> *, 
                             std::complex<double> *, double *, std::complex<double> *);

// Applies Hamiltonian operator to one orbital
//
//  INPUT
//    kptr   = kpoint object
//    psi    = the orbital
//    vtot   = total local potential on wavefunction grid
//    nv     = Non-local potential applied to this orbital
//  OUTPUT
//    h_psi  = H|psi>
//
template <typename KpointType>
double ApplyHamiltonian (Kpoint<KpointType> *kptr, KpointType * __restrict__ psi, KpointType * __restrict__ h_psi, double * __restrict__ vtot, KpointType * __restrict__ nv)
{
    int pbasis = kptr->pbasis;
    double fd_diag;

    int density = 1;
    int dimx = kptr->G->get_PX0_GRID(density);
    int dimy = kptr->G->get_PY0_GRID(density);
    int dimz = kptr->G->get_PZ0_GRID(density);
    double gridhx = kptr->G->get_hxgrid(density);
    double gridhy = kptr->G->get_hygrid(density);
    double gridhz = kptr->G->get_hzgrid(density);
    fd_diag = ApplyAOperator<KpointType>(psi, h_psi, dimx, dimy, dimz, gridhx, gridhy, gridhz, ct.kohn_sham_fd_order, kptr->kp.kvec);

    // Factor of -0.5 and add in potential terms
    KpointType tmag(0.5*kptr->kp.kmag);
    for(int idx = 0;idx < pbasis;idx++){ 
        h_psi[idx] = -0.5 * h_psi[idx] + nv[idx] + (vtot[idx] + tmag)*psi[idx];
    }

    return fd_diag;
}
