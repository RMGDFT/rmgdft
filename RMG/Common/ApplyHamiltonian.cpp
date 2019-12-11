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
#include "rmg_complex.h"


template double ApplyHamiltonian<double,float>(Kpoint<double> *, float *, float *, double *, double *, double *);
template double ApplyHamiltonian<double,double>(Kpoint<double> *, double *, double *, double *, double *, double *);

template double ApplyHamiltonian<std::complex<double>,std::complex<float>>(Kpoint<std::complex<double>> *, std::complex<float> *, 
                             std::complex<float> *, double *, double *, std::complex<double> *);
template double ApplyHamiltonian<std::complex<double>,std::complex<double>>(Kpoint<std::complex<double>> *, std::complex<double> *, 
                             std::complex<double> *, double *, double *, std::complex<double> *);

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
template <typename KpointType, typename CalcType>
double ApplyHamiltonian (Kpoint<KpointType> *kptr, CalcType * __restrict__ psi, CalcType * __restrict__ h_psi, double * __restrict__ vtot, double *vxc_psi, KpointType * __restrict__ nv)
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
    fd_diag = ApplyAOperator<CalcType>(psi, h_psi, dimx, dimy, dimz, gridhx, gridhy, gridhz, ct.kohn_sham_fd_order, kptr->kp.kvec);

    // Factor of -0.5 and add in potential terms
    double tmag(0.5*kptr->kp.kmag);
    for(int idx = 0;idx < pbasis;idx++){ 
        h_psi[idx] = -0.5 * h_psi[idx] + nv[idx] + (vtot[idx] + tmag)*psi[idx];
    }


    if(ct.noncoll)
    {
        ApplyAOperator<CalcType>(&psi[pbasis], &h_psi[pbasis], dimx, dimy, dimz, gridhx, gridhy, gridhz, ct.kohn_sham_fd_order, kptr->kp.kvec);
        for(int idx = 0;idx < pbasis;idx++){ 
            h_psi[idx+pbasis] = -0.5 * h_psi[idx+pbasis] + nv[idx+pbasis] + (vtot[idx] + tmag)*psi[idx+pbasis];
        }

        double *vxc_x = &vxc_psi[pbasis];
        double *vxc_y = &vxc_psi[2*pbasis];
        double *vxc_z = &vxc_psi[3*pbasis];

        // Needed for all of the template variations. Non-complex variants are never actually used but are needed to keep
        // the compiler from throwing errors.
        typedef typename std::conditional_t< std::is_same<CalcType, double>::value, std::complex<double>,
                         std::conditional_t< std::is_same<CalcType, std::complex<double>>::value, std::complex<double>,
                         std::conditional_t< std::is_same<CalcType, std::complex<float>>::value, std::complex<float>, std::complex<float> >>> nctype_t;
        nctype_t *a_psi_C = (nctype_t *)h_psi;
        nctype_t *psi_C = (nctype_t *)psi;

        for(int idx = 0; idx < pbasis; idx++)
        {
            a_psi_C[idx] += psi_C[idx] * std::complex<double>(vxc_z[idx], 0.0);
            a_psi_C[idx] += psi_C[idx+pbasis] * std::complex<double>(vxc_x[idx], -vxc_y[idx]);
            a_psi_C[idx + pbasis] += - psi_C[idx + pbasis] * std::complex<double>(vxc_z[idx], 0.0);
            a_psi_C[idx + pbasis] += psi_C[idx] * std::complex<double>(vxc_x[idx], vxc_y[idx]);
        }

    }

    return fd_diag;
}
