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


template double ApplyHamiltonian<double,float>(Kpoint<double> *, int, float *, float *, double *, double *, double *, bool);
template double ApplyHamiltonian<double,double>(Kpoint<double> *, int, double *, double *, double *, double *, double *, bool);

template double ApplyHamiltonian<std::complex<double>,std::complex<float>>(Kpoint<std::complex<double>> *, int, std::complex<float> *, 
                             std::complex<float> *, double *, double *, std::complex<double> *, bool);
template double ApplyHamiltonian<std::complex<double>,std::complex<double>>(Kpoint<std::complex<double>> *, int, std::complex<double> *, 
                             std::complex<double> *, double *, double *, std::complex<double> *, bool);

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
double ApplyHamiltonian (Kpoint<KpointType> *kptr, int istate, CalcType * __restrict__ psi, CalcType * __restrict__ h_psi, double * __restrict__ vtot, double *vxc_psi, KpointType * __restrict__ nv, bool potential_acceleration)
{
    int pbasis = kptr->pbasis;
    double fd_diag;
    double *veff = NULL;

    potential_acceleration = potential_acceleration && (ct.potential_acceleration_constant_step > 0.0);
    potential_acceleration = potential_acceleration & (ct.scf_steps > 0);

    int density = 1;
    int dimx = kptr->G->get_PX0_GRID(density);
    int dimy = kptr->G->get_PY0_GRID(density);
    int dimz = kptr->G->get_PZ0_GRID(density);
    double gridhx = kptr->G->get_hxgrid(density);
    double gridhy = kptr->G->get_hygrid(density);
    double gridhz = kptr->G->get_hzgrid(density);
    fd_diag = ApplyAOperator<CalcType>(psi, h_psi, dimx, dimy, dimz, gridhx, gridhy, gridhz, ct.kohn_sham_fd_order, kptr->kp.kvec);

    if(potential_acceleration) {
        int active_threads = ct.MG_THREADS_PER_NODE;
        if(ct.mpi_queue_mode && (active_threads > 1)) active_threads--;
        int offset = (istate / kptr->dvh_skip) * pbasis;
        int my_pe_x, my_pe_y, my_pe_z;
        kptr->G->pe2xyz(pct.gridpe, &my_pe_x, &my_pe_y, &my_pe_z);
        int my_pe_offset = my_pe_x % pct.coalesce_factor;
        veff = &kptr->dvh[offset*pct.coalesce_factor + my_pe_offset*pbasis];
    }
    else {
        veff = vtot;
    }

    // Factor of -0.5 and add in potential terms
    double tmag(0.5*kptr->kp.kmag);
    for(int idx = 0;idx < pbasis;idx++){ 
        h_psi[idx] = -0.5 * h_psi[idx] + nv[idx] + (veff[idx] + tmag)*psi[idx];
    }


    if(ct.noncoll)
    {
        ApplyAOperator<CalcType>(&psi[pbasis], &h_psi[pbasis], dimx, dimy, dimz, gridhx, gridhy, gridhz, ct.kohn_sham_fd_order, kptr->kp.kvec);
        for(int idx = 0;idx < pbasis;idx++){ 
            h_psi[idx+pbasis] = -0.5 * h_psi[idx+pbasis] + nv[idx+pbasis] + (veff[idx] + tmag)*psi[idx+pbasis];
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
