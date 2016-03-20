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


template void ApplyHamiltonian<double>(Kpoint<double> *, int, double *, double *, double *);
template void ApplyHamiltonian<std::complex<double> >(Kpoint<std::complex<double>> *, int, std::complex<double> *, double *, std::complex<double> *);

// Applies Hamiltonian operator to one orbital
//
//  INPUT
//    kptr   = kpoint object
//    istate = index of orbital from zero
//    vtot   = total local potential on wavefunction grid
//    nv     = Non-local potential applied to this orbital
//  OUTPUT
//    h_psi  = H|psi>
//
template <typename KpointType>
void ApplyHamiltonian (Kpoint<KpointType> *kptr, int istate, KpointType *h_psi, double *vtot, KpointType *nv)
{
    KpointType *psi = kptr->Kstates[istate].psi;
    int pbasis = Rmg_G->get_P0_BASIS(Rmg_G->default_FG_RATIO);

    // Apply Laplacian to psi
    ApplyLaplacian (psi, h_psi, ct.kohn_sham_fd_order, "Fine");

    // Factor of -0.5 and add in potential terms
    for(int idx = 0;idx < pbasis;idx++) h_psi[idx] = -0.5 * h_psi[idx] + nv[idx] + vtot[idx]*psi[idx];

    // if complex orbitals apply gradient to orbital and compute dot products
    std::complex<double> *kdr = new std::complex<double>[pbasis]();

    if(typeid(KpointType) == typeid(std::complex<double>)) {

        KpointType *gx = new KpointType[pbasis];
        KpointType *gy = new KpointType[pbasis];
        KpointType *gz = new KpointType[pbasis];

        ApplyGradient (psi, gx, gy, gz, APP_CI_EIGHT, "Coarse");

        std::complex<double> I_t(0.0, 1.0);
        for(int idx = 0;idx < pbasis;idx++) {

            kdr[idx] = -I_t * (kptr->kvec[0] * (std::complex<double>)gx[idx] +
                               kptr->kvec[1] * (std::complex<double>)gy[idx] +
                               kptr->kvec[2] * (std::complex<double>)gz[idx]);
        }

        delete [] gz;
        delete [] gy;
        delete [] gx;

        for(int idx=0;idx < pbasis;idx++) h_psi[idx] + kdr[idx] + 0.5*kptr->kmag;
    }

}
