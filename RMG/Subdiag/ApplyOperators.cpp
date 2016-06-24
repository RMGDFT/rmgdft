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
#include "rmg_alloc.h"
#include "rmg_error.h"
#include "State.h"
#include "rmgthreads.h"
#include "RmgTimer.h"
#include "RmgThread.h"
#include "GlobalSums.h"
#include "Subdiag.h"

#include "prototypes.h"
#include "common_prototypes.h"
#include "common_prototypes1.h"
#include "transition.h"

template void ApplyOperators<double>(Kpoint<double> *, int, double *, double *, double *, double *, double *);
template void ApplyOperators<std::complex<double> >(Kpoint<std::complex<double>> *, int, std::complex<double> *, std::complex<double> *, double *, 
              std::complex<double> *, std::complex<double> *);

// Applies A and B operators to one wavefunction
template <typename KpointType>
void ApplyOperators (Kpoint<KpointType> *kptr, int istate, KpointType *a_psi, KpointType *b_psi, double *vtot, KpointType *nv, KpointType *Bns)
{
    BaseGrid *G = kptr->G;
    Lattice *L = &Rmg_L;
    State<KpointType> *sp = &kptr->Kstates[istate];
    KpointType *psi = kptr->Kstates[istate].psi;

    double vel = L->get_omega() / (G->get_NX_GRID(1) * G->get_NY_GRID(1) * G->get_NZ_GRID(1));

    int dimx = G->get_PX0_GRID(1);
    int dimy = G->get_PY0_GRID(1);
    int dimz = G->get_PZ0_GRID(1);
    int pbasis = dimx * dimy * dimz;
    int sbasis = (dimx + 5) * (dimy + 5) * (dimz + 5);  // large enough for up to 10th order central FD operator

    bool potential_acceleration = ((ct.potential_acceleration_constant_step > 0.0) || (ct.potential_acceleration_poisson_step > 0.0));
    potential_acceleration = potential_acceleration & (ct.scf_steps > 0);


    // Apply A operator to psi
    ApplyAOperator (psi, a_psi, "Coarse");


    // Apply B operator to psi
    ApplyBOperator (psi, b_psi, "Coarse");

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

    }


    // Generate 2*V*psi
    KpointType *sg_twovpsi_t = new KpointType[sbasis];
    if(potential_acceleration) {
        int offset = (sp->istate / kptr->dvh_skip) * pbasis;
        CPP_genvpsi (psi, sg_twovpsi_t, &kptr->dvh[offset], (void *)kdr, kptr->kmag, dimx, dimy, dimz);
    }
    else {
        CPP_genvpsi (psi, sg_twovpsi_t, vtot, (void *)kdr, kptr->kmag, dimx, dimy, dimz);
    }

    // B operating on 2*V*psi stored in work
    KpointType *work_t = new KpointType[sbasis];
    ApplyBOperator (sg_twovpsi_t, work_t, "Coarse");

    for(int idx = 0; idx < pbasis; idx++) {

        a_psi[idx] = work_t[idx] - a_psi[idx];

    }

    // Add in non-local which has already had B applied in AppNls
    for(int idx = 0; idx < pbasis; idx++) {

        a_psi[idx] += TWO * nv[idx];

    }

    for (int idx = 0; idx < pbasis; idx++) {

        a_psi[idx] = 0.5 * vel * a_psi[idx];

    }

    // Add in already applied Bns to b_psi for US
    if(!ct.norm_conserving_pp) {

        for(int idx = 0; idx < pbasis; idx++) {

            b_psi[idx] += Bns[idx];

        }

    }

    delete [] kdr;
    delete [] work_t;
    delete [] sg_twovpsi_t;

} // end ApplyOperators


