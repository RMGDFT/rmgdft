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
#include "rmg_error.h"
#include "State.h"
#include "rmgthreads.h"
#include "RmgTimer.h"
#include "RmgThread.h"
#include "GlobalSums.h"
#include "Subdiag.h"

#include "prototypes_rmg.h"
#include "common_prototypes.h"
#include "common_prototypes1.h"
#include "transition.h"

template void ApplyOperators<double>(Kpoint<double> *, int, double *, double *, double *, double *, double *, double *);
template void ApplyOperators<std::complex<double> >(Kpoint<std::complex<double>> *, int, std::complex<double> *, 
        std::complex<double> *, double *, double *,
        std::complex<double> *, std::complex<double> *);

// Applies A and B operators to one wavefunction
    template <typename KpointType>
void ApplyOperators (Kpoint<KpointType> *kptr, int istate, KpointType *a_psi, KpointType *b_psi, double *vtot, double *vxc_psi,  KpointType *nv, KpointType *Bns)
{
    // We want a clean exit if user terminates early
    CheckShutdown();

    BaseGrid *G = kptr->G;
    Lattice *L = &Rmg_L;
    State<KpointType> *sp = &kptr->Kstates[istate];
    KpointType *psi = kptr->Kstates[istate].psi;

    double vel = L->get_omega() / (G->get_NX_GRID(1) * G->get_NY_GRID(1) * G->get_NZ_GRID(1));

    int dimx = G->get_PX0_GRID(1);
    int dimy = G->get_PY0_GRID(1);
    int dimz = G->get_PZ0_GRID(1);
    int pbasis = dimx * dimy * dimz;
    int pbasis_noncoll = pbasis * ct.noncoll_factor;

    bool potential_acceleration = (ct.potential_acceleration_constant_step > 0.0);
    potential_acceleration = potential_acceleration & (ct.scf_steps > 0);


    // Apply A operator to psi
    ApplyAOperator (psi, a_psi);
    for(int idx = 0;idx < pbasis;idx++) b_psi[idx] = psi[idx];
    if(ct.noncoll)
    {
        ApplyAOperator (&psi[pbasis], &a_psi[pbasis]);
        for(int idx = 0;idx < pbasis;idx++) b_psi[idx+pbasis] = psi[idx+pbasis];
    }



    // if complex orbitals apply gradient to orbital and compute dot products
    std::complex<double> *kdr = NULL;
    if(typeid(KpointType) == typeid(std::complex<double>)) {

        kdr = new std::complex<double>[pbasis *ct.noncoll_factor]();
        std::complex<double> I_t(0.0, 1.0);
        KpointType *gx = new KpointType[pbasis];
        KpointType *gy = new KpointType[pbasis];
        KpointType *gz = new KpointType[pbasis];

        ApplyGradient (psi, gx, gy, gz, APP_CI_EIGHT, "Coarse");

        for(int idx = 0;idx < pbasis;idx++) {

            kdr[idx] = -I_t * (kptr->kp.kvec[0] * (std::complex<double>)gx[idx] +
                    kptr->kp.kvec[1] * (std::complex<double>)gy[idx] +
                    kptr->kp.kvec[2] * (std::complex<double>)gz[idx]);
        }

        if(ct.noncoll)
        {
            ApplyGradient (&psi[pbasis], gx, gy, gz, APP_CI_EIGHT, "Coarse");

            for(int idx = 0;idx < pbasis;idx++) {

                kdr[idx+pbasis] = -I_t * (kptr->kp.kvec[0] * (std::complex<double>)gx[idx] +
                        kptr->kp.kvec[1] * (std::complex<double>)gy[idx] +
                        kptr->kp.kvec[2] * (std::complex<double>)gz[idx]);
            }
        }

        delete [] gz;
        delete [] gy;
        delete [] gx;

    }


    // Generate 2*V*psi
    KpointType *sg_twovpsi_t = new KpointType[pbasis];
    double *veff= NULL;
    if(potential_acceleration) {
        int active_threads = ct.MG_THREADS_PER_NODE;
        if(ct.mpi_queue_mode && (active_threads > 1)) active_threads--;
        int offset = (sp->istate / kptr->dvh_skip) * pbasis;
        int my_pe_x, my_pe_y, my_pe_z;
        kptr->G->pe2xyz(pct.gridpe, &my_pe_x, &my_pe_y, &my_pe_z);
        int my_pe_offset = my_pe_x % pct.coalesce_factor;
        veff = &kptr->dvh[offset*pct.coalesce_factor + my_pe_offset*pbasis];
    }
    else {
        veff = vtot;
    }


    CPP_genvpsi (psi, sg_twovpsi_t, veff, (void *)kdr, kptr->kp.kmag, dimx, dimy, dimz);

        // For central FD B is just the identity
    for(int idx = 0; idx < pbasis; idx++) a_psi[idx] = sg_twovpsi_t[idx] - a_psi[idx];
   
      
      if(ct.noncoll)
    {
        std::complex<double> *kdr_down = kdr + pbasis;
        CPP_genvpsi (&psi[pbasis], sg_twovpsi_t, veff, (void *)kdr_down, kptr->kp.kmag, dimx, dimy, dimz);
        for(int idx = 0; idx < pbasis; idx++) a_psi[idx + pbasis] = sg_twovpsi_t[idx] - a_psi[idx + pbasis];

        double *vxc_x = &vxc_psi[pbasis];
        double *vxc_y = &vxc_psi[2*pbasis];
        double *vxc_z = &vxc_psi[3*pbasis];
        std::complex<double> *a_psi_C = (std::complex<double> *)a_psi;

        for(int idx = 0; idx < pbasis; idx++) 
        {
            a_psi_C[idx] += 2.0 * psi[idx] * vxc_z[idx];
            a_psi_C[idx] += 2.0 * psi[idx+pbasis] * std::complex<double>(vxc_x[idx], -vxc_y[idx]); 
            a_psi_C[idx + pbasis] += -2.0 * psi[idx + pbasis] * vxc_z[idx];
            a_psi_C[idx + pbasis] += 2.0 * psi[idx] * std::complex<double>(vxc_x[idx], vxc_y[idx]); 
        } 

    }
    // Add in non-local which has already had B applied in AppNls
    for(int idx = 0; idx < pbasis_noncoll; idx++) a_psi[idx] += TWO * nv[idx];
    // Scale correctly.
    for (int idx = 0; idx < pbasis_noncoll; idx++) a_psi[idx] = 0.5 * vel * a_psi[idx];

    // Add in already applied Bns to b_psi for US
    if(!ct.norm_conserving_pp) for(int idx = 0; idx < pbasis_noncoll; idx++) b_psi[idx] += Bns[idx];

    if(kdr) delete [] kdr;
    delete [] sg_twovpsi_t;

} // end ApplyOperators


