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
#include "TradeImages.h"
#include "FiniteDiff.h"
#include "Mgrid.h"
#include "BlasWrappers.h"
#include "const.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "common_prototypes.h"
#include "common_prototypes1.h"
#include "rmg_alloc.h"
#include "rmg_error.h"
#include "RmgTimer.h"
#include "GlobalSums.h"
#include "Kpoint.h"
#include "packfuncs.h"
#include "transition.h"


std::mutex vtot_sync_mutex;

template void PotentialAcceleration<double,float>(Kpoint<double> *, State<double> *, double *, double *, float *, double *);
template void PotentialAcceleration<double,double>(Kpoint<double> *, State<double> *, double *, double *, double *, double *);
template void PotentialAcceleration<std::complex<double>, std::complex<float> >(Kpoint<std::complex<double> > *, State<std::complex<double> > *, double *, double *, std::complex<float> *, std::complex<double> *);
template void PotentialAcceleration<std::complex<double>, std::complex<double> >(Kpoint<std::complex<double> > *, State<std::complex<double> > *, double *, double *, std::complex<double> *, std::complex<double> *);



template <typename OrbitalType, typename CalcType>
void PotentialAcceleration(Kpoint<OrbitalType> *kptr, State<OrbitalType> *sp, double *vtot_psi, double *nvtot_psi, CalcType *tmp_psi_t, OrbitalType *saved_psi)
{

    int pbasis = kptr->pbasis;


    // Save approximate potential used for this orbital and update potential for future orbitals
    if(!(sp->istate % kptr->dvh_skip)) {
        int offset = (sp->istate / kptr->dvh_skip) * pbasis;
        for(int i = 0;i <pbasis;i++) {
            kptr->dvh[i + offset] = vtot_psi[i];
        }
    }


    double t1 = 1.8 * ct.potential_acceleration_constant_step;
    if(sp->occupation[0] < 0.5) t1 = 0.0;

    vtot_sync_mutex.lock();
    for(int idx = 0;idx <pbasis;idx++) {
       vtot_psi[idx] = vtot_psi[idx] + t1 * PI * sp->occupation[0] * std::real(tmp_psi_t[idx] * std::conj((tmp_psi_t[idx] - (CalcType)saved_psi[idx])));
    }
    vtot_sync_mutex.unlock();

}

