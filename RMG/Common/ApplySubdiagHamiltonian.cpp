/*
 *
 * Copyright 2025 The RMG Project Developers. See the COPYRIGHT file 
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
#include <boost/pool/pool.hpp>
#include "TradeImages.h"
#include "FiniteDiff.h"
#include "Mgrid.h"
#include "RmgSumAll.h"
#include "BlasWrappers.h"
#include "const.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "common_prototypes.h"
#include "common_prototypes1.h"
#include "rmg_error.h"
#include "GlobalSums.h"
#include "Kpoint.h"
#include "packfuncs.h"
#include "transition.h"
#include "BaseThread.h"
#include "MpiQueue.h"
#include "GatherScatter.h"
#include "Solvers.h"
#include "rmg_complex.h"


template void ApplySubdiagHamiltonian<double,float>(Kpoint<double> *, State<double> *, double *, double *, double *, double *);
template void ApplySubdiagHamiltonian<double,double>(Kpoint<double> *, State<double> *, double *, double *, double *, double *);
template void ApplySubdiagHamiltonian<std::complex<double>, std::complex<float> >(Kpoint<std::complex<double>> *, State<std::complex<double> > *, double *,
        double *, std::complex<double> *, std::complex<double> *);
template void ApplySubdiagHamiltonian<std::complex<double>, std::complex<double> >(Kpoint<std::complex<double>> *, State<std::complex<double> > *, double *,
        double *, std::complex<double> *, std::complex<double> *);


template <typename OrbitalType, typename CalcType>
void ApplySubdiagHamiltonian (Kpoint<OrbitalType> *kptr, State<OrbitalType> * sp, double * vtot_psi, double *vxc_psi, OrbitalType *nv, OrbitalType *ns)
{
    BaseThread *Thread = BaseThread::getBaseThread(0);
    int tid = Thread->get_thread_tid();

    BaseGrid *G = kptr->G;
    int dimx = G->get_PX0_GRID(1) * pct.coalesce_factor;
    int dimy = G->get_PY0_GRID(1);
    int dimz = G->get_PZ0_GRID(1);

    int pbasis = dimx * dimy * dimz;
    int pbasis_noncoll = pbasis * ct.noncoll_factor;

    size_t aratio = sizeof(OrbitalType) / sizeof(CalcType);

    // Per thread memory pool in Kpoint class. If you need more than 32 chunks modify
    // pool creation code in Kpoint.cpp. Done like this to enable changing the type
    // of memory allocated in one place.
    int pool_blocks = 0;
    boost::pool<rmg_user_allocator> *p = kptr->kalloc[tid];
    CalcType *res2_t = (CalcType *)p->ordered_malloc(1);pool_blocks++;
    CalcType *work1_t = (CalcType *)p->ordered_malloc(4);pool_blocks+=4;
    CalcType *tmp_psi_t  = (CalcType *)p->ordered_malloc(1);pool_blocks++;
    CalcType *hr0_t  =  (CalcType *)p->ordered_malloc(4);pool_blocks+=4;
    OrbitalType *nv_t  = (OrbitalType *)p->ordered_malloc(aratio);pool_blocks+=aratio;
    // Copy double precision psi into correct precison array
    GatherPsi(G, pbasis_noncoll, sp->istate, kptr->orbital_storage, tmp_psi_t, pct.coalesce_factor);

    // Copy nv into local array
    if(ct.coalesce_states)
        GatherPsi(G, pbasis_noncoll, sp->istate, nv, nv_t, pct.coalesce_factor);
    else
        GatherPsi(G, pbasis_noncoll, 0, nv, nv_t, 1);

    // For USPP copy double precision ns into correct precision temp array. For NCPP ns=psi. */
    if(ct.norm_conserving_pp)
    {
        for(int ix=0;ix < pbasis_noncoll;ix++) work1_t[ix] = tmp_psi_t[ix];
    }
    else
    {
        GatherPsi(G, pbasis_noncoll, sp->istate, kptr->ns, work1_t, pct.coalesce_factor);
    }

    ApplyHamiltonian<OrbitalType,CalcType> (
          kptr, sp,
          sp->istate,
          tmp_psi_t, // coalesced psi
          hr0_t,   // hamiltonian applied to  |kro>
          vtot_psi,
          vxc_psi,
          nv_t,
          false);

    // Copy single precision orbital back to double precision
    ScatterPsi(G, pbasis_noncoll, sp->istate, hr0_t, 
               &kptr->orbital_storage[kptr->nstates*kptr->pbasis], pct.coalesce_factor);

    p->free(res2_t, pool_blocks);

} // end ApplySubdiagHamiltonian


