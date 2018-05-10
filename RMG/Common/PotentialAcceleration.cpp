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
#include <atomic>
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
#include "rmgthreads.h"
#include <boost/atomic.hpp>


std::mutex vtot_sync_mutex;
std::mutex dvtot_sync_mutex;

template void PotentialAcceleration<double,float>(Kpoint<double> *, State<double> *, double *, double *, float *, double *);
template void PotentialAcceleration<double,double>(Kpoint<double> *, State<double> *, double *, double *, double *, double *);
template void PotentialAcceleration<std::complex<double>, std::complex<float> >(Kpoint<std::complex<double> > *, State<std::complex<double> > *, double *, double *, std::complex<float> *, std::complex<double> *);
template void PotentialAcceleration<std::complex<double>, std::complex<double> >(Kpoint<std::complex<double> > *, State<std::complex<double> > *, double *, double *, std::complex<double> *, std::complex<double> *);


std::atomic_int counter{0};
std::atomic_int skipper{0};

void PotentialAccelerationReset(int skip)
{
    counter.store(skip);
}

void PotentialAccelerationWait(int istate)
{
//if(pct.gridpe==0)printf("HHH0 %d  %d  %d  %d\n",pct.gridpe,istate,counter.load(),start);
    while(istate >= counter.load()){;}
//if(pct.gridpe==0)printf("HHH1 %d  %d  %d  %d\n",pct.gridpe,istate,counter.load(),start);
}

template <typename OrbitalType, typename CalcType>
void PotentialAcceleration(Kpoint<OrbitalType> *kptr, State<OrbitalType> *sp, double *vtot_psi, double *dvtot_psi, CalcType *tmp_psi_t, OrbitalType *saved_psi)
{
    BaseThread *T = BaseThread::getBaseThread(0);
    int tid = T->get_thread_tid();
    SCF_THREAD_CONTROL *s = (SCF_THREAD_CONTROL *)T->get_pptr(tid);
    int active_threads = 1;
    if(ct.MG_THREADS_PER_NODE > 1) active_threads = s->extratag1;
    int base_state = s->extratag2;
    int pbasis = kptr->pbasis * pct.coalesce_factor;
    int skip = kptr->dvh_skip;
    if(ct.coalesce_states) skip = active_threads * pct.coalesce_factor;
    int offset = (sp->istate / skip) * pbasis;

    double t1 = 1.8 * ct.potential_acceleration_constant_step;
    if(sp->occupation[0] < 0.5) t1 = 0.0;

    vtot_sync_mutex.lock();
    for(int idx = 0;idx <pbasis;idx++) {
       kptr->dvh[idx + offset] += t1 * PI * sp->occupation[0] * std::real(tmp_psi_t[idx] * std::conj((tmp_psi_t[idx] - (CalcType)saved_psi[idx])));
    }
    vtot_sync_mutex.unlock();

    if(sp->istate == (counter.load() - 1))
    {
        while(skipper.load() != skip/pct.coalesce_factor-1) {;}
        if(pct.coalesce_factor > 1)
        {
//          GlobalSums (&kptr->dvh[offset], pbasis, pct.coalesced_local_comm);
            MPI_Allreduce(MPI_IN_PLACE, &kptr->dvh[offset], pbasis, MPI_DOUBLE, MPI_SUM, pct.coalesced_local_comm);
        }
        for(int i = 0;i <pbasis;i++)
        {
            kptr->dvh[i + offset] += vtot_psi[i];
            vtot_psi[i] = kptr->dvh[i + offset];
        }
        skipper.store(0);

        // Adjust for end point boundary when coalescing is true
        if(pct.coalesce_factor > 1)
        {
            int tthreads = active_threads;
            while(tthreads > 1)
            {
                int icheck = base_state + 2*tthreads*pct.coalesce_factor;
                if(icheck > kptr->nstates)
                {
                    skip -= pct.coalesce_factor;
                    tthreads--;
                }
                else
                {
                    break;
                }
            }
        }
        counter.fetch_add(skip);
    }
    else
    {
        skipper++;
    }

}

