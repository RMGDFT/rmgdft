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
#include <omp.h>
#include <cmath>
#include <float.h>
#include "FiniteDiff.h"
#include "const.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "rmgthreads.h"
#include "RmgTimer.h"
#include "RmgThread.h"
#include "GlobalSums.h"
#include "RmgSumAll.h"
#include "Kpoint.h"
#include "RmgGemm.h"
#include "Mgrid.h"
#include "RmgException.h"
#include "Subdiag.h"
#include "Solvers.h"
#include "GpuAlloc.h"
#include "ErrorFuncs.h"
#include "RmgParallelFft.h"
#include "TradeImages.h"
#include "packfuncs.h"

#include "transition.h"
#include "blas.h"

#if GPU_ENABLED
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cublas_v2.h>
#endif



template void DavPreconditioner<double>(Kpoint<double> *, double *, double, double *, double *, int, double);
template void DavPreconditioner<std::complex<double> >(Kpoint<std::complex<double>> *, std::complex<double> *, 
                               double, double *, double *, int, double);

template void DavPreconditionerOne<double>(Kpoint<double> *, double *, double, double, double *, double);
template void DavPreconditionerOne<std::complex<double> >(Kpoint<std::complex<double>> *, std::complex<double> *, 
                               double, double, double *, double);

template <typename OrbitalType>
void DavPreconditioner (Kpoint<OrbitalType> *kptr, OrbitalType *res, double fd_diag, double *eigs, 
                        double *vtot, int notconv, double avg_potential)
{

    BaseThread *T = BaseThread::getBaseThread(0);

    int active_threads = ct.MG_THREADS_PER_NODE;
    if(ct.mpi_queue_mode) active_threads--;
    if(active_threads < 1) active_threads = 1;

    int dimx = kptr->G->get_PX0_GRID(1);
    int dimy = kptr->G->get_PY0_GRID(1);
    int dimz = kptr->G->get_PZ0_GRID(1);

    double *nvtot = new double[2*(dimx + 2)*(dimy + 2)*(dimz + 2)];
    CPP_pack_ptos (nvtot, vtot, dimx, dimy, dimz);
    for(int idx = 0;idx <(dimx + 2)*(dimy + 2)*(dimz + 2);idx++) nvtot[idx] = -nvtot[idx];

    int istop = notconv / active_threads;
    istop = istop * active_threads;

    for(int st1=0;st1 < istop;st1+=active_threads) {

          SCF_THREAD_CONTROL thread_control;

          for(int ist = 0;ist < active_threads;ist++) {
              thread_control.job = HYBRID_DAV_PRECONDITIONER;
              thread_control.vtot = nvtot;
              thread_control.p1 = (void *)kptr;
              thread_control.p2 = (void *)&res[(st1 + ist) * kptr->pbasis];
              thread_control.avg_potential = avg_potential;
              thread_control.fd_diag = fd_diag;
              thread_control.eig = eigs[st1 + ist];
              thread_control.basetag = st1 + ist;
              QueueThreadTask(ist, thread_control);
          }

          // Thread tasks are set up so wake them
          if(!ct.mpi_queue_mode) T->run_thread_tasks(active_threads);

    }

    if(ct.mpi_queue_mode) T->run_thread_tasks(active_threads, Rmg_Q);

    // Process any remaining states in serial fashion
    for(int st1 = istop;st1 < notconv;st1++) {
        DavPreconditionerOne (kptr, &res[st1 * kptr->pbasis], fd_diag, eigs[st1], nvtot, avg_potential);
    }

    delete [] nvtot;
}

template <typename OrbitalType>
void DavPreconditionerOne (Kpoint<OrbitalType> *kptr, OrbitalType *res, double fd_diag, double eig, 
                        double *vtot, double avg_potential)
{
    BaseGrid *G = kptr->G;
    TradeImages *T = kptr->T;
    Lattice *L = kptr->L;
    Mgrid MG(L, T);
    int pre[MAX_MG_LEVELS] = { 2, 8, 8, 20, 20, 20, 20, 20 };
    int post[MAX_MG_LEVELS] = { 2, 2, 2, 2, 2, 2, 2, 2 };
    int levels = ct.eig_parm.levels;
    double Zfac = 2.0 * ct.max_zvalence;
    double tstep = 0.666666666666;

    int dimx = G->get_PX0_GRID(1);
    int dimy = G->get_PY0_GRID(1);
    int dimz = G->get_PZ0_GRID(1);
    int pbasis = kptr->pbasis;

    double hxgrid = G->get_hxgrid(1);
    double hygrid = G->get_hygrid(1);
    double hzgrid = G->get_hzgrid(1);

    OrbitalType *work_t = (OrbitalType *)malloc(8*(dimx + 2)*(dimy + 2)*(dimz + 2) * sizeof(OrbitalType));
    OrbitalType *work1_t = &work_t[2*(dimx + 2)*(dimy + 2)*(dimz + 2)];
    OrbitalType *work2_t = &work_t[4*(dimx + 2)*(dimy + 2)*(dimz + 2)];


    // Apply preconditioner
    double t1 = 0.0;
    for(int idx = 0;idx <pbasis;idx++) t1 += std::real(res[idx]);

    GlobalSums (&t1, 1, pct.grid_comm);
    t1 /= (double)(G->get_NX_GRID(1) * G->get_NY_GRID(1) * G->get_NZ_GRID(1));

    // neutralize cell
    for(int idx = 0;idx <pbasis;idx++) res[idx] -= OrbitalType(t1);
    if(typeid(OrbitalType) == typeid(double))
    {
        CPP_pack_ptos_convert ((float *)work1_t, (double *)res, dimx, dimy, dimz);
        MG.mgrid_solv<float>((float *)work2_t, (float *)work1_t, (float *)work_t,
                    dimx, dimy, dimz, hxgrid, hygrid, hzgrid,
                    0, G->get_neighbors(), levels, pre, post, 1,
                    //tstep, 1.0*Zfac, -avg_potential, NULL,     // which one is best?
                    tstep, 1.0, 0.0, vtot,
                    G->get_NX_GRID(1), G->get_NY_GRID(1), G->get_NZ_GRID(1),
                    G->get_PX_OFFSET(1), G->get_PY_OFFSET(1), G->get_PZ_OFFSET(1),
                    G->get_PX0_GRID(1), G->get_PY0_GRID(1), G->get_PZ0_GRID(1), ct.boundaryflag);
        CPP_pack_stop_convert((float *)work2_t, (double *)res, dimx, dimy, dimz);
    }
    else
    {
        CPP_pack_ptos (work1_t, res, dimx, dimy, dimz);
        MG.mgrid_solv (work2_t, work1_t, work_t,
                    dimx, dimy, dimz, hxgrid, hygrid, hzgrid,
                    0, G->get_neighbors(), levels, pre, post, 1,
                    tstep, 1.0*Zfac, -avg_potential, NULL,     // which one is best?
                    //tstep, 1.0*Zfac, 0.0, nvtot,
                    G->get_NX_GRID(1), G->get_NY_GRID(1), G->get_NZ_GRID(1),
                    G->get_PX_OFFSET(1), G->get_PY_OFFSET(1), G->get_PZ_OFFSET(1),
                    G->get_PX0_GRID(1), G->get_PY0_GRID(1), G->get_PZ0_GRID(1), ct.boundaryflag);
        CPP_pack_stop (work2_t, res, dimx, dimy, dimz);
    }

    for(int idx = 0;idx <pbasis;idx++) res[idx] += eig * t1;;

    free(work_t);
}
