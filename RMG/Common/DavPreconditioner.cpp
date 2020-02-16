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
#include "GatherScatter.h"
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

template void DavPreconditionerOne<double>(Kpoint<double> *, int, double *, double, double, double *, double);
template void DavPreconditionerOne<std::complex<double> >(Kpoint<std::complex<double>> *, int, std::complex<double> *, 
                               double, double, double *, double);

template <typename OrbitalType>
void DavPreconditioner (Kpoint<OrbitalType> *kptr, OrbitalType *res, double fd_diag, double *eigs, 
                        double *vtot, int notconv, double avg_potential)
{

    BaseThread *T = BaseThread::getBaseThread(0);

    int active_threads = ct.MG_THREADS_PER_NODE;
    if(ct.mpi_queue_mode) active_threads--;
    if(active_threads < 1) active_threads = 1;

    int my_pe_x, my_pe_y, my_pe_z;
    kptr->G->pe2xyz(pct.gridpe, &my_pe_x, &my_pe_y, &my_pe_z);
    int my_pe_offset = my_pe_x % pct.coalesce_factor;


    int dimx = kptr->G->get_PX0_GRID(1);
    int dimy = kptr->G->get_PY0_GRID(1);
    int dimz = kptr->G->get_PZ0_GRID(1);

    double *tvtot = vtot;
    if(pct.coalesce_factor > 1)
    {
        tvtot = new double[dimx*dimy*dimz*pct.coalesce_factor];
        GatherGrid(kptr->G, kptr->pbasis, vtot, tvtot);
    }

    double *nvtot = new double[pct.coalesce_factor*4*(dimx + 2)*(dimy + 2)*(dimz + 2)];
    CPP_pack_ptos (nvtot, tvtot, dimx*pct.coalesce_factor, dimy, dimz);
    for(int idx = 0;idx <(pct.coalesce_factor*dimx + 2)*(dimy + 2)*(dimz + 2);idx++) nvtot[idx] = -nvtot[idx];

    int istop = notconv / (active_threads * pct.coalesce_factor);
    istop = istop * (active_threads * pct.coalesce_factor);
    kptr->T->set_coalesce_factor(pct.coalesce_factor);
    int istart = my_pe_offset*active_threads;

    for(int st1=0;st1 < istop;st1+=active_threads*pct.coalesce_factor) {

          SCF_THREAD_CONTROL thread_control;

          for(int ist = 0;ist < active_threads;ist++) {
              thread_control.job = HYBRID_DAV_PRECONDITIONER;
              thread_control.vtot = nvtot;
              thread_control.p1 = (void *)kptr;
              thread_control.p2 = (void *)res;
              thread_control.avg_potential = avg_potential;
              thread_control.fd_diag = fd_diag;
              thread_control.eig = eigs[st1 + ist + istart];
              thread_control.basetag = st1 + ist + istart;
              thread_control.extratag1 = active_threads;
              thread_control.extratag2 = st1;
              QueueThreadTask(ist, thread_control);
          }

          // Thread tasks are set up so wake them
          if(!ct.mpi_queue_mode) T->run_thread_tasks(active_threads);

    }

    if(ct.mpi_queue_mode) T->run_thread_tasks(active_threads, Rmg_Q);

    kptr->T->set_coalesce_factor(1);
    // Process any remaining states in serial fashion
    for(int st1 = istop;st1 < notconv;st1++) {
//        DavPreconditionerOne (kptr, st1, res, fd_diag, eigs[st1], vtot, avg_potential);
    }

    if(pct.coalesce_factor > 1) delete [] tvtot;
    delete [] nvtot;
}

template <typename OrbitalType>
void DavPreconditionerOne (Kpoint<OrbitalType> *kptr, int st, OrbitalType *res, double fd_diag, double eig, 
                        double *vtot, double avg_potential)
{
    // We want a clean exit if user terminates early
    CheckShutdown();

    BaseGrid *G = kptr->G;
    TradeImages *T = kptr->T;
    Lattice *L = kptr->L;
    Mgrid MG(L, T);
    int pre[MAX_MG_LEVELS] = { 2, 8, 8, 20, 20, 20, 20, 20 };
    int post[MAX_MG_LEVELS] = { 2, 2, 2, 2, 2, 2, 2, 2 };
    int levels = ct.eig_parm.levels;
    double Zfac = 2.0 * ct.max_zvalence;
    double tstep = 0.666666666666;

    int coalesce_factor = T->get_coalesce_factor();
    int dimx = G->get_PX0_GRID(1) * coalesce_factor;
    int dimy = G->get_PY0_GRID(1);
    int dimz = G->get_PZ0_GRID(1);
    int pbasis = dimx * dimy * dimz;
    int pbasis_noncoll = pbasis;

    double hxgrid = G->get_hxgrid(1);
    double hygrid = G->get_hygrid(1);
    double hzgrid = G->get_hzgrid(1);

    OrbitalType *work_t = (OrbitalType *)malloc(8*(dimx + 2)*(dimy + 2)*(dimz + 2) * sizeof(OrbitalType));
    OrbitalType *work1_t = &work_t[4*(dimx + 2)*(dimy + 2)*(dimz + 2)];
    OrbitalType *work2_t = &work_t[6*(dimx + 2)*(dimy + 2)*(dimz + 2)];

//printf("DDDDDD  %d\n",st);fflush(NULL);
    // Apply preconditioner
    for(int is=0;is < ct.noncoll_factor;is++)
    {
        // Typedefs to map different data types to correct MG template.
        typedef typename std::conditional_t< std::is_same<OrbitalType, double>::value, float,
                         std::conditional_t< std::is_same<OrbitalType, std::complex<double>>::value, std::complex<float>,
                         std::conditional_t< std::is_same<OrbitalType, std::complex<float>>::value, std::complex<float>, float> > > mgtype_t;
        typedef typename std::conditional_t< std::is_same<OrbitalType, double>::value, double,
                         std::conditional_t< std::is_same<OrbitalType, std::complex<double>>::value, std::complex<double>,
                         std::conditional_t< std::is_same<OrbitalType, std::complex<float>>::value, std::complex<float>, float> > > convert_type_t;

        GatherPsi(G, pbasis, st, res, work2_t, coalesce_factor);
        double t1 = 0.0;
        for(int idx = 0;idx <pbasis;idx++) t1 += std::real(work2_t[idx]);

        if(coalesce_factor==1)
        {
            GlobalSums (&t1, 1, pct.grid_comm);
        }
        else
        {
            GlobalSums (&t1, 1, pct.coalesced_grid_comm);
        }

        t1 /= (double)(G->get_NX_GRID(1) * G->get_NY_GRID(1) * G->get_NZ_GRID(1));

        // neutralize cell
        for(int idx = 0;idx <pbasis;idx++) work2_t[idx] -= OrbitalType(t1);

        CPP_pack_ptos_convert ((mgtype_t *)work1_t, (convert_type_t *)work2_t, dimx, dimy, dimz);
        MG.mgrid_solv<mgtype_t>((mgtype_t *)work2_t, (mgtype_t *)work1_t, (mgtype_t *)work_t,
                    dimx, dimy, dimz, hxgrid, hygrid, hzgrid,
                    0, levels, pre, post, 1,
                    //tstep, 1.0*Zfac, -avg_potential, NULL,     // which one is best?
                    tstep, 1.0, 0.0, vtot,
                    G->get_NX_GRID(1), G->get_NY_GRID(1), G->get_NZ_GRID(1),
                    G->get_PX_OFFSET(1), G->get_PY_OFFSET(1), G->get_PZ_OFFSET(1),
                    coalesce_factor*G->get_PX0_GRID(1), G->get_PY0_GRID(1), G->get_PZ0_GRID(1), ct.boundaryflag);
        CPP_pack_stop_convert((mgtype_t *)work2_t, (convert_type_t *)work1_t, dimx, dimy, dimz);

        for(int idx = 0;idx <pbasis;idx++) work1_t[idx] += eig * t1;;
        ScatterPsi(G, pbasis, st, work1_t, res, coalesce_factor);
    }
    free(work_t);
}
