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
#include "rmg_reduce.h"
#include "rmg_sum_all.h"
#include "Kpoint.h"
#include "rmg_gemm.h"
#include "rmg_mgrid.h"
#include "RmgException.h"
#include "Subdiag.h"
#include "Solvers.h"
#include "GpuAlloc.h"

#include "RmgParallelFft.h"
#include "TradeImages.h"
#include "packfuncs.h"

#include "transition.h"
#include "blas.h"

#if CUDA_ENABLED
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cublas_v2.h>
#endif

#include "init_var.h"
#include "prototypes_on.h"

void Preconditioner (double *res, int num_states)
{

    double gamma = get_gamma_precond(vtot_c, states[0].eig[0]);
    BaseThread *T = BaseThread::getBaseThread(0);

    RmgTimer RT("Precond");

    int active_threads = rmg::get_active_threads();
    int istop = num_states / active_threads;
    istop = istop * active_threads;
    if(active_threads == 1) istop = 0;

    Rmg_T->set_coalesce_factor(1);
    for(int st1=0;st1 < istop;st1+=active_threads) {

          SCF_THREAD_CONTROL thread_control;

          for(int ist = 0;ist < active_threads;ist++) {
              thread_control.job = HYBRID_ON_PRECOND;
              thread_control.p1 = (void *)res;
              thread_control.istate = st1 + ist;
              thread_control.basetag = st1 + ist;
              thread_control.gamma = gamma;
              QueueThreadTask(ist, thread_control);
          }

          // Thread tasks are set up so wake them
          if(!ct.mpi_queue_mode) T->run_thread_tasks(active_threads);

    }

    if(ct.mpi_queue_mode) T->run_thread_tasks(active_threads, Rmg_Q);

    // Process any remaining states in serial fashion
    if(istop < num_states)
    {
        for(int st = istop;st < num_states;st++) {
            PreconditionerOne(res, st, gamma);
        }
    }

}

void make_masks(rmg::mgrid &MG, TradeImages *T, rmg::grid *G, int dimx, int dimy, int dimz, int maxlevels, int st)
{
    int dx1, dy1, dz1; 
    int dx2, dy2, dz2; 
    int ixoff, iyoff, izoff;

    MG.boundary_masks[0].resize((dimx+2)*(dimy+2)*(dimz+2));
    rmg::pack_ptos (MG.boundary_masks[0].data(), LocalOrbital->boundary[st].data(), dimx, dimy, dimz);
    T->trade_images (MG.boundary_masks[0].data(), dimx, dimy, dimz, FULL_TRADE);

    dx1 = dimx;
    dy1 = dimy;
    dz1 = dimz;

    for(int i=1;i <= maxlevels;i++)
    {
        dx2 = MG.MG_SIZE (dx1, i-1, G->get_NX_GRID(1), G->get_PX_OFFSET(1), dimx, &ixoff, ct.boundaryflag);
        dy2 = MG.MG_SIZE (dy1, i-1, G->get_NY_GRID(1), G->get_PY_OFFSET(1), dimy, &iyoff, ct.boundaryflag);
        dz2 = MG.MG_SIZE (dz1, i-1, G->get_NZ_GRID(1), G->get_PZ_OFFSET(1), dimz, &izoff, ct.boundaryflag);
        MG.boundary_masks[i].resize((dx2+2)*(dy2+2)*(dz2+2));
        MG.mg_restrict (MG.boundary_masks[i-1].data(), MG.boundary_masks[i].data(), dx1, dy1, dz1, dx2, dy2, dz2, ixoff, iyoff, izoff);
        T->trade_images (MG.boundary_masks[i].data(), dx2, dy2, dz2, FULL_TRADE);
        dx1 = dx2;
        dy1 = dy2;
        dz1 = dz2;
    }
}

void PreconditionerOne (double *res, int st, double gamma)
{

    rmg::grid *G = Rmg_G;
    TradeImages *T =Rmg_T;
    Lattice *L = &Rmg_L;
    rmg::mgrid MG(L, T, Rmg_G, 1, 0.0);
    FiniteDiff FD(L);
    MG.on_flag = true;

    int levels = ct.eig_parm.levels;
    int dimx = G->get_PX0_GRID(1);
    int dimy = G->get_PY0_GRID(1);
    int dimz = G->get_PZ0_GRID(1);
    int pbasis = G->get_P0_BASIS(1);

    double hxgrid = G->get_hxgrid(1);
    double hygrid = G->get_hygrid(1);
    double hzgrid = G->get_hzgrid(1);

    double *work_t = new double[8*(dimx + 2)*(dimy + 2)*(dimz + 2)];
    double *work1_t = &work_t[4*(dimx + 2)*(dimy + 2)*(dimz + 2)];
    double *work2_t = &work_t[6*(dimx + 2)*(dimy + 2)*(dimz + 2)];
    double *pot = new double[pct.coalesce_factor*4*(dimx + 2)*(dimy + 2)*(dimz + 2)];

    double *res_t = new double[dimx * dimy * dimz];
    double *res_t2 = new double[dimx * dimy * dimz];
    double beta = 0.5;

    double one = 1.0;
    int ione = 1;
    for(int idx = 0; idx < pbasis; idx++)
        res_t[idx] = gamma * res[idx + st * pbasis];

    // Make multigrid boundary masks
    make_masks(MG, T, G, dimx, dimy, dimz, levels, st);

    /* Pre smoothing cycles */
    for (int cycles = 0; cycles <= ct.eig_parm.gl_pre; cycles++)
    {
        RmgTimer *RT = new RmgTimer("Precond: app_cil");
        double diag = CPP_app_cil_driver (&Rmg_L, Rmg_T, res_t, res_t2, dimx, dimy, dimz,
                hxgrid, hygrid, hzgrid, APP_CI_FOURTH);

        delete RT;
        daxpy(&pbasis, &one, &res[st * pbasis], &ione, res_t2, &ione);
        LocalOrbital->ApplyBoundary(res_t2, st);
        // We don't update with the residual here if it's time for a multigrid iteration
        if(cycles == ct.eig_parm.gl_pre) break;
        double t1 = -ct.eig_parm.gl_step / diag;
        daxpy(&pbasis, &t1, res_t2, &ione, res_t, &ione);
        LocalOrbital->ApplyBoundary(res_t, st);
    }

    // Multigrid V cycle
    if (levels)
    {
        rmg::pack_ptos (work1_t, res_t2, dimx, dimy, dimz);
        rmg::pack_ptos (pot, LocalOrbital->pot_precond[st].data(), dimx, dimy, dimz);
        RmgTimer *RT= new RmgTimer("Precond: mgrid");
        MG.mgrid_solv<double>(work2_t, work1_t, work_t,
                dimx, dimy, dimz,
                0, levels, 1.0, pot,
                G->get_PX0_GRID(1), G->get_PY0_GRID(1), G->get_PZ0_GRID(1));

        delete RT;
        rmg::pack_stop(work2_t, res_t2, dimx, dimy, dimz);

        double t1 = -1.;
        daxpy(&pbasis, &t1, res_t2, &ione, res_t, &ione);
        LocalOrbital->ApplyBoundary(res_t, st);

    }

    /* Post smoothing cycles */
    for (int cycles = 0; cycles < ct.eig_parm.gl_pst; cycles++)
    {
        RmgTimer *RT = new RmgTimer("Precond: app_cil");
        double diag = CPP_app_cil_driver (&Rmg_L, Rmg_T, res_t, res_t2, dimx, dimy, dimz,
                hxgrid, hygrid, hzgrid, APP_CI_FOURTH);
        delete RT;
        daxpy(&pbasis, &one, &res[st * pbasis], &ione, res_t2, &ione);
        LocalOrbital->ApplyBoundary(res_t2, st);
        double t1 = -ct.eig_parm.gl_step / diag;
        daxpy(&pbasis, &t1, res_t2, &ione, res_t, &ione);
        LocalOrbital->ApplyBoundary(res_t, st);
    }



    for(int idx = 0; idx < pbasis; idx++) res[idx + st * pbasis] = beta * res_t[idx];

    delete []work_t;
    delete []res_t;
    delete []res_t2;
    delete [] pot;
}


