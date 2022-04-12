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

    int active_threads = ct.MG_THREADS_PER_NODE;
    if(ct.mpi_queue_mode) active_threads--;
    if(active_threads < 1) active_threads = 1;

    int istop = num_states / active_threads;
    istop = istop * active_threads;

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

void PreconditionerOne (double *res, int st, double gamma)
{

    BaseGrid *G = Rmg_G;
    TradeImages *T =Rmg_T;
    Lattice *L = &Rmg_L;
    Mgrid MG(L, T);
    int pre[MAX_MG_LEVELS] = { 2, 2, 2, 2, 20, 20, 20, 20 };
    int post[MAX_MG_LEVELS] = { 2, 2, 2, 2, 2, 2, 2, 2 };
    int levels = ct.eig_parm.levels;
    double Zfac = 2.0 * ct.max_zvalence;
    double tstep = 0.666666666666;

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

    double *res_t = new double[dimx * dimy * dimz];
    double *res_t2 = new double[dimx * dimy * dimz];

    double beta = 0.5;
    int nits = ct.eig_parm.gl_pre + ct.eig_parm.gl_pst;

    double one = 1.0;
    int ione = 1;
    for(int idx = 0; idx < pbasis; idx++)
        res_t[idx] = gamma * res[idx + st * pbasis];

    /* Smoothing cycles */
    for (int cycles = 0; cycles <= nits; cycles++)
    {

        RmgTimer *RT = new RmgTimer("Precond: app_cil");
        double diag = CPP_app_cil_driver (&Rmg_L, Rmg_T, res_t, res_t2, dimx, dimy, dimz,
                hxgrid, hygrid, hzgrid, APP_CI_FOURTH);
        delete RT;
        daxpy(&pbasis, &one, &res[st * pbasis], &ione, res_t2, &ione);
        LocalOrbital->ApplyBoundary(res_t2, st);
        //for(int idx = 0; idx < pbasis; idx++) if (!LocalOrbital->mask[st * pbasis + idx])
        //    res_t2[idx] = 0.0;

        double t1; 
        /* Now either smooth the wavefunction or do a multigrid cycle */
        if (cycles == ct.eig_parm.gl_pre)
        {
            //CPP_pack_ptos_convert ((float *)work1_t, (double *)res_t2, dimx, dimy, dimz);
            CPP_pack_ptos (work1_t, res_t2, dimx, dimy, dimz);
            //MG.mgrid_solv<float>((float *)work2_t, (float *)work1_t, (float *)work_t,
            RT= new RmgTimer("Precond: mgrid");
            MG.mgrid_solv<double>(work2_t, work1_t, work_t,
                    dimx, dimy, dimz, hxgrid, hygrid, hzgrid,
                    0, levels, pre, post, 1,
                    tstep, 1.0*Zfac, 0.1, NULL,     // which one is best?
                    //tstep, 1.0, 0.0, vtot,
                    G->get_NX_GRID(1), G->get_NY_GRID(1), G->get_NZ_GRID(1),
                    G->get_PX_OFFSET(1), G->get_PY_OFFSET(1), G->get_PZ_OFFSET(1),
                    G->get_PX0_GRID(1), G->get_PY0_GRID(1), G->get_PZ0_GRID(1), ct.boundaryflag);
            delete RT;
            //CPP_pack_stop_convert((float *)work2_t, (double *)res_t2, dimx, dimy, dimz);
            CPP_pack_stop(work2_t, res_t2, dimx, dimy, dimz);

            t1 = -1.;

        }
        else
        {
            double t5 = diag - Zfac;
            t5 = -1.0 / t5;
            t1 = ct.eig_parm.gl_step * t5;

        }                       /* end if cycles == ct.eig_parm.gl_pre */

        daxpy(&pbasis, &t1, res_t2, &ione, res_t, &ione);

        LocalOrbital->ApplyBoundary(res_t, st);
//        for(int idx = 0; idx < pbasis; idx++) if (!LocalOrbital->mask[st * pbasis + idx])
 //           res_t[idx] = 0.0;

    }

    for(int idx = 0; idx < pbasis; idx++) res[idx + st * pbasis] = beta * res_t[idx];

    delete []work_t;
    delete []res_t;
    delete []res_t2;
}


