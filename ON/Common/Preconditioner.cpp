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



void Preconditioner (double *res, int num_states)
{

    BaseGrid *G = Rmg_G;
    TradeImages *T =Rmg_T;
    Lattice *L = &Rmg_L;
    Mgrid MG(L, T);
    int pre[MAX_MG_LEVELS] = { 2, 8, 8, 20, 20, 20, 20, 20 };
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

    double *work_t = new double[8*(dimx + 2)*(dimy + 2)*(dimz + 2) * sizeof(double)];
    double *work1_t = &work_t[4*(dimx + 2)*(dimy + 2)*(dimz + 2)];
    double *work2_t = &work_t[6*(dimx + 2)*(dimy + 2)*(dimz + 2)];


    for(int st = 0 ; st < num_states; st++)
    {
        CPP_pack_ptos (work1_t, &res[st*pbasis], dimx, dimy, dimz);
        MG.mgrid_solv (work2_t, work1_t, work_t,
                dimx, dimy, dimz, hxgrid, hygrid, hzgrid,
                0, levels, pre, post, 1,
                tstep, 1.0*Zfac, 0.0, NULL,     // which one is best?
                //tstep, 1.0, 0.0, vtot,
                G->get_NX_GRID(1), G->get_NY_GRID(1), G->get_NZ_GRID(1),
                G->get_PX_OFFSET(1), G->get_PY_OFFSET(1), G->get_PZ_OFFSET(1),
                G->get_PX0_GRID(1), G->get_PY0_GRID(1), G->get_PZ0_GRID(1), ct.boundaryflag);
        CPP_pack_stop (work2_t, &res[st*pbasis], dimx, dimy, dimz);
    }


    free(work_t);
}
