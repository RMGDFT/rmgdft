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
#include "pfft.h"

#include "transition.h"
#include "blas.h"

#if GPU_ENABLED
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cublas_v2.h>
#endif



template void DavPreconditioner<double>(Kpoint<double> *, double *, double *, double, double *, double *, int, double);
template void DavPreconditioner<std::complex<double> >(Kpoint<std::complex<double>> *, std::complex<double> *, std::complex<double> *, 
                               double, double *, double *, int, double);


template <typename OrbitalType>
void DavPreconditioner (Kpoint<OrbitalType> *kptr, OrbitalType *psi, OrbitalType *res, double fd_diag, double *eigs, 
                        double *vtot, int notconv, double avg_potential)
{
    BaseGrid *G = kptr->G;
    TradeImages *T = kptr->T;
    Lattice *L = kptr->L;
    Mgrid MG(L, T);
    int pre[8] = { 2, 2, 4, 20, 20, 20, 20, 20 };
    int post[8] = { 2, 2, 2, 2, 2, 2, 2, 2 };
    int levels = ct.eig_parm.levels;
    double Zfac = 2.0 * ct.max_zvalence;
    double tstep = 0.8;

    int ixoff, iyoff, izoff;
    int dimx = G->get_PX0_GRID(1);
    int dimy = G->get_PY0_GRID(1);
    int dimz = G->get_PZ0_GRID(1);
    int dx2 = MG.MG_SIZE (dimx, 0, G->get_NX_GRID(1), G->get_PX_OFFSET(1), G->get_PX0_GRID(1), &ixoff, ct.boundaryflag);
    int dy2 = MG.MG_SIZE (dimy, 0, G->get_NY_GRID(1), G->get_PY_OFFSET(1), G->get_PY0_GRID(1), &iyoff, ct.boundaryflag);
    int dz2 = MG.MG_SIZE (dimz, 0, G->get_NZ_GRID(1), G->get_PZ_OFFSET(1), G->get_PZ0_GRID(1), &izoff, ct.boundaryflag);

    OrbitalType *s_diag = kptr->s_diag;
    OrbitalType *vnl_diag = kptr->vnl_diag;

    int pbasis = kptr->pbasis;

    double hxgrid = G->get_hxgrid(1);
    double hygrid = G->get_hygrid(1);
    double hzgrid = G->get_hzgrid(1);

    OrbitalType *work_t = new OrbitalType[4*(dimx + 2)*(dimy + 2)*(dimz + 2)];
    OrbitalType *work1_t = new OrbitalType[4*(dimx + 2)*(dimy + 2)*(dimz + 2)];
    OrbitalType *work2_t = new OrbitalType[4*(dimx + 2)*(dimy + 2)*(dimz + 2)];
    double *nvtot = new double[4*(dimx + 2)*(dimy + 2)*(dimz + 2)];

    for(int idx = 0;idx <pbasis;idx++) nvtot[idx] = -vtot[idx];

    // Apply preconditioner
    for(int st1=0;st1 < notconv;st1++) {

        double t1 = 0.0;
        for(int idx = 0;idx <pbasis;idx++) t1 += std::real(res[st1*pbasis + idx]);

        t1 = real_sum_all(t1, pct.grid_comm) / (double)(G->get_NX_GRID(1) * G->get_NY_GRID(1) * G->get_NZ_GRID(1));
        // neutralize cell
        for(int idx = 0;idx <pbasis;idx++) res[st1*pbasis + idx] -= OrbitalType(t1);

        CPP_pack_ptos (work1_t, &res[st1*pbasis], dimx, dimy, dimz);
        MG.mgrid_solv (work2_t, work1_t, work_t,
                    dimx, dimy, dimz, hxgrid, hygrid, hzgrid,
                    0, G->get_neighbors(), levels, pre, post, 1,
                    tstep, 1.0*Zfac, -avg_potential, NULL,     // which one is best?
                    //tstep, 1.0*Zfac, 0.0, nvtot,
                    G->get_NX_GRID(1), G->get_NY_GRID(1), G->get_NZ_GRID(1),
                    G->get_PX_OFFSET(1), G->get_PY_OFFSET(1), G->get_PZ_OFFSET(1),
                    G->get_PX0_GRID(1), G->get_PY0_GRID(1), G->get_PZ0_GRID(1), ct.boundaryflag);
        CPP_pack_stop (work2_t, &res[st1*pbasis], dimx, dimy, dimz);

        for(int idx = 0;idx <pbasis;idx++) res[st1*pbasis + idx] += eigs[st1] * t1;;

#if 0

        double eps1 = 0.0;
        double eps2 = 0.0;
        for(int idx = 0;idx < pbasis;idx++) {
            double a1 = std::real(fd_diag + vtot[idx] + vnl_diag[idx] - eigs[st1]*s_diag[idx]);
            a1 = 1.0 / a1;
            double r = std::real(res[st1*pbasis + idx]);
            double x = std::real(psi[st1*pbasis + idx]);
            eps1 += r*x*a1;
            eps2 += x*x*a1;
        }
        double eps = eps1 / eps2;

        for(int idx = 0;idx < pbasis;idx++) {

            double a1 = std::real(fd_diag + vtot[idx] + vnl_diag[idx] - eigs[st1]*s_diag[idx]);
            a1 = 1.0 / a1;
            double r = std::real(res[st1*pbasis + idx]);
            double x = std::real(psi[st1*pbasis + idx]);
            //res[st1*pbasis + idx] = a1*(-r + eps*x);
            //psi[st1*pbasis + idx] = a1*(-r);
        }

 CPP_app_smooth_test (&res[st1*pbasis], work_t, dimx, dimy, dimz);
 for(int idx=0;idx < pbasis;idx++)res[st1*pbasis + idx] = work_t[idx];

        CPP_pack_ptos (work_t, &res[st1*pbasis], dimx, dimy, dimz);
        T->trade_images (work_t, dimx, dimy, dimz, FULL_TRADE);
        CPP_app_smooth (work_t, work1_t, dimx, dimy, dimz);
        T->trade_images (work1_t, dimx, dimy, dimz, FULL_TRADE);
        CPP_app_smooth (work1_t, work_t, dimx, dimy, dimz);
        CPP_pack_stop (work_t, &res[st1*pbasis], dimx, dimy, dimz);

#endif
    }

    delete [] nvtot;
    delete [] work2_t;
    delete [] work1_t;
    delete [] work_t;

}
