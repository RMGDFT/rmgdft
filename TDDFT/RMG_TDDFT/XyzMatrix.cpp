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
#include "rmg_alloc.h"
#include "rmg_error.h"
#include "rmgthreads.h"
#include "RmgTimer.h"
#include "RmgThread.h"
#include "GlobalSums.h"
#include "Kpoint.h"
#include "RmgGemm.h"
#include "Subdiag.h"
#include "GpuAlloc.h"
#include "ErrorFuncs.h"
#include "blas.h"
#include "RmgParallelFft.h"

#include "common_prototypes.h"
#include "common_prototypes1.h"
#include "transition.h"
#include "prototypes_tddft.h"

#if CUDA_ENABLED
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cublas_v2.h>
#endif


template void XyzMatrix<double>(Kpoint<double> *, double *, int n, int m, int l);
template void XyzMatrix<std::complex<double> >(Kpoint<std::complex<double>> *, std::complex<double> *, int n, int m, int l);

template <typename KpointType>
void XyzMatrix (Kpoint<KpointType> *kptr, KpointType *Aij, int n, int m, int l)
{

    BaseGrid *G = kptr->G;
    Lattice *L = kptr->L;

    int num_states = kptr->nstates;
    int pbasis = kptr->pbasis;
    double vel = L->get_omega() / ((double)(G->get_NX_GRID(1) * G->get_NY_GRID(1) * G->get_NZ_GRID(1)));

    static KpointType *tmp_arrayT;

    int factor = 1;
    if(!ct.is_gamma) factor = 2;

    char *trans_t = "t";
    char *trans_n = "n";
    char *trans_c = "c";
    char *trans_a;
    if(typeid(KpointType) == typeid(std::complex<double>)) {
         trans_a = trans_c;
    }
    else {
        trans_a = trans_t;
    }   


    tmp_arrayT = (KpointType *)RmgMallocHost(pbasis*num_states * sizeof(KpointType));
    
    int ix, iy, iz;
    Rmg_G->pe2xyz (pct.gridpe, &ix, &iy, &iz);
    double hxgrid = Rmg_G->get_hxgrid(1);
    double hygrid = Rmg_G->get_hygrid(1);
    double hzgrid = Rmg_G->get_hzgrid(1);

    int px0_grid = Rmg_G->get_PX0_GRID(1);
    int py0_grid = Rmg_G->get_PY0_GRID(1);
    int pz0_grid = Rmg_G->get_PZ0_GRID(1);
    double xoff = ix * px0_grid * hxgrid;
    double yoff = iy * py0_grid * hygrid;
    double zoff = iz * pz0_grid * hzgrid;

    double xtal[3], xcrt[3];
    for(int i = 0; i < px0_grid; i++)
        for(int j = 0; j < py0_grid; j++)
            for(int k = 0; k < pz0_grid; k++)
            {

                xtal[0] = xoff + i * hxgrid;
                xtal[1] = yoff + j * hygrid;
                xtal[2] = zoff + k * hzgrid;
                Rmg_L.to_cartesian(xtal, xcrt);
                
                double xyz = std::pow(xcrt[0], n) * std::pow(xcrt[1], m) * std::pow(xcrt[2], l);
                
                int idx = i * py0_grid * pz0_grid + j * pz0_grid + k;
                for (int st1 = 0; st1 < num_states; st1++)
                {
                    tmp_arrayT[st1 * pbasis + idx] = kptr->Kstates[st1].psi[idx] * xyz;
                } 
            }

    /* tmp_arrayT:   V|psi> */

    // Compute A matrix
    KpointType alpha(vel);
    KpointType beta(0.0);
    RmgGemm(trans_a, trans_n, num_states, num_states, pbasis, alpha, kptr->orbital_storage, pbasis, tmp_arrayT, 
            pbasis, beta, Aij, num_states);

    MPI_Allreduce(MPI_IN_PLACE, (double *)Aij, num_states * num_states * factor, MPI_DOUBLE, MPI_SUM, pct.grid_comm);

    RmgFreeHost(tmp_arrayT);
}

