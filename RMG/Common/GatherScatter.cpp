/*
 *
 * Copyright 2018 The RMG Project Developers. See the COPYRIGHT file 
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
#include "const.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "common_prototypes.h"
#include "common_prototypes1.h"
#include "rmg_error.h"
#include "RmgTimer.h"
#include "transition.h"
#include "BaseThread.h"
#include "GatherScatter.h"
#include "MpiQueue.h"

template void GatherGrid<double, float>(BaseGrid *, int, int, double *, float *);
template void GatherGrid<double, double>(BaseGrid *, int, int, double *, double *);
template void GatherGrid<std::complex<double>, std::complex<float> >(BaseGrid *, int, int, std::complex<double> *, std::complex<float> *);
template void GatherGrid<std::complex<double>, std::complex<double> >(BaseGrid *, int, int, std::complex<double> *, std::complex<double> *);
template void ScatterGrid<double, float>(BaseGrid *, int, int, double *, float *);
template void ScatterGrid<double, double>(BaseGrid *, int, int, double *, double *);
template void ScatterGrid<float, double>(BaseGrid *, int, int, float *, double *);
template void ScatterGrid<std::complex<double>, std::complex<float> >(BaseGrid *, int, int, std::complex<double> *, std::complex<float> *);
template void ScatterGrid<std::complex<double>, std::complex<double> >(BaseGrid *, int, int, std::complex<double> *, std::complex<double> *);
template void ScatterGrid<std::complex<float>, std::complex<double> >(BaseGrid *, int, int, std::complex<float> *, std::complex<double> *);

void CopyAndConvert(int n, double *A, float *B)
{
    for(int idx = 0;idx < n;idx++)
        B[idx] = (float)A[idx];
}

void CopyAndConvert(int n, double *A, double *B)
{
    for(int idx = 0;idx < n;idx++)
        B[idx] = A[idx];
}

void CopyAndConvert(int n, std::complex<double> *A, std::complex<float> *B)
{
    for(int idx = 0;idx < n;idx++)
        B[idx] = (std::complex<float>)A[idx];
}
void CopyAndConvert(int n, std::complex<double> *A, std::complex<double> *B)
{
    for(int idx = 0;idx < n;idx++)
        B[idx] = (std::complex<double>)A[idx];
}

void CopyAndConvert(int n, float *A, double *B)
{
    for(int idx = 0;idx < n;idx++)
        B[idx] = (double)A[idx];
}

void CopyAndConvert(int n, std::complex<float> *A, std::complex<double> *B)
{
    for(int idx = 0;idx < n;idx++)
        B[idx] = (std::complex<double>)A[idx];
}




template <typename OrbitalType, typename CalcType>
void GatherGrid(BaseGrid *G, int n, int istate, OrbitalType *A, CalcType *B)
{
    if(!ct.coalesce_states)
    {
        CopyAndConvert(n, A, B);
        return;
    }

    BaseThread *T = BaseThread::getBaseThread(0);
    int tid = T->get_thread_tid();
    if(tid < 0) tid = 0;
    int sblock = tid * pct.coalesce_factor;

    int chunksize = n / pct.coalesce_factor;
    int offset = istate % pct.coalesce_factor;
    int my_pe_x, target_pe_x, pe_y, pe_z;
    G->pe2xyz(pct.gridpe, &my_pe_x, &pe_y, &pe_z);
    int base_pe_x = (my_pe_x / pct.coalesce_factor) * pct.coalesce_factor;

    mpi_queue_item_t *qitems_r = new mpi_queue_item_t[pct.coalesce_factor];
    mpi_queue_item_t *qitems_s = new mpi_queue_item_t[pct.coalesce_factor];


    if(typeid(OrbitalType) == typeid(CalcType))
    {
        // States are coalesced so we have to get the remote parts of istate
        for(int i=0;i < pct.coalesce_factor;i++)
        {
            // Queue receives
            if(i != offset)
            {
            } 
        }
        // Next we send the parts of states that other MPI procs require
        for(int st=istate;st < istate + pct.coalesce_factor;st++)
        {
            // Queue sends
            if(st != istate)
            {
            } 
        }

        // Now do the local part
        CopyAndConvert(chunksize, A, &B[offset*chunksize]);
    }

    delete [] qitems_s;
    delete [] qitems_r;

}
template <typename OrbitalType, typename CalcType>
void ScatterGrid(BaseGrid *G, int n, int istate, OrbitalType *A, CalcType *B)
{
    if(!ct.coalesce_states)
    {
        CopyAndConvert(n, A, B);
        return;
    }
}

