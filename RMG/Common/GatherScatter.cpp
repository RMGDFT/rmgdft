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

template void GatherPsi<double, float>(BaseGrid *, int, int, double *, float *);
template void GatherPsi<double, double>(BaseGrid *, int, int, double *, double *);
template void GatherPsi<std::complex<double>, std::complex<float> >(BaseGrid *, int, int, std::complex<double> *, std::complex<float> *);
template void GatherPsi<std::complex<double>, std::complex<double> >(BaseGrid *, int, int, std::complex<double> *, std::complex<double> *);
template void ScatterPsi<float, double>(BaseGrid *, int, int, float *, double *);
template void ScatterPsi<double, double>(BaseGrid *, int, int, double *, double *);
template void ScatterPsi<std::complex<double>, std::complex<double> >(BaseGrid *, int, int, std::complex<double> *, std::complex<double> *);
template void ScatterPsi<std::complex<float>, std::complex<double> >(BaseGrid *, int, int, std::complex<float> *, std::complex<double> *);
template void GatherGrid<double>(BaseGrid *, int, double *, double *);
template void GatherGrid<std::complex<double>>(BaseGrid *, int, std::complex<double> *, std::complex<double> *);
template void GatherEigs<double>(Kpoint<double> *);
template void GatherEigs<std::complex<double>> (Kpoint<std::complex<double>> *);


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
void GatherPsi(BaseGrid *G, int n, int istate, OrbitalType *A, CalcType *B)
{
    int chunksize = n / pct.coalesce_factor;
    int my_pe_x, pe_y, pe_z;
    G->pe2xyz(pct.gridpe, &my_pe_x, &pe_y, &pe_z);
    int pe_offset = my_pe_x % pct.coalesce_factor;
    CopyAndConvert(chunksize, &A[istate*chunksize], &B[pe_offset*chunksize]);
    if(!ct.coalesce_states) return;

    BaseThread *T = BaseThread::getBaseThread(0);
    int tid = T->get_thread_tid();
    int active_threads = ct.MG_THREADS_PER_NODE;
    if(ct.mpi_queue_mode) active_threads--;
    if(active_threads < 1) active_threads = 1;
    int base_istate = istate / (active_threads * pct.coalesce_factor);
    base_istate *= (active_threads * pct.coalesce_factor);
    CalcType *sbuf = new CalcType[n];

    std::atomic_bool is_completed_r[MAX_CFACTOR];
    std::atomic_bool is_completed_s[MAX_CFACTOR];
    std::atomic_int group_count{0};
    group_count.fetch_add(2*(pct.coalesce_factor-1), std::memory_order_seq_cst);

    mpi_queue_item_t qitems_r[MAX_CFACTOR];
    mpi_queue_item_t qitems_s[MAX_CFACTOR];

    for(int it = 0;it < MAX_CFACTOR;it++)
    {
        qitems_r[it].is_completed = &is_completed_r[it];
        qitems_s[it].is_completed = &is_completed_s[it];
        qitems_r[it].group_count = &group_count;
        qitems_s[it].group_count = &group_count;

    }


    // States are coalesced so we have to get the remote parts of istate
    for(int i=0;i < pct.coalesce_factor;i++)
    {
        // Queue receives
        if(i != pe_offset)
        {
            qitems_r[i].comm = T->get_unique_comm(istate);
            qitems_r[i].is_unpacked = false;
            qitems_r[i].is_completed->store(false);
            // The actual tag passed to mpi consists of the state index shifted left 5 bits and
            // the passed tag parameter. Most MPI implementations use 32 bits for the tag so this
            // gives 2^27 bits for the state index. Crays use 2^21 bits for the tag which gives
            // 2^16 states while the MPI spec only requires 2^16 bits which gives 2^11 states.
            qitems_r[i].mpi_tag = (istate<<5);
            qitems_r[i].target = Rmg_T->target_node[i - pe_offset + MAX_CFACTOR][1][1];
            qitems_r[i].type = RMG_MPI_IRECV;
            qitems_r[i].buf = (void *)&B[i*chunksize];
            qitems_r[i].buflen = sizeof(CalcType)*chunksize;

            // Push it onto the queue
            Rmg_Q->queue[tid]->push(qitems_r[i]);

        }
        else
        {
            qitems_r[i].is_unpacked = true;
            qitems_r[i].is_completed->store(true);
        } 
    }

    // Next we send the parts of states that other MPI procs require
    for(int i=0;i < pct.coalesce_factor;i++)
    {
        // Queue sends
        int remote_istate = base_istate + i * active_threads + istate % active_threads;
        if(istate != remote_istate)
        {
            qitems_s[i].comm = T->get_unique_comm(remote_istate);
            qitems_s[i].is_unpacked = false;
            qitems_s[i].is_completed->store(false);
            qitems_s[i].mpi_tag = (remote_istate<<5);
            qitems_s[i].target = Rmg_T->target_node[i - pe_offset + MAX_CFACTOR][1][1];
            qitems_s[i].type = RMG_MPI_ISEND;
#if GPU_ENABLED
            if(typeid(OrbitalType) == typeid(CalcType))
                cudaMemcpy(&sbuf[i*chunksize], &A[remote_istate*chunksize], chunksize*sizeof(OrbitalType), cudaMemcpyDefault);
            else
                CopyAndConvert(chunksize, &A[remote_istate*chunksize], &sbuf[i*chunksize]);
#else
            CopyAndConvert(chunksize, &A[remote_istate*chunksize], &sbuf[i*chunksize]);
#endif
//            qitems_s[i].buf = (void *)&A[remote_istate*chunksize];
            qitems_s[i].buf = (void *)&sbuf[i*chunksize];
            qitems_s[i].buflen = sizeof(CalcType)*chunksize;

            // Push it onto the queue
            Rmg_Q->queue[tid]->push(qitems_s[i]);
        } 
        else
        {
            qitems_s[i].is_unpacked = true;
            qitems_s[i].is_completed->store(true);
        } 
    }

    Rmg_Q->waitgroup(group_count);
    delete [] sbuf;

}

template <typename CalcType, typename OrbitalType>
void ScatterPsi(BaseGrid *G, int n, int istate, CalcType *A, OrbitalType *B)
{
    int chunksize = n / pct.coalesce_factor;
    int my_pe_x, pe_y, pe_z;
    G->pe2xyz(pct.gridpe, &my_pe_x, &pe_y, &pe_z);
    int pe_offset = my_pe_x % pct.coalesce_factor;
    CopyAndConvert(chunksize, &A[pe_offset*chunksize], &B[istate*chunksize]);
    if(!ct.coalesce_states) return;

    BaseThread *T = BaseThread::getBaseThread(0);
    int tid = T->get_thread_tid();
    int active_threads = ct.MG_THREADS_PER_NODE;
    if(ct.mpi_queue_mode) active_threads--;
    if(active_threads < 1) active_threads = 1;
    int base_istate = istate / (active_threads * pct.coalesce_factor);
    base_istate *= (active_threads * pct.coalesce_factor);
    CalcType *rbuf = new CalcType[n];


    std::atomic_bool is_completed_r[MAX_CFACTOR];
    std::atomic_bool is_completed_s[MAX_CFACTOR];
    std::atomic_int group_count{0};
    group_count.fetch_add(2*(pct.coalesce_factor-1), std::memory_order_seq_cst);

    mpi_queue_item_t qitems_r[MAX_CFACTOR];
    mpi_queue_item_t qitems_s[MAX_CFACTOR];

    for(int it = 0;it < MAX_CFACTOR;it++)
    {
        qitems_r[it].is_completed = &is_completed_r[it];
        qitems_s[it].is_completed = &is_completed_s[it];
        qitems_r[it].group_count = &group_count;
        qitems_s[it].group_count = &group_count;

    }

    // States are coalesced so we have to get the remote parts of istate
    for(int i=0;i < pct.coalesce_factor;i++)
    {
        // Queue receives
        int remote_istate = base_istate + i * active_threads + istate % active_threads;
        if(istate != remote_istate)
        {
            qitems_r[i].comm = T->get_unique_comm(remote_istate);
            qitems_r[i].is_unpacked = false;
            qitems_r[i].is_completed->store(false);
            qitems_r[i].mpi_tag = (remote_istate<<5);
            qitems_r[i].target = Rmg_T->target_node[i - pe_offset + MAX_CFACTOR][1][1];
            qitems_r[i].type = RMG_MPI_IRECV;
//            qitems_r[i].buf = (void *)&B[remote_istate*chunksize];
            qitems_r[i].buf = (void *)&rbuf[i*chunksize];
            qitems_r[i].buflen = sizeof(CalcType)*chunksize;

            // Push it onto the queue
            Rmg_Q->queue[tid]->push(qitems_r[i]);

        }
        else
        {
            qitems_r[i].is_unpacked = true;
            qitems_r[i].is_completed->store(true);
        }
    }

    // Next we send the parts of states that other MPI procs require
    for(int i=0;i < pct.coalesce_factor;i++)
    {
        // Queue sends
        if(i != pe_offset)
        {
            qitems_s[i].comm = T->get_unique_comm(istate);
            qitems_s[i].is_unpacked = false;
            qitems_s[i].is_completed->store(false);
            qitems_s[i].mpi_tag = (istate<<5);
            qitems_s[i].target = Rmg_T->target_node[i - pe_offset + MAX_CFACTOR][1][1];
            qitems_s[i].type = RMG_MPI_ISEND;
            qitems_s[i].buf = (void *)&A[i*chunksize];
            qitems_s[i].buflen = sizeof(CalcType)*chunksize;

            // Push it onto the queue
            Rmg_Q->queue[tid]->push(qitems_s[i]);
        }
        else
        {
            qitems_s[i].is_unpacked = true;
            qitems_s[i].is_completed->store(true);
        }
    }

    Rmg_Q->waitgroup(group_count);

    for(int i=0;i < pct.coalesce_factor;i++)
    {

        int remote_istate = base_istate + i * active_threads + istate % active_threads;
        if(istate != remote_istate)
        {

            CopyAndConvert(chunksize, &rbuf[i*chunksize], &B[remote_istate*chunksize]);

        }
    }

    delete [] rbuf;

}

// This is used to generate an array that represents a coalesced common domain. A good
// example being vtot_psi.
// n = non coalesced matrix size
// A = non coalesced matrix
// B = coalesced matrix (size = pct.coalesce_factor*n)
template <typename DataType>
void GatherGrid(BaseGrid *G, int n, DataType *A, DataType *B)
{
    if(!ct.coalesce_states || (pct.coalesce_factor == 1))
    {
        for(int i=0;i < n;i++) B[i] = A[i];
        return;
    }

    int my_pe_x, pe_y, pe_z;
    G->pe2xyz(pct.gridpe, &my_pe_x, &pe_y, &pe_z);
    int pe_offset = my_pe_x % pct.coalesce_factor;
    int length = pct.coalesce_factor * n;

    for(int i=0;i < length;i++) B[i] = 0.0;

    for(int i=0;i < n;i++) B[i + pe_offset*n] = A[i];

    if(typeid(DataType) == typeid(std::complex<double>)) length *= 2;
    MPI_Allreduce(MPI_IN_PLACE, B, length, MPI_DOUBLE, MPI_SUM, pct.coalesced_local_comm);

}

// Called after MgridSubspace has finished in order to sync eigenvalues to PE 0.
template <typename OrbitalType>
void GatherEigs(Kpoint<OrbitalType> *kptr)
{
    int my_pe_x, pe_y, pe_z;
    Rmg_G->pe2xyz(pct.gridpe, &my_pe_x, &pe_y, &pe_z);
    int pe_offset = my_pe_x % pct.coalesce_factor;
    int active_threads = ct.MG_THREADS_PER_NODE;
    if(ct.mpi_queue_mode) active_threads--;
    if(active_threads < 1) active_threads = 1;
    int istop = kptr->nstates / (active_threads * pct.coalesce_factor);
    istop = istop * active_threads * pct.coalesce_factor;
    int istart = pe_offset*active_threads;


    double *eigs = new double[kptr->nstates];
    double *neigs = new double[kptr->nstates];

    for(int i = 0;i < kptr->nstates;i++) eigs[i] = 0.0;

    for(int st1=0;st1 < istop;st1+=active_threads*pct.coalesce_factor) 
    {
        for(int ist = 0;ist < active_threads;ist++)
        {
            eigs[kptr->Kstates[st1 + ist + istart].istate] = kptr->Kstates[st1 + ist + istart].feig[0];
        }
    } 

    MPI_Reduce(eigs, neigs, kptr->nstates, MPI_DOUBLE, MPI_SUM , 0, pct.coalesced_local_comm);

    for(int st1=0;st1 < kptr->nstates;st1++) kptr->Kstates[st1].feig[0] = neigs[st1];

    delete [] neigs;
    delete [] eigs;
}
