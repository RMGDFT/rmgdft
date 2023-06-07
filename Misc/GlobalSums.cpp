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

#include "BaseThread.h"
#include "RmgTimer.h"
#include "rmg_error.h"
#include "main.h"
#include "transition.h"
#include "GlobalSums.h"
#include "transition.h"
#include <typeinfo>
#include <complex>

template void GlobalSums<float>(float*, int, MPI_Comm);
template void GlobalSums<double>(double*, int, MPI_Comm);
template void GlobalSums<std::complex<double> >(std::complex <double>*, int, MPI_Comm);

static double *fixed_vector1 = NULL;
static double *fixed_vector2 = NULL;
#define MAX_FIXED_VECTOR 2048


static MPI_Comm *coalesced_comm_pool;
static MPI_Comm *coalesced_local_comm_pool;

size_t block_size = 67103864;

void GlobalSumsInit(void) {
    int retval;
    BaseThread *T = BaseThread::getBaseThread(0);

    retval = MPI_Alloc_mem(2 * sizeof(double) * T->get_threads_per_node() * MAX_FIXED_VECTOR , MPI_INFO_NULL, &fixed_vector1);
    if(retval != MPI_SUCCESS) {
        rmg_error_handler(__FILE__, __LINE__, "Error in MPI_Alloc_mem.\n");
    }
    retval = MPI_Alloc_mem(2 * sizeof(double) * T->get_threads_per_node() * MAX_FIXED_VECTOR , MPI_INFO_NULL, &fixed_vector2);
    if(retval != MPI_SUCCESS) {
        rmg_error_handler(__FILE__, __LINE__, "Error in MPI_Alloc_mem.\n");
    }

    coalesced_comm_pool = new MPI_Comm[20*ct.MG_THREADS_PER_NODE + 1];
    coalesced_local_comm_pool = new MPI_Comm[20*ct.MG_THREADS_PER_NODE + 1];
    for(int thread = 0;thread < 20*ct.MG_THREADS_PER_NODE + 1;thread++)
    {
    }

    if(ct.coalesce_states)
    {
        for(int thread = 0;thread < 20*ct.MG_THREADS_PER_NODE + 1;thread++)
        {
            MPI_Comm_dup(pct.coalesced_local_comm, &coalesced_local_comm_pool[thread]);
        MPI_Comm_dup(pct.coalesced_grid_comm, &coalesced_comm_pool[thread]);
        }
    }
}

MPI_Comm get_unique_coalesced_comm(int istate)
{
   int comm_index = istate % (20*ct.MG_THREADS_PER_NODE + 1);
   return coalesced_comm_pool[comm_index];
}

MPI_Comm get_unique_coalesced_local_comm(int istate)
{
   int comm_index = istate % (20*ct.MG_THREADS_PER_NODE + 1);
   return coalesced_local_comm_pool[comm_index];
}

template <typename RmgType> void GlobalSums (RmgType * vect, int length, MPI_Comm comm)
{
    RmgTimer RT0("GlobalSums");
    BaseThread *T = BaseThread::getBaseThread(0);

    RmgType *v1, *v2;

    v1 = (RmgType *)fixed_vector1;
    v2 = (RmgType *)fixed_vector2;

    if(ct.mpi_queue_mode && T->is_loop_over_states())
    {
        int tid = T->get_thread_tid();
        if(tid < 0) tid = 0;
        int istate = T->get_thread_basetag();
        std::atomic_bool is_completed;
        is_completed.store(false);

        mpi_queue_item_t qi;
        qi.is_completed = &is_completed;
        qi.type = RMG_MPI_SUM;
        if(ct.coalesce_states && (pct.coalesce_factor > 1))
        {
            int comm_index = istate % (20*ct.MG_THREADS_PER_NODE + 1);
            qi.comm = coalesced_comm_pool[comm_index];
        }
        else
        {
            qi.comm = T->get_unique_comm(istate);
        }

        qi.buflen = length;
        qi.buf = (void *)vect;

        if(typeid(RmgType) == typeid(int)) qi.datatype = MPI_INT;
        if(typeid(RmgType) == typeid(float)) qi.datatype = MPI_FLOAT;
        if(typeid(RmgType) == typeid(double)) qi.datatype = MPI_DOUBLE;
//        std::atomic_thread_fence(std::memory_order_seq_cst);
        //Rmg_Q->push(tid, qi);
        Rmg_Q->queue[tid]->push(qi);
        while(!is_completed.load(std::memory_order_acquire)){;}
//        std::atomic_thread_fence(std::memory_order_seq_cst);
        return;
    }

    int tid = T->get_thread_tid();

    if(tid < 0) {

        if(typeid(RmgType) == typeid(int))
            MPI_Allreduce(MPI_IN_PLACE, vect, length, MPI_INT, MPI_SUM, comm);

        if(typeid(RmgType) == typeid(float))
            MPI_Allreduce(MPI_IN_PLACE, vect, length, MPI_FLOAT, MPI_SUM, comm);

        if(typeid(RmgType) == typeid(double))
            MPI_Allreduce(MPI_IN_PLACE, vect, length, MPI_DOUBLE, MPI_SUM, comm);

        if(typeid(RmgType) == typeid(std::complex<double>))
            MPI_Allreduce(MPI_IN_PLACE, vect, length, MPI_DOUBLE_COMPLEX, MPI_SUM, comm);

    }
    else {

        if(length < MAX_FIXED_VECTOR) {

            for(int idx = 0;idx < length;idx++)
                v1[length * tid + idx] = vect[idx];
            T->thread_barrier_wait(false);

            if(tid == 0)
                MPI_Allreduce(v1, v2, length * T->get_threads_per_node(), MPI_DOUBLE, MPI_SUM, comm);

            T->thread_barrier_wait(false);
            for(int idx = 0;idx < length;idx++) vect[idx] = v2[length * tid + idx];

            return;
        }

        rmg_error_handler(__FILE__, __LINE__, "You have attempted to do a large threaded global sums. Please reconsider what you are trying to accomplish.\n");

    }


} // end GlobalSums


// Used by subdiag routines to get around integer size limitations for large Allreduce operations.
// Should only be called from a non-threaded region.
void BlockAllreduce(double *mat, size_t count, MPI_Comm comm)
{
    BaseThread *T = BaseThread::getBaseThread(0);

    if(T->is_loop_over_states())
        rmg_error_handler(__FILE__, __LINE__, "BlockAllReduce cannot be called from a threaded region.\n");

    size_t blocks = count / block_size;
    size_t rem = count % block_size;
    double *tptr = mat;
    for(size_t ib=0;ib < blocks;ib++)
    {
        MPI_Allreduce(MPI_IN_PLACE, tptr, block_size, MPI_DOUBLE, MPI_SUM, comm);
        tptr += block_size;
    }
    if(rem)
        MPI_Allreduce(MPI_IN_PLACE, tptr, rem, MPI_DOUBLE, MPI_SUM, comm);
}

void BlockAllreduce(float *mat, size_t count, MPI_Comm comm)
{
    BaseThread *T = BaseThread::getBaseThread(0);

    if(T->is_loop_over_states())
        rmg_error_handler(__FILE__, __LINE__, "BlockAllReduce cannot be called from a threaded region.\n");

    size_t blocks = count / block_size;
    size_t rem = count % block_size;
    float *tptr = mat;
    for(size_t ib=0;ib < blocks;ib++)
    {
        MPI_Allreduce(MPI_IN_PLACE, tptr, block_size, MPI_FLOAT, MPI_SUM, comm);
        tptr += block_size;
    }
    if(rem)
        MPI_Allreduce(MPI_IN_PLACE, tptr, rem, MPI_FLOAT, MPI_SUM, comm);
}

void BlockAllreduce(std::complex<double> *mat, size_t count, MPI_Comm comm)
{
    BaseThread *T = BaseThread::getBaseThread(0);

    if(T->is_loop_over_states())
        rmg_error_handler(__FILE__, __LINE__, "BlockAllReduce cannot be called from a threaded region.\n");

    size_t blocks = (2 * count) / block_size;
    size_t rem = (2 * count) % block_size;
    double *tptr = (double *)mat;
    for(size_t ib=0;ib < blocks;ib++)
    {
        MPI_Allreduce(MPI_IN_PLACE, tptr, block_size, MPI_DOUBLE, MPI_SUM, comm);
        tptr += block_size;
    }
    if(rem)
        MPI_Allreduce(MPI_IN_PLACE, tptr, rem, MPI_DOUBLE, MPI_SUM, comm);
}

void BlockAllreduce(std::complex<float> *mat, size_t count, MPI_Comm comm)
{
    BaseThread *T = BaseThread::getBaseThread(0);

    if(T->is_loop_over_states())
        rmg_error_handler(__FILE__, __LINE__, "BlockAllReduce cannot be called from a threaded region.\n");

    size_t blocks = (2 * count) / block_size;
    size_t rem = (2 * count) % block_size;
    float *tptr = (float *)mat;
    for(size_t ib=0;ib < blocks;ib++)
    {
        MPI_Allreduce(MPI_IN_PLACE, tptr, block_size, MPI_FLOAT, MPI_SUM, comm);
        tptr += block_size;
    }
    if(rem)
        MPI_Allreduce(MPI_IN_PLACE, tptr, rem, MPI_FLOAT, MPI_SUM, comm);
}


