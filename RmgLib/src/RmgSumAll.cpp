/*
 *
 * Copyright (c) 1995, Emil Briggs
 * Copyright (C) 1998  Emil Briggs, Charles Brabec, Mark Wensell, 
 *                     Dan Sullivan, Chris Rapcewicz, Jerzy Bernholc
 * Copyright (C) 2001  Emil Briggs, Wenchang Lu,
 *                     Marco Buongiorno Nardelli,Charles Brabec, 
 *                     Mark Wensell,Dan Sullivan, Chris Rapcewicz,
 *                     Jerzy Bernholc
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the <organization> nor the
 *       names of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.

 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * 
*/

#include "BaseThread.h"
#include "RmgSumAll.h"
#include <typeinfo>


std::mutex RmgSumAllLock;
volatile int vector_state = 0;
double RmgSumAllVector[2*MAX_RMG_THREADS];
double RecvBuf[2*MAX_RMG_THREADS];

template int RmgSumAll<int>(int, MPI_Comm);
template float RmgSumAll<float>(float, MPI_Comm);
template double RmgSumAll<double>(double, MPI_Comm);

template <typename RmgType>
RmgType RmgSumAll (RmgType x, MPI_Comm comm)
{
    BaseThread *T = BaseThread::getBaseThread(0);
    int tid;
    RmgType inreg;
    RmgType outreg;
    RmgType *inbuf = (RmgType *)RmgSumAllVector;
    RmgType *outbuf = (RmgType *)RecvBuf;

    tid = T->get_thread_tid();
    if(tid < 0) {

        inreg = x;
        
        if(typeid(RmgType) == typeid(int))
            MPI_Allreduce (&inreg, &outreg, 1, MPI_INT, MPI_SUM, comm);

        if(typeid(RmgType) == typeid(float))
            MPI_Allreduce (&inreg, &outreg, 1, MPI_FLOAT, MPI_SUM, comm);

        if(typeid(RmgType) == typeid(double))
            MPI_Allreduce (&inreg, &outreg, 1, MPI_DOUBLE, MPI_SUM, comm);

        return outreg;

    }

    RmgSumAllLock.lock();
        vector_state = 1;
        inbuf[tid] = x;
    RmgSumAllLock.unlock();

    // Wait until everyone gets here
    T->thread_barrier_wait(true);

    RmgSumAllLock.lock();
        if(vector_state == 1) {

            if(typeid(RmgType) == typeid(int))
                MPI_Allreduce(inbuf, outbuf, T->get_threads_per_node(), MPI_INT, MPI_SUM, comm);

            if(typeid(RmgType) == typeid(float))
                MPI_Allreduce(inbuf, outbuf, T->get_threads_per_node(), MPI_FLOAT, MPI_SUM, comm);

            if(typeid(RmgType) == typeid(double))
                MPI_Allreduce(inbuf, outbuf, T->get_threads_per_node(), MPI_DOUBLE, MPI_SUM, comm);

            vector_state = 0;
        }
    RmgSumAllLock.unlock();

    return outbuf[tid]; 
}

