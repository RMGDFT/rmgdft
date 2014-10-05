#include "BaseThread.h"
#include "RmgTimer.h"
#include "rmg_error.h"
#include "GlobalSums.h"
#include <typeinfo>
#include <complex>

template void GlobalSums<float>(float*, int, MPI_Comm);
template void GlobalSums<double>(double*, int, MPI_Comm);
template void GlobalSums<std::complex<double> >(std::complex <double>*, int, MPI_Comm);

static double *fixed_vector1 = NULL;
static double *fixed_vector2 = NULL;
#define MAX_FIXED_VECTOR 2048




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
}


template <typename RmgType> void GlobalSums (RmgType * vect, int length, MPI_Comm comm)
{
    int tid;
    RmgTimer RT0("GlobalSums");
    BaseThread *T = BaseThread::getBaseThread(0);
    RmgType *v1, *v2;

    v1 = (RmgType *)fixed_vector1;
    v2 = (RmgType *)fixed_vector2;

    tid = T->get_thread_tid();
    if(tid < 0) {

        if(typeid(RmgType) == typeid(int))
            MPI_Allreduce(MPI_IN_PLACE, vect, length, MPI_INT, MPI_SUM, comm);

        if(typeid(RmgType) == typeid(float))
            MPI_Allreduce(MPI_IN_PLACE, vect, length, MPI_FLOAT, MPI_SUM, comm);

        if(typeid(RmgType) == typeid(double))
            MPI_Allreduce(MPI_IN_PLACE, vect, length, MPI_DOUBLE, MPI_SUM, comm);

    }
    else {

        if(length < MAX_FIXED_VECTOR) {

            for(int idx = 0;idx < length;idx++)
                v1[length * tid + idx] = vect[idx];
            T->thread_barrier_wait();

            if(tid == 0)
                MPI_Allreduce(v1, v2, length * T->get_threads_per_node(), MPI_DOUBLE, MPI_SUM, comm);

            T->thread_barrier_wait();
            for(int idx = 0;idx < length;idx++) vect[idx] = v2[length * tid + idx];

            return;
        }

        rmg_error_handler(__FILE__, __LINE__, "You have attempted to do a large threaded global sums. Please reconsider what you are trying to accomplish.\n");

    }


} // end GlobalSums

