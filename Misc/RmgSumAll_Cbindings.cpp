
#include "BaseThread.h"
#include "RmgSumAll.h"
#include <typeinfo>


extern "C" int int_sum_all(int x, MPI_Comm comm)
{
    return RmgSumAll<int>(x, comm);
}

extern "C" double double_sum_all(double x, MPI_Comm comm)
{
    return RmgSumAll<double>(x, comm);
}

extern "C" float float_sum_all(float x, MPI_Comm comm)
{
    return RmgSumAll<float>(x, comm);
}

extern "C" double real_sum_all(double x, MPI_Comm comm)
{
    return RmgSumAll<double>(x, comm);
}
