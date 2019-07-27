
#include "BaseThread.h"
#include "RmgSumAll.h"
#include <typeinfo>


int int_sum_all(int x, MPI_Comm comm)
{
    return RmgSumAll<int>(x, comm);
}

double double_sum_all(double x, MPI_Comm comm)
{
    return RmgSumAll<double>(x, comm);
}

float float_sum_all(float x, MPI_Comm comm)
{
    return RmgSumAll<float>(x, comm);
}

double real_sum_all(double x, MPI_Comm comm)
{
    return RmgSumAll<double>(x, comm);
}
