
#include "BaseThread.h"
#include "rmg_sum_all.h"
#include <typeinfo>


int int_sum_all(int x, MPI_Comm comm)
{
    return rmg::sum_all<int>(x, comm);
}

double double_sum_all(double x, MPI_Comm comm)
{
    return rmg::sum_all<double>(x, comm);
}

float float_sum_all(float x, MPI_Comm comm)
{
    return rmg::sum_all<float>(x, comm);
}

double real_sum_all(double x, MPI_Comm comm)
{
    return rmg::sum_all<double>(x, comm);
}
