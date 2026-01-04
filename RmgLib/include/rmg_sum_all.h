#pragma once

#ifndef RMG_rmg_sum_all_h
#define RMG_rmg_sum_all_h

#include <mpi.h>
namespace rmg
{
    template <typename T>
    T sum_all (T x, MPI_Comm comm);
    template <typename T>
    T max_all (T x, MPI_Comm comm);
} // end namespace rmg

#endif
