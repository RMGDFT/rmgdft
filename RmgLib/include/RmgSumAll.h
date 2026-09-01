#ifndef RMG_RmgSumAll_h
#define RMG_RmgSumAll_h

#include <mpi.h>

template <typename RmgType>
RmgType RmgSumAll (RmgType x, MPI_Comm comm);
template <typename RmgType>
RmgType RmgMaxAll (RmgType x, MPI_Comm comm);


#endif
