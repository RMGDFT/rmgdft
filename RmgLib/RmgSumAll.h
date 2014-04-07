#ifndef RMG_BlasWrappers_h
#define RMG_BlasWrappers_h

#include <mpi.h>

template <typename RmgType>
RmgType RmgSumAll (RmgType x, MPI_Comm comm);


#endif
