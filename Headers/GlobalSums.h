#ifndef RMG_GlobalSums_h
#define RMG_GlobalSums_h

#include <mpi.h>

template <typename RmgType>
void GlobalSums (RmgType * vect, int length, MPI_Comm comm);
void GlobalSumsInit(void);
MPI_Comm get_unique_coalesced_comm(int istate);
MPI_Comm get_unique_coalesced_local_comm(int istate);


#endif

