#ifndef RMG_GlobalSums_h
#define RMG_GlobalSums_h

#include <mpi.h>

template <typename RmgType>
void GlobalSums (RmgType * vect, int length, MPI_Comm comm);
void GlobalSumsInit(void);


#endif

