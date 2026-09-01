#ifndef RMG_GlobalSums_h
#define RMG_GlobalSums_h

#include <mpi.h>

template <typename RmgType>
void GlobalSums (RmgType * vect, int length, MPI_Comm comm);
void GlobalSumsInit(void);
MPI_Comm get_unique_coalesced_comm(int istate);
MPI_Comm get_unique_coalesced_local_comm(int istate);

// Used by subdiag routines to get around integer size limitations for large Allreduce operations.
// Should only be called from a non-threaded region.
void BlockAllreduce(double *mat, size_t count, MPI_Comm comm);
void BlockAllreduce(float *mat, size_t count, MPI_Comm comm);
#ifdef __cplusplus
#include <complex>
void BlockAllreduce(std::complex<double> *mat, size_t count, MPI_Comm comm);
void BlockAllreduce(std::complex<float> *mat, size_t count, MPI_Comm comm);
#endif


#endif

