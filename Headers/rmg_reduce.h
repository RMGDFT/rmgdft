#pragma once
#ifndef RMG_rmg_reduce_h
#define RMG_rmg_reduce_h

#include <mpi.h>
#include <complex>

MPI_Comm get_unique_coalesced_comm(int istate);
MPI_Comm get_unique_coalesced_local_comm(int istate);

// Used by subdiag routines to get around integer size limitations for large Allreduce operations.
// Should only be called from a non-threaded region.
namespace rmg {
    void init_reduce(void);
    template <typename RmgType>
    void allreduce (RmgType * vect, int length, MPI_Comm comm);
    void block_allreduce(double *mat, size_t count, MPI_Comm comm);
    void block_allreduce(float *mat, size_t count, MPI_Comm comm);
    void block_allreduce(std::complex<double> *mat, size_t count, MPI_Comm comm);
    void block_allreduce(std::complex<float> *mat, size_t count, MPI_Comm comm);
    template <typename RmgType>
    void reduce (RmgType * vect, int length, int root, MPI_Comm comm);
    void block_reduce(double *mat, size_t count, int root, MPI_Comm comm);
    void block_reduce(float *mat, size_t count, int root, MPI_Comm comm);
    void block_reduce(std::complex<double> *mat, size_t count, int root, MPI_Comm comm);
    void block_reduce(std::complex<float> *mat, size_t count, MPI_Comm comm);
}

#endif

