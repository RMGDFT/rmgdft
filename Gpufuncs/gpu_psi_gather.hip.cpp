/*
 *
 * Copyright 2026 The RMG Project Developers. See the COPYRIGHT file 
 * at the top-level directory of this distribution or in the current
 * directory.
 * 
 * This file is part of RMG. 
 * RMG is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * any later version.
 *
 * RMG is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
*/

#if HIP_ENABLED
#include "Gpufuncs.h"
#include <hip/hip_runtime.h>
#include <algorithm>

// Used for GPU aware MPI in the PsiUpdate diagonalization routine.
__global__ void assemble_block_matrix(
    double* block_matrix,
    const double* distAij,
    size_t this_block_size_row,
    int nstates,
    int num_blocks,
    int nb,
    int mb,
    int mxllda,
    int ib,
    int myrow,
    int mycol,
    int nprow,
    int npcol)
{
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;

    size_t size_mat = this_block_size_row * nstates;

    if (idx >= size_mat)
        return;

    // Decode the output matrix element.
    size_t i = idx / nstates;
    size_t j = idx % nstates;

    // Default behavior: zero-initialize the output.
    double value = 0.0;

    // Determine which column block contains this element.
    int jb = static_cast<int>(j / nb);
    int j_local = static_cast<int>(j % nb);

    if (jb < num_blocks &&
        myrow == jb % nprow &&
        mycol == ib % npcol)
    {
        int this_block_size_col =
            std::min(nb, nstates - nb * jb);

        if (j_local < this_block_size_col)
        {
            size_t src =
                ((static_cast<size_t>(ib / npcol) * mb + i)
                    * mxllda)
                + static_cast<size_t>(jb / nprow) * nb
                + j_local;

            value = distAij[src];
        }
    }

    block_matrix[idx] = value;
}


__global__ void assemble_block_matrix(
    std::complex<double> *block_matrix,
    const std::complex<double> *distAij,
    size_t this_block_size_row,
    int nstates,
    int num_blocks,
    int nb,
    int mb,
    int mxllda,
    int ib,
    int myrow,
    int mycol,
    int nprow,
    int npcol)
{
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;

    size_t size_mat = this_block_size_row * nstates;

    if (idx >= size_mat)
        return;

    // Decode the output matrix element.
    size_t i = idx / nstates;
    size_t j = idx % nstates;

    // Default behavior: zero-initialize the output.
    std::complex<double> value = 0.0;

    // Determine which column block contains this element.
    int jb = static_cast<int>(j / nb);
    int j_local = static_cast<int>(j % nb);

    if (jb < num_blocks &&
        myrow == jb % nprow &&
        mycol == ib % npcol)
    {
        int this_block_size_col =
            std::min(nb, nstates - nb * jb);

        if (j_local < this_block_size_col)
        {
            size_t src =
                ((static_cast<size_t>(ib / npcol) * mb + i)
                    * mxllda)
                + static_cast<size_t>(jb / nprow) * nb
                + j_local;

            value = distAij[src];
        }
    }

    block_matrix[idx] = value;
}


template <typename T>
void gpu_psi_gather(T *block_matrix, T *distAij, int this_block_size_row, int nstates, int num_blocks, 
                    int nb, int mb, int mxllda, int ib, int myrow, int mycol, int nprow, int npcol)
{

    size_t size_mat = (size_t)this_block_size_row * (size_t)nstates;

    int threads = 256;
    int blocks = (size_mat + threads - 1) / threads;

    hipLaunchKernelGGL(
        assemble_block_matrix,
        dim3(blocks),
        dim3(threads),
        0,
        0,
        block_matrix,
        distAij,
        this_block_size_row,
        nstates,
        num_blocks,
        nb,
        mb,
        mxllda,
        ib,
        myrow,
        mycol,
        nprow,
        npcol
    );

}
template void gpu_psi_gather(double *, double *, int, int, int, int, int, int, int, int, int, int, int);
template void gpu_psi_gather(std::complex<double> *, std::complex<double> *, int, int, int, int, int, int, int, int, int, int, int);

#endif
