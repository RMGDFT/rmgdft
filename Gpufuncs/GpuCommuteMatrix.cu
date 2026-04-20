/*
 *
 * Copyright 2024 The RMG Project Developers. See the COPYRIGHT file
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

/*
   This function is used in the TDDFT code to perform 
   dP = C = A * B, dP[M, num_rows] is tiledMM of C global (M,M)
   outopyt dP = alpha *[AB - BA]
   for(int i = 0; i < num_rows; i++)
        for(int j = 0; j < M; j++)
        {
            dP[i * M + j] = alpha * (dP[i * M + j] - std::conj(C[j * M + i + my_rank * num_rows]));
        }

*/

#if CUDA_ENABLED

#include <cuda.h>
#include "Gpufuncs.h"


#define TILE_DIM 32
#define BLOCK_ROWS 4

__global__ void commutematrix_simple(int M,
                                     int num_rows,
                                     int row_offset,
                                     cuDoubleComplex alpha,
                                     cuDoubleComplex* dP,
                                     const cuDoubleComplex* C)
{
    int j = blockIdx.x * blockDim.x + threadIdx.x; // column
    int i = blockIdx.y * blockDim.y + threadIdx.y; // local row

    if (i < num_rows && j < M) {
        int global_i = i + row_offset;

        int idx = i * M + j;

        cuDoubleComplex a = dP[idx];

        // direct access (no transpose buffer)
        cuDoubleComplex c = C[j * M + global_i];

        // conjugate + update
        dP[idx] = cuCmul(alpha, cuCsub(a, cuConj(c)));
    }
}

void GpuCommuteMatrix(int M, int num_rows, int my_rank, std::complex<double> alpha, std::complex<double> *dP, std::complex<double> *C)
{

    dim3 block(TILE_DIM, BLOCK_ROWS);
    dim3 grid((M + TILE_DIM - 1) / TILE_DIM,
          (M + TILE_DIM - 1) / TILE_DIM);
    commutematrix_simple<<<grid, block>>>(M, num_rows, my_rank * num_rows, 
             make_cuDoubleComplex(alpha.real(), alpha.imag()),
            (cuDoubleComplex *)dP, (cuDoubleComplex *)C);
}
#endif

