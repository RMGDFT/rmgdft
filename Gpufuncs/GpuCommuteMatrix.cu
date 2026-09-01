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


const int nTPB = 128; // Only 1 block but 128 threads per block

__global__ void commutematrix_simple(int M,
                                     int num_rows,
                                     int row_offset,
                                     cuDoubleComplex alpha,
                                     cuDoubleComplex* dP,
                                     const cuDoubleComplex* C)
{
  //for(size_t j = 0;j < M; j += nTPB)
    size_t j = blockIdx.x*nTPB + threadIdx.x;
    
    if(j < M)
    {
        for(int i = 0; i < num_rows; i++)
        {
            size_t idx = i * M + j;
            size_t idx_g =j * M + (i + row_offset);
            cuDoubleComplex a = dP[idx];
            cuDoubleComplex c = C[idx_g];
            // conjugate + update
            dP[idx] = cuCmul(alpha, cuCsub(a, hipConj(c)));
        }
    }

}

void GpuCommuteMatrix(int M, int num_rows, int my_rank, std::complex<double> alpha, std::complex<double> *dP, std::complex<double> *C)
{

    dim3 block(TILE_DIM, BLOCK_ROWS);
    dim3 grid((M + TILE_DIM - 1) / TILE_DIM,
          (M + TILE_DIM - 1) / TILE_DIM);
    commutematrix_simple<<<nblocks, nTPB>>>(M, num_rows, my_rank * num_rows, 
             make_cuDoubleComplex(alpha.real(), alpha.imag()),
            (cuDoubleComplex *)dP, (cuDoubleComplex *)C);
}
#endif

