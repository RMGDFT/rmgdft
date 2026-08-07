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
   transpose of  tiledMM 
   for(int i = 0; i < num_rows; i++)
        for(int j = 0; j < M; j++)
        {
            C[i * M + j] = C_glob[j * M + i + my_rank * num_rows]));
        }

*/

#if HIP_ENABLED


#include <complex>
#include <hip/hip_ext.h>
#include "Gpufuncs.h"


#include <hip/hip_runtime.h>
#include <hip/hip_complex.h>

const int nTPB = 128; // Only 1 block but 128 threads per block

__global__ void tiledM_transpose_simple(int M,
                                     int num_rows,
                                     int row_offset,
                                     double* C,
                                     const double* C_glob)
{
  //for(size_t j = 0;j < M; j += nTPB)
    size_t j = blockIdx.x*nTPB + threadIdx.x;
    
    if(j < M)
    {
        for(int i = 0; i < num_rows; i++)
        {
            size_t idx = i * M + j;
            size_t idx_g =j * M + (i + row_offset);
            C[idx] = C_glob[idx_g];
        }
    }

}

void GpuTiledM_transpose(int M, int num_rows, int my_rank, double *C, double *C_glob)
{

    int nblocks = M / nTPB + 1;
    tiledM_transpose_simple<<<nblocks, nTPB>>>(M, num_rows, my_rank * num_rows, 
            C, C_glob);
}
#endif

