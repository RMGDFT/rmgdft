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
   This function is used in the TDDFT code to perform the reduction below on the GPU.

        for(st1 = 0; st1 < numst; st1++)
            for(idx = 0; idx < pbasis; idx++)
                rho_temp[idx] += psi[st1 * pbasis + idx] * xpsi[st1 * pbasis + idx];

*/

#if CUDA_ENABLED
#include <cuda.h>
#include <complex>
#include "Gpufuncs.h"

const int nTPB = 128; // Only 1 block but 128 threads per block

void GpuRhomatrixConvert(double *rho_matrix_dev, double *rho_matrix, double *occ_dev, int numst, int myrank, int nprocs)
{

    size_t tile_size = numst * numst/nprocs;
    for(int i = 0; i < numst/nprocs; i++)
    {
        for(int j = 0; j < numst; j++)
        {
            if(myrank * numst/nprocs + i == j) 
            {
                rho_matrix_dev[myrank * tile_size + i * numst + j] = rho_matrix[i * numst + j] - occ_dev[j];
            }
            else
            {
                rho_matrix_dev[myrank * tile_size + i * numst + j] = rho_matrix[i * numst + j];
            }
        }
    }

}

void GpuRhomatrixConvert(float *rho_matrix_dev, double *rho_matrix, double *occ_dev, int numst, int myrank, int nprocs)
{

    size_t tile_size = numst * numst/nprocs;
    for(int i = 0; i < numst/nprocs; i++)
    {
        for(int j = 0; j < numst; j++)
        {
            if(myrank * numst/nprocs + i == j) 
            {
                rho_matrix_dev[myrank * tile_size + i * numst + j] = rho_matrix[i * numst + j] - occ_dev[j];
            }
            else
            {
                rho_matrix_dev[myrank * tile_size + i * numst + j] = rho_matrix[i * numst + j];
            }
        }
    }

}


void GpuRhomatrixConvert(double *rho_matrix_dev, std::complex<double> *rho_matrix, double *occ_dev, int numst, int myrank, int nprocs)
{

    size_t tile_size = numst * numst/nprocs;
    for(int i = 0; i < numst/nprocs; i++)
    {
        for(int j = 0; j < numst; j++)
        {
            if(myrank * numst/nprocs + i == j) 
            {
                rho_matrix_dev[myrank * tile_size + i * numst + j] = std::real(rho_matrix[i * numst + j]) - occ_dev[j];
            }
            else
            {
                rho_matrix_dev[myrank * tile_size + i * numst + j] = std::real(rho_matrix[i * numst + j]);
            }
        }
    }

}

void GpuRhomatrixConvert(float *rho_matrix_dev, std::complex<double> *rho_matrix, double *occ_dev, int numst, int myrank, int nprocs)
{

    size_t tile_size = numst * numst/nprocs;
    for(int i = 0; i < numst/nprocs; i++)
    {
        for(int j = 0; j < numst; j++)
        {
            if(myrank * numst/nprocs + i == j) 
            {
                rho_matrix_dev[myrank * tile_size + i * numst + j] = std::real(rho_matrix[i * numst + j]) - occ_dev[j];
            }
            else
            {
                rho_matrix_dev[myrank * tile_size + i * numst + j] = std::real(rho_matrix[i * numst + j]);
            }
        }
    }

}


void GpuRhomatrixConvert(std::complex<double> *rho_matrix_dev, std::complex<double> *rho_matrix, double *occ_dev, int numst, int myrank, int nprocs)
{

    size_t tile_size = numst * numst/nprocs;
    for(int i = 0; i < numst/nprocs; i++)
    {
        for(int j = 0; j < numst; j++)
        {
            if(myrank * numst/nprocs + i == j) 
            {
                rho_matrix_dev[myrank * tile_size + i * numst + j] = rho_matrix[i * numst + j] - occ_dev[j];
            }
            else
            {
                rho_matrix_dev[myrank * tile_size + i * numst + j] = rho_matrix[i * numst + j];
            }
        }
    }

}

void GpuRhomatrixConvert(std::complex<float> *rho_matrix_dev, std::complex<double> *rho_matrix, double *occ_dev, int numst, int myrank, int nprocs)
{

    size_t tile_size = numst * numst/nprocs;
    for(int i = 0; i < numst/nprocs; i++)
    {
        for(int j = 0; j < numst; j++)
        {
            if(myrank * numst/nprocs + i == j) 
            {
                rho_matrix_dev[myrank * tile_size + i * numst + j] = rho_matrix[i * numst + j] - occ_dev[j];
            }
            else
            {
                rho_matrix_dev[myrank * tile_size + i * numst + j] = rho_matrix[i * numst + j];
            }
        }
    }

}


#endif

