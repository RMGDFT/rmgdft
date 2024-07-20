#include "hip/hip_runtime.h"
/*  This file is part of cuTranspose.

    cuTranspose is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    cuTranspose is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with cuTranspose.  If not, see <http://www.gnu.org/licenses/>.

    Copyright 2016 Ibai Gurrutxaga, Javier Muguerza, Jose L. Jodra.
*/

/********************************************
 * Includes                                 *
 ********************************************/
#include "gputranspose.h"
#include "transpose_kernels.h"

/********************************************
 * Public functions                         *
 ********************************************/

__global__
void dev_copy( float*       out,
               const float* in,
               int           np0,
               int           np1,
               int           elements_per_thread )
{
	int x, y, z,
	    ind,
	    i;

	x = threadIdx.x + TILE_SIZE * blockIdx.x;
	y = threadIdx.y + TILE_SIZE * blockIdx.y;
	z = blockIdx.z;

	if( x >= np0 || y >= np1 )
		return;

	ind = x + (y + z * np1) * np0;

	for( i = 0;
	     i < TILE_SIZE && y + i < np1;
	     i += TILE_SIZE / elements_per_thread )
	{
		out[ind + i*np0] = in[ind + i*np0];
	}
}


__global__
void dev_copy( double*       out,
               const double* in,
               int           np0,
               int           np1,
               int           elements_per_thread )
{
	int x, y, z,
	    ind,
	    i;

	x = threadIdx.x + TILE_SIZE * blockIdx.x;
	y = threadIdx.y + TILE_SIZE * blockIdx.y;
	z = blockIdx.z;

	if( x >= np0 || y >= np1 )
		return;

	ind = x + (y + z * np1) * np0;

	for( i = 0;
	     i < TILE_SIZE && y + i < np1;
	     i += TILE_SIZE / elements_per_thread )
	{
		out[ind + i*np0] = in[ind + i*np0];
	}
}

/********************************************
 * Public functions                         *
 ********************************************/

__global__
void dev_transpose_021_ept1( float*       out,
                             const float* in,
                             int           np0,
                             int           np1,
                             int           np2 )
{
	__shared__ float tile[TILE_SIZE][TILE_SIZE + 1];

	int x, y, z,
	    ind_in,
	    ind_out;

	int lx = threadIdx.x,
	    ly = threadIdx.y,
	    bx = blockIdx.x,
	    by = blockIdx.y;

	x = lx + TILE_SIZE * bx;
	y = ly + TILE_SIZE * by;
	z = blockIdx.z;

	ind_in = x + (y + z * np1) * np0;
	ind_out = x + (z + y * np2) * np0;

	if( x < np0 && y < np1 )
	{
		tile[lx][ly] = in[ind_in];
	}

	__syncthreads();

	if( x < np0 && y < np1	 )
	{
		out[ind_out] = tile[lx][ly];
	}
}


__global__
void dev_transpose_021_ept2( float*       out,
                             const float* in,
                             int           np0,
                             int           np1,
                             int           np2 )
{
	__shared__ float tile[TILE_SIZE][TILE_SIZE + 1];

	int x, y, z,
	    ind_in,
	    ind_out;

	int lx = threadIdx.x,
	    ly = threadIdx.y,
	    bx = blockIdx.x,
	    by = blockIdx.y;

	x = lx + TILE_SIZE * bx;
	y = ly + TILE_SIZE * by;
	z = blockIdx.z;

	ind_in = x + (y + z * np1) * np0;
	ind_out = x + (z + y * np2) * np0;

	if( x < np0 && y < np1 )
	{
		tile[lx][ly] = in[ind_in];
		if( y + 8 < np1 )
		{
			tile[lx][ly +  8] = in[ind_in +  8*np0];
		}
	}

	__syncthreads();

	if( x < np0 && y < np1 )
	{
		out[ind_out] = tile[lx][ly];
		if( y + 8 < np1 )
		{
			out[ind_out +  8*np0*np2] = tile[lx][ly + 8];
		}
	}
}


__global__
void dev_transpose_021_ept4( float*       out,
                             const float* in,
                             int           np0,
                             int           np1,
                             int           np2 )
{
__shared__ float tile[TILE_SIZE][TILE_SIZE + 1];

	int x, y, z,
	    ind_in,
	    ind_out;

	int lx = threadIdx.x,
	    ly = threadIdx.y,
	    bx = blockIdx.x,
	    by = blockIdx.y;

	x = lx + TILE_SIZE * bx;
	y = ly + TILE_SIZE * by;
	z = blockIdx.z;

	ind_in = x + (y + z * np1) * np0;
	ind_out = x + (z + y * np2) * np0;

	if( x < np0 && y < np1 )
	{
		tile[lx][ly] = in[ind_in];
		if( y + 4 < np1 )
		{
			tile[lx][ly +  4] = in[ind_in +  4*np0];
			if( y + 8 < np1 )
			{
				tile[lx][ly +  8] = in[ind_in +  8*np0];
				if( y + 12 < np1 )
				{
					tile[lx][ly +  12] = in[ind_in +  12*np0];
				}
			}
		}
	}

	__syncthreads();

	if( x < np0 && y < np1 )
	{
		out[ind_out] = tile[lx][ly];
		if( y + 4 < np1 )
		{
			out[ind_out +  4*np0*np2] = tile[lx][ly + 4];
			if( y + 8 < np1 )
			{
				out[ind_out +  8*np0*np2] = tile[lx][ly + 8];
				if( y + 12 < np1 )
				{
					out[ind_out +  12*np0*np2] = tile[lx][ly + 12];
				}
			}
		}
	}
}


__global__
void dev_transpose_021_in_place( float* data,
                                 int     np0,
                                 int     np1 )
{
	__shared__ float inf_tile[TILE_SIZE][TILE_SIZE + 1];
	__shared__ float sup_tile[TILE_SIZE][TILE_SIZE + 1];

	int x, y, z,
	    ind_inf,
	    ind_sup;

	int lx = threadIdx.x,
	    ly = threadIdx.y,
	    bx = blockIdx.x,
	    by = blockIdx.y;

	z = blockIdx.z;
	y = ly + TILE_SIZE * by;

	if( z > y ) // Thread in diagonal or upper triangle
		return;

	x = lx + TILE_SIZE * bx;

	if( x < np0 && y < np1 )
	{
		ind_inf = x + (y + z * np1) * np0;
		ind_sup = x + (z + y * np1) * np0;
		inf_tile[lx][ly] = data[ind_inf];
		sup_tile[lx][ly] = data[ind_sup];
	}

	__syncthreads();

	if( x < np0 && y < np1 )
	{
		data[ind_sup] = inf_tile[lx][ly];
		data[ind_inf] = sup_tile[lx][ly];
	}
}

/********************************************
 * Includes                                 *
 ********************************************/

/********************************************
 * Public functions                         *
 ********************************************/
/**
 * Kernel that performs a 1,0,2 transpose (order inversion) out of place.
 *
 * The grid of threads must be a 3-dimensional grid of 2-dimensional blocks.
 * The x dimension must match the innermost dimension of the input grid, and
 * the y dimension must match the innermost dimension of the transposed grid.
 * \param [out] out The output array.
 * \param [in] in The input array.
 * \param [in] np0 The size of the first dimension of the input array.
 * \param [in] np1 The size of the second dimension of the input array.
 * \param [in] np2 The size of the third dimension of the input array.
 */

 __global__
 void dev_transpose_102_ept1( float*       out,
                              const float* in,
                              int           np0,
                              int           np1,
                              int           np2 )
{
	__shared__ float tile[TILE_SIZE][TILE_SIZE + 1];

	int x_in, y_in, z,
	    x_out, y_out,
	    ind_in,
	    ind_out;

	int lx = threadIdx.x,
	    ly = threadIdx.y,
	    bx = blockIdx.x,
	    by = blockIdx.y;

	x_in = lx + TILE_SIZE * bx;
	y_in = ly + TILE_SIZE * by;

	z = blockIdx.z;

	x_out = ly + TILE_SIZE * bx;
	y_out = lx + TILE_SIZE * by;


	ind_in = x_in + (y_in + z * np1) * np0;
	ind_out = y_out + (x_out + z * np0) * np1;

	if( x_in < np0 && y_in < np1 )
	{
		tile[lx][ly] = in[ind_in];
	}

	__syncthreads();

	if( x_out < np0 && y_out < np1 )
	{
		out[ind_out] = tile[ly][lx];
	}
}

/**
 * Kernel that performs a 1,0,2 transpose (order inversion) out of place.
 *
 * The grid of threads must be a 3-dimensional grid of 2-dimensional blocks.
 * The x dimension must match the innermost dimension of the input grid, and
 * the y dimension must match the innermost dimension of the transposed grid.
 * \param [out] out The output array.
 * \param [in] in The input array.
 * \param [in] np0 The size of the first dimension of the input array.
 * \param [in] np1 The size of the second dimension of the input array.
 * \param [in] np2 The size of the third dimension of the input array.
 */

__global__
void dev_transpose_102_ept2( float*       out,
                             const float* in,
                             int           np0,
                             int           np1,
                             int           np2 )
{
	__shared__ float tile[TILE_SIZE][TILE_SIZE + 1];

	int x_in, y_in, z,
	    x_out, y_out,
	    ind_in,
	    ind_out;

	int lx = threadIdx.x,
	    ly = threadIdx.y,
	    bx = blockIdx.x,
	    by = blockIdx.y;

	x_in = lx + TILE_SIZE * bx;
	y_in = ly + TILE_SIZE * by;

	z = blockIdx.z;

	x_out = ly + TILE_SIZE * bx;
	y_out = lx + TILE_SIZE * by;


	ind_in = x_in + (y_in + z * np1) * np0;
	ind_out = y_out + (x_out + z * np0) * np1;

	if( x_in < np0 && y_in < np1 )
	{
		tile[lx][ly] = in[ind_in];
		if( y_in + 8 < np1 )
		{
			tile[lx][ly +  8] = in[ind_in +  8*np0];
		}
	}

	__syncthreads();

	if( x_out < np0 && y_out < np1 )
	{
		out[ind_out] = tile[ly][lx];
		if( x_out + 8 < np0 )
		{
			out[ind_out +  8*np1] = tile[ly + 8][lx];
		}
	}
}

/**
 * Kernel that performs a 1,0,2 transpose (order inversion) out of place.
 *
 * The grid of threads must be a 3-dimensional grid of 2-dimensional blocks.
 * The x dimension must match the innermost dimension of the input grid, and
 * the y dimension must match the innermost dimension of the transposed grid.
 * \param [out] out The output array.
 * \param [in] in The input array.
 * \param [in] np0 The size of the first dimension of the input array.
 * \param [in] np1 The size of the second dimension of the input array.
 * \param [in] np2 The size of the third dimension of the input array.
 */

__global__
void dev_transpose_102_ept4( float*       out,
                             const float* in,
                             int           np0,
                             int           np1,
                             int           np2 )
{
	__shared__ float tile[TILE_SIZE][TILE_SIZE + 1];

	int x_in, y_in, z,
	    x_out, y_out,
	    ind_in,
	    ind_out;

	int lx = threadIdx.x,
	    ly = threadIdx.y,
	    bx = blockIdx.x,
	    by = blockIdx.y;

	x_in = lx + TILE_SIZE * bx;
	y_in = ly + TILE_SIZE * by;

	z = blockIdx.z;

	x_out = ly + TILE_SIZE * bx;
	y_out = lx + TILE_SIZE * by;


	ind_in = x_in + (y_in + z * np1) * np0;
	ind_out = y_out + (x_out + z * np0) * np1;

	if( x_in < np0 && y_in < np1 )
	{
		tile[lx][ly] = in[ind_in];
		if( y_in + 4 < np1 )
		{
			tile[lx][ly +  4] = in[ind_in +  4*np0];
			if( y_in + 8 < np1 )
			{
				tile[lx][ly +  8] = in[ind_in +  8*np0];
				if( y_in + 12 < np1 )
				{
					tile[lx][ly +  12] = in[ind_in +  12*np0];
				}
			}
		}
	}

	__syncthreads();

	if( x_out < np0 && y_out < np1 )
	{
		out[ind_out] = tile[ly][lx];
		if( x_out + 4 < np0 )
		{
			out[ind_out +  4*np1] = tile[ly + 4][lx];
			if( x_out + 8 < np0 )
			{
				out[ind_out +  8*np1] = tile[ly + 8][lx];
				if( x_out + 12 < np0 )
				{
					out[ind_out +  12*np1] = tile[ly + 12][lx];
				}
			}
		}
	}
}

/**
 * Kernel that performs a 1,0,2 transpose (order inversion) in place.
 *
 * The grid of threads must be a 3-dimensional grid of 2-dimensional blocks.
 * The x dimension must match the innermost dimension of the input grid, and
 * the y dimension must match the innermost dimension of the transposed grid.
 * \param [in, out] in The data array.
 * \param [in] np0 The size of the first and last dimensions of the input array.
 * \param [in] np2 The size of the third dimension of the input array.
 */

__global__
void dev_transpose_102_in_place( float* data,
                                 int     np0,
                                 int     np2 )
{
	__shared__ float inf_tile[TILE_SIZE][TILE_SIZE + 1];
	__shared__ float sup_tile[TILE_SIZE][TILE_SIZE + 1];

	int x_inf, y_inf, z,
	    x_sup, y_sup,
	    ind_inf,
	    ind_sup;

	int lx = threadIdx.x,
	    ly = threadIdx.y,
	    bx = blockIdx.x,
	    by = blockIdx.y;

	if( bx > by ) // Block in upper triangle
		return;

	x_inf = lx + TILE_SIZE * bx;
	y_inf = ly + TILE_SIZE * by;
	z = blockIdx.z;

	if( x_inf < np0 && y_inf < np0 )
	{
		ind_inf = x_inf + (y_inf + z * np0) * np0;
		inf_tile[lx][ly] = data[ind_inf];
	}

	x_sup = ly + TILE_SIZE * bx;
	y_sup = lx + TILE_SIZE * by;
	if( bx < by ) // Block in lower triangle
	{
		if( x_sup < np0 && y_sup < np0 )
		{
			ind_sup = y_sup + (x_sup + z * np0) * np0;
			sup_tile[lx][ly] = data[ind_sup];
		}
	}
	else // Block in diagonal
		ind_sup = ind_inf;

	__syncthreads();

	if( x_sup < np0 && y_sup < np0 )
	{
		data[ind_sup] = inf_tile[ly][lx];
	}

	if( bx < by )
	{
		if( x_inf < np0 && y_inf < np0 )
		{
			data[ind_inf] = sup_tile[ly][lx];
		}
	}
}


/********************************************
 * Public functions                         *
 ********************************************/

__global__
void dev_transpose_120_ept1( float*       out,
                             const float* in,
                             int           np0,
                             int           np1,
                             int           np2 )
{

	__shared__ float tile[TILE_SIZE][TILE_SIZE + 1];

	int x_in, y_in, z,
	    x_out, y_out,
	    ind_in,
	    ind_out;

	int lx = threadIdx.x,
	    ly = threadIdx.y,
	    bx = blockIdx.x,
	    by = blockIdx.y;

	x_in = lx + TILE_SIZE * bx;
	y_in = ly + TILE_SIZE * by;

	z = blockIdx.z;

	x_out = ly + TILE_SIZE * bx;
	y_out = lx + TILE_SIZE * by;

	ind_in = x_in + (y_in + z * np1) * np0;
	ind_out = y_out + (z + x_out * np2) * np1;

	if( x_in < np0 && y_in < np1 )
	{
		tile[lx][ly] = in[ind_in];
	}

	__syncthreads();

	if( y_out < np1 && x_out < np0 )
	{
		out[ind_out] = tile[ly][lx];
	}
}


__global__
void dev_transpose_120_ept2( float*       out,
                             const float* in,
                             int           np0,
                             int           np1,
                             int           np2 )
{

	__shared__ float tile[TILE_SIZE][TILE_SIZE + 1];

	int x_in, y_in, z,
	    x_out, y_out,
	    ind_in,
	    ind_out;

	int lx = threadIdx.x,
	    ly = threadIdx.y,
	    bx = blockIdx.x,
	    by = blockIdx.y;

	x_in = lx + TILE_SIZE * bx;
	y_in = ly + TILE_SIZE * by;

	z = blockIdx.z;

	x_out = ly + TILE_SIZE * bx;
	y_out = lx + TILE_SIZE * by;

	ind_in = x_in + (y_in + z * np1) * np0;
	ind_out = y_out + (z + x_out * np2) * np1;

	if( x_in < np0 && y_in < np1 )
	{
		tile[lx][ly] = in[ind_in];
		if( y_in + 8 < np1 )
		{
			tile[lx][ly +  8] = in[ind_in +  8*np0];
		}
	}

	__syncthreads();

	if( y_out < np1 && x_out < np0 )
	{
		out[ind_out] = tile[ly][lx];
		if( x_out + 8 < np0 )
		{
			out[ind_out +  8*np1*np2] = tile[ly + 8][lx];
		}
	}

}


__global__
void dev_transpose_120_ept4( float*       out,
                             const float* in,
                             int           np0,
                             int           np1,
                             int           np2 )
{

	__shared__ float tile[TILE_SIZE][TILE_SIZE + 1];

	int x_in, y_in, z,
	    x_out, y_out,
	    ind_in,
	    ind_out;

	int lx = threadIdx.x,
	    ly = threadIdx.y,
	    bx = blockIdx.x,
	    by = blockIdx.y;

	x_in = lx + TILE_SIZE * bx;
	y_in = ly + TILE_SIZE * by;

	z = blockIdx.z;

	x_out = ly + TILE_SIZE * bx;
	y_out = lx + TILE_SIZE * by;

	ind_in = x_in + (y_in + z * np1) * np0;
	ind_out = y_out + (z + x_out * np2) * np1;

	if( x_in < np0 && y_in < np1 )
	{
		tile[lx][ly] = in[ind_in];
		if( y_in + 4 < np1 )
		{
			tile[lx][ly +  4] = in[ind_in +  4*np0];
			if( y_in + 8 < np1 )
			{
				tile[lx][ly +  8] = in[ind_in +  8*np0];
				if( y_in + 12 < np1 )
				{
					tile[lx][ly +  12] = in[ind_in +  12*np0];
				}
			}
		}
	}

	__syncthreads();

	if( y_out < np1 && x_out < np0 )
	{
		out[ind_out] = tile[ly][lx];
		if( x_out + 4 < np0 )
		{
			out[ind_out +  4*np1*np2] = tile[ly + 4][lx];
			if( x_out + 8 < np0 )
			{
				out[ind_out +  8*np1*np2] = tile[ly + 8][lx];
				if( x_out + 12 < np0 )
				{
					out[ind_out +  12*np1*np2] = tile[ly + 12][lx];
				}
			}
		}
	}
}


__global__
void dev_transpose_120_in_place( float* data,
                                 int     np0 )
{
	__shared__ float cube13[BRICK_SIZE][BRICK_SIZE][BRICK_SIZE + 1];
	__shared__ float cube2[BRICK_SIZE][BRICK_SIZE][BRICK_SIZE + 1];

	int x1, y1, z1,
	    x2, y2, z2,
	    x3, y3, z3,
	    ind1, ind2, ind3;

	int lx = threadIdx.x,
	    ly = threadIdx.y,
	    lz = threadIdx.z,
	    bx = blockIdx.x,
	    by = blockIdx.y,
	    bz = blockIdx.z;

	int diagonal = (bx == by && by == bz);

	if( bx > by || bx > bz ||
	    ((bx == by || bx == bz) && by > bz) )
		return;

	x1 = lx + BRICK_SIZE * bx;
	y1 = ly + BRICK_SIZE * by;
	z1 = lz + BRICK_SIZE * bz;

	x2 = ly + BRICK_SIZE * bx;
	y2 = lx + BRICK_SIZE * by;
	z2 = lz + BRICK_SIZE * bz;

	x3 = lz + BRICK_SIZE * bx;
	y3 = ly + BRICK_SIZE * by;
	z3 = lx + BRICK_SIZE * bz;

	ind1 = x1 + (y1 + z1 * np0) * np0;
	ind2 = y2 + (z2 + x2 * np0) * np0;
	ind3 = z3 + (x3 + y3 * np0) * np0;

	// Swap lx and ly to avoid the synchronization commented below.
	if( x1 < np0 && y1 < np0 && z1 < np0 )
		cube13[ly][lx][lz] = data[ind1];
	if( x2 < np0 && y2 < np0 && z2 < np0 )
		if( ! diagonal )
			cube2[lx][ly][lz] = data[ind2];

	__syncthreads();

	if( diagonal )
	{
		if( x1 < np0 && y1 < np0 && z1 < np0 )
			data[ind1] = cube13[lx][lz][ly];
	}
	else
	{
		if( x2 < np0 && y2 < np0 && z2 < np0 )
			data[ind2] = cube13[lx][ly][lz];

		//__syncthreads();
		if( x3 < np0 && y3 < np0 && z3 < np0 )
		{
			cube13[lx][ly][lz] = data[ind3];
			data[ind3] = cube2[ly][lz][lx];
		}
		__syncthreads();

		if( x1 < np0 && y1 < np0 && z1 < np0 )
			data[ind1] = cube13[lz][ly][lx];
	}
}


/********************************************
 * Public functions                         *
 ********************************************/

__global__
void dev_transpose_201_ept1( float*       out,
                             const float* in,
                             int           np0,
                             int           np1,
                             int           np2 )
{
__shared__ float tile[TILE_SIZE][TILE_SIZE + 1];

	int x_in, y, z_in,
	    x_out, z_out,
	    ind_in,
	    ind_out;

	int lx = threadIdx.x,
	    ly = threadIdx.y,
	    bx = blockIdx.x,
	    by = blockIdx.y;

	x_in = lx + TILE_SIZE * bx;
	z_in = ly + TILE_SIZE * by;

	y = blockIdx.z;

	x_out = ly + TILE_SIZE * bx;
	z_out = lx + TILE_SIZE * by;


	ind_in = x_in + (y + z_in * np1) * np0;
	ind_out = z_out + (x_out + y * np0) * np2;

	if( x_in < np0 && z_in < np2 )
	{
		tile[lx][ly] = in[ind_in];
	}

	__syncthreads();

	if( x_out < np0 && z_out < np2 )
	{
		out[ind_out] = tile[ly][lx];
	}
}


__global__
void dev_transpose_201_ept2( float*       out,
                              const float* in,
                              int           np0,
                              int           np1,
                              int           np2 )
{
__shared__ float tile[TILE_SIZE][TILE_SIZE + 1];

	int x_in, y, z_in,
	    x_out, z_out,
	    ind_in,
	    ind_out;

	int lx = threadIdx.x,
	    ly = threadIdx.y,
	    bx = blockIdx.x,
	    by = blockIdx.y;

	x_in = lx + TILE_SIZE * bx;
	z_in = ly + TILE_SIZE * by;

	y = blockIdx.z;

	x_out = ly + TILE_SIZE * bx;
	z_out = lx + TILE_SIZE * by;


	ind_in = x_in + (y + z_in * np1) * np0;
	ind_out = z_out + (x_out + y * np0) * np2;

	if( x_in < np0 && z_in < np2 )
	{
		tile[lx][ly] = in[ind_in];
		if( z_in + 8 < np2 )
		{
			tile[lx][ly +  8] = in[ind_in +  8*np0*np1];
		}
	}

	__syncthreads();

	if( x_out < np0 && z_out < np2 )
	{
		out[ind_out] = tile[ly][lx];
		if( x_out + 8 < np0 )
		{
			out[ind_out +  8*np2] = tile[ly + 8][lx];
		}
	}
}


__global__
void dev_transpose_201_ept4( float*       out,
                             const float* in,
                             int           np0,
                             int           np1,
                             int           np2 )
{
	__shared__ float tile[TILE_SIZE][TILE_SIZE + 1];

	int x_in, y, z_in,
	    x_out, z_out,
	    ind_in,
	    ind_out;

	int lx = threadIdx.x,
	    ly = threadIdx.y,
	    bx = blockIdx.x,
	    by = blockIdx.y;

	x_in = lx + TILE_SIZE * bx;
	z_in = ly + TILE_SIZE * by;

	y = blockIdx.z;

	x_out = ly + TILE_SIZE * bx;
	z_out = lx + TILE_SIZE * by;


	ind_in = x_in + (y + z_in * np1) * np0;
	ind_out = z_out + (x_out + y * np0) * np2;

	if( x_in < np0 && z_in < np2 )
	{
		tile[lx][ly] = in[ind_in];
		if( z_in + 4 < np2 )
		{
			tile[lx][ly +  4] = in[ind_in +  4*np0*np1];
			if( z_in + 8 < np2 )
			{
				tile[lx][ly +  8] = in[ind_in +  8*np0*np1];
				if( z_in + 12 < np2 )
				{
					tile[lx][ly +  12] = in[ind_in +  12*np0*np1];
				}
			}
		}
	}

	__syncthreads();

	if( x_out < np0 && z_out < np2 )
	{
		out[ind_out] = tile[ly][lx];
		if( x_out + 4 < np0 )
		{
			out[ind_out +  4*np2] = tile[ly + 4][lx];
			if( x_out + 8 < np0 )
			{
				out[ind_out +  8*np2] = tile[ly + 8][lx];
				if( x_out + 12 < np0 )
				{
					out[ind_out +  12*np2] = tile[ly + 12][lx];
				}
			}
		}
	}
}


__global__
void dev_transpose_201_in_place( float* data,
                                 int     np0 )
{
	__shared__ float cube12[BRICK_SIZE][BRICK_SIZE][BRICK_SIZE + 1];
	__shared__ float cube3[BRICK_SIZE][BRICK_SIZE][BRICK_SIZE + 1];

	int x1, y1, z1,
	    x2, y2, z2,
	    x3, y3, z3,
	    ind1, ind2, ind3;

	int lx = threadIdx.x,
	    ly = threadIdx.y,
	    lz = threadIdx.z,
	    bx = blockIdx.x,
	    by = blockIdx.y,
	    bz = blockIdx.z;

	int diagonal = (bx == by && by == bz);

	if( bx > by || bx > bz ||
	    ((bx == by || bx == bz) && by > bz) )
		return;

	x1 = lx + BRICK_SIZE * bx;
	y1 = ly + BRICK_SIZE * by;
	z1 = lz + BRICK_SIZE * bz;

	x2 = ly + BRICK_SIZE * bx;
	y2 = lx + BRICK_SIZE * by;
	z2 = lz + BRICK_SIZE * bz;

	x3 = lz + BRICK_SIZE * bx;
	y3 = ly + BRICK_SIZE * by;
	z3 = lx + BRICK_SIZE * bz;

	ind1 = x1 + (y1 + z1 * np0) * np0;
	ind2 = y2 + (z2 + x2 * np0) * np0;
	ind3 = z3 + (x3 + y3 * np0) * np0;

	// Swap lx and lz to avoid the synchronization commented below.
	if( x1 < np0 && y1 < np0 && z1 < np0 )
		cube12[lz][ly][lx] = data[ind1];
	if( x3 < np0 && y3 < np0 && z3 < np0 )
		if( ! diagonal )
			cube3[lx][ly][lz] = data[ind3];

	__syncthreads();

	if( diagonal )
	{
		if( x1 < np0 && y1 < np0 && z1 < np0 )
			data[ind1] = cube12[lx][lz][ly];
	}
	else
	{
		if( x3 < np0 && y3 < np0 && z3 < np0 )
			data[ind3] = cube12[lx][ly][lz];

		//__syncthreads();
		if( x2 < np0 && y2 < np0 && z2 < np0 )
		{
			cube12[lx][ly][lz] = data[ind2];
			data[ind2] = cube3[lz][lx][ly];
		}
		__syncthreads();

		if( x1 < np0 && y1 < np0 && z1 < np0 )
			data[ind1] = cube12[ly][lx][lz];
	}
}

 
/********************************************
 * Public functions                         *
 ********************************************/

__global__
void dev_transpose_210_ept1( float*       out,
                             const float* in,
                             int           np0,
                             int           np1,
                             int           np2 )
{

	__shared__ float tile[TILE_SIZE][TILE_SIZE + 1];
	
	int x_in, y, z_in,
	    x_out, z_out,
	    ind_in,
	    ind_out;
	
	int lx = threadIdx.x,
	    ly = threadIdx.y,
	    bx = blockIdx.x,
	    by = blockIdx.y;

	x_in = lx + TILE_SIZE * bx;
	z_in = ly + TILE_SIZE * by;

	y = blockIdx.z;

	x_out = ly + TILE_SIZE * bx;
	z_out = lx + TILE_SIZE * by;


	ind_in = x_in + (y + z_in * np1) * np0;
	ind_out = z_out + (y + x_out * np1) * np2;

	if( x_in < np0 && z_in < np2 )
	{
			tile[lx][ly] = in[ind_in];
	}	

	__syncthreads();

	if( z_out < np2 && x_out < np0 )
	{
		out[ind_out] = tile[ly][lx];
	}

}


__global__
void dev_transpose_210_ept2( float*       out,
                             const float* in,
                             int           np0,
                             int           np1,
                             int           np2 )
{

	__shared__ float tile[TILE_SIZE][TILE_SIZE + 1];
	
	int x_in, y, z_in,
	    x_out, z_out,
	    ind_in,
	    ind_out;
	
	int lx = threadIdx.x,
	    ly = threadIdx.y,
	    bx = blockIdx.x,
	    by = blockIdx.y;
		
	x_in = lx + TILE_SIZE * bx;
	z_in = ly + TILE_SIZE * by;

	y = blockIdx.z;

	x_out = ly + TILE_SIZE * bx;
	z_out = lx + TILE_SIZE * by;

	ind_in = x_in + (y + z_in * np1) * np0;
	ind_out = z_out + (y + x_out * np1) * np2;

	if( x_in < np0 && z_in < np2 )
	{
		tile[lx][ly] = in[ind_in];
		if( z_in + 8 < np2 )
		{
			tile[lx][ly +  8] = in[ind_in +  8*np0*np1];
		}
	}
	
	__syncthreads();

	if( z_out < np2 && x_out < np0 )
	{
		out[ind_out] = tile[ly][lx];
		if( x_out + 8 < np0 )
		{
			out[ind_out +  8*np2*np1] = tile[ly + 8][lx];
		}
	}

}


__global__
void dev_transpose_210_ept4( float*       out,
                             const float* in,
                             int           np0,
                             int           np1,
                             int           np2 )
{

	__shared__ float tile[TILE_SIZE][TILE_SIZE + 1];
	
	int x_in, y, z_in,
	    x_out, z_out,
	    ind_in,
	    ind_out;
	
	int lx = threadIdx.x,
	    ly = threadIdx.y,
	    bx = blockIdx.x,
	    by = blockIdx.y;
		
	x_in = lx + TILE_SIZE * bx;
	z_in = ly + TILE_SIZE * by;

	y = blockIdx.z;

	x_out = ly + TILE_SIZE * bx;
	z_out = lx + TILE_SIZE * by;

	ind_in = x_in + (y + z_in * np1) * np0;
	ind_out = z_out + (y + x_out * np1) * np2;

	if( x_in < np0 && z_in < np2 )
	{
		tile[lx][ly] = in[ind_in];
		if( z_in + 4 < np2 )
		{
			tile[lx][ly +  4] = in[ind_in +  4*np0*np1];
			if( z_in + 8 < np2 )
			{
				tile[lx][ly +  8] = in[ind_in +  8*np0*np1];
				if( z_in + 12 < np2 )
				{
					tile[lx][ly + 12] = in[ind_in + 12*np0*np1];
				}
			}
		}
	}
	
	__syncthreads();
	
	if( z_out < np2 && x_out < np0 )
	{
		out[ind_out] = tile[ly][lx];
		if( x_out + 4 < np0 )
		{
			out[ind_out +  4*np2*np1] = tile[ly +  4][lx];
			if( x_out + 8 < np0 )
			{
				out[ind_out +  8*np2*np1] = tile[ly +  8][lx];
				if( x_out + 12 < np0 )
				{
					out[ind_out + 12*np2*np1] = tile[ly + 12][lx];
				}
			}
		}
	}
}


__global__
void dev_transpose_210_in_place( float* data,
                                 int     np0,
                                 int     np1 )
{
	__shared__ float inf_tile[TILE_SIZE][TILE_SIZE + 1];
	__shared__ float sup_tile[TILE_SIZE][TILE_SIZE + 1];
	
	int x_inf, y, z_inf,
	    x_sup, z_sup,
	    ind_inf,
	    ind_sup;

	int lx = threadIdx.x,
	    ly = threadIdx.y,
	    bx = blockIdx.x,
	    by = blockIdx.y;

	if( bx > by ) // Block in upper triangle
		return;

	x_inf = lx + TILE_SIZE * bx;
	z_inf = ly + TILE_SIZE * by;
	y = blockIdx.z;

	if( x_inf < np0 && z_inf < np0 )
	{
		ind_inf = x_inf + (y + z_inf * np1) * np0;
		inf_tile[lx][ly] = data[ind_inf];
	}

	x_sup = ly + TILE_SIZE * bx;
	z_sup = lx + TILE_SIZE * by;
	if( bx < by ) // Block in lower triangle
	{
		if( x_sup < np0 && z_sup < np0 )
		{
			ind_sup = z_sup + (y + x_sup * np1) * np0;
			sup_tile[lx][ly] = data[ind_sup];
		}
	}
	else // Block in diagonal
		ind_sup = ind_inf;

	__syncthreads();

	if( x_sup < np0 && z_sup < np0 )
	{
		data[ind_sup] = inf_tile[ly][lx];
	}

	if( bx < by )
	{
		if( x_inf < np0 && z_inf < np0 )
		{
			data[ind_inf] = sup_tile[ly][lx];
		}
	}
}
/*  This file is part of cuTranspose.

    cuTranspose is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    cuTranspose is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with cuTranspose.  If not, see <http://www.gnu.org/licenses/>.

    Copyright 2016 Ibai Gurrutxaga, Javier Muguerza, Jose L. Jodra.
*/

/********************************************
 * Includes                                 *
 ********************************************/
#include <stdio.h>
#include <stdlib.h>

/********************************************
 * Private function prototypes              *
 ********************************************/
static void set_grid_dims( const int* size,
                           int        d2,
                           dim3*      block_size,
                           dim3*      num_blocks,
                           int        elements_per_thread);
static void set_grid_dims_cube( const int* size,
                                dim3*      block_size,
                                dim3*      num_blocks,
                                int        elements_per_thread );
static int valid_parameters( int        in_place,
                             const int* size,
                             const int* permutation,
                             int        elements_per_thread );

/********************************************
 * Exported functions                         *
 ********************************************/
/**
 * Perform a transposition of a 3D array in a GPU. The first dimension
 * is considered the innermost one. If the \p input and \ output pointers
 * are equal an in-place transposition will be performed.
 *
 * \param[out] output A pointer to the allocated memory space in the device
 * where the transposed array will be stored.
 * \param[in] input A pointer to device memory where the input data is stored.
 * \param[in] size A 3 element array with the size of the input data
 * on each dimension.
 * \param[in] permutation An array with a permutation specifying the
 * transposition to be performed.
 * \param[in] elements_per_thread The number of elements that a GPU thread must transpose.
 * It will be ignored if \p in_place is true.
 * \return 0 on success and -1 otherwise.
 */
int cut_transpose3d(float *, const float *, const int*, const int*, int, gpuStream_t);
//int cut_transpose3d(std::complex<float> *, const std::complex<float> *, const int*, const int*, int, gpuStream_t);
//int cut_transpose3d(double *, const double *, const int*, const int*, int, gpuStream_t);
//int cut_transpose3d(std::complex<double> *, const std::complex<double> *, const int*, const int*, int, gpuStream_t);


int cut_transpose3d( float*       output,
                     const float* input,
                     const int*    size,
                     const int*    permutation,
                     int           elements_per_thread,
                     gpuStream_t stream)

{
	dim3 num_blocks,
	     block_size;
	int in_place = output == input;

	if( !valid_parameters( in_place, size, permutation, elements_per_thread ) )
		return -1;

	if(( permutation[0] == 1 && permutation[1] == 2 && permutation[2] == 0 && in_place ) ||
	   ( permutation[0] == 2 && permutation[1] == 0 && permutation[2] == 1 && in_place ))
		set_grid_dims_cube( size,
		                    &block_size,
		                    &num_blocks,
		                    elements_per_thread );
	else
		set_grid_dims( size,
		               permutation[0],
		               &block_size,
		               &num_blocks,
		               elements_per_thread );

	if( permutation[0] == 0 && permutation[1] == 1 && permutation[2] == 2 )
	{
			dev_copy<<< num_blocks, block_size >>>( output,
			                                        input,
			                                        size[0],
			                                        size[1], 
			                                        elements_per_thread );
	}
	else if( permutation[0] == 0 && permutation[1] == 2 && permutation[2] == 1 )
	{
		if( in_place )
		{
			if( size[1] != size[2] )
			{
				fprintf( stderr, "This in-place transposition requires equal dimensions for "
								 "the second and third dimensions.\n" );
				return -1;
			}
			dev_transpose_021_in_place<<< num_blocks, block_size >>>( output,
			                                                          size[0],
			                                                          size[1] );
		}
		else {
			switch(elements_per_thread)
			{
			case 1:
				dev_transpose_021_ept1<<< num_blocks, block_size >>>( output,
				                                                      input,
				                                                      size[0],
				                                                      size[1],
				                                                      size[2]);
				break;
			case 2:
				dev_transpose_021_ept2<<< num_blocks, block_size >>>( output,
				                                                      input,
				                                                      size[0],
				                                                      size[1],
				                                                      size[2]);
				break;
			case 4:
				dev_transpose_021_ept4<<< num_blocks, block_size >>>( output,
				                                                      input,
				                                                      size[0],
				                                                      size[1],
				                                                      size[2]);
				break;
			}
		}
	}
	else if( permutation[0] == 1 && permutation[1] == 0 && permutation[2] == 2 )
	{
		if( in_place )
		{
			if( size[0] != size[1] )
			{
				fprintf( stderr, "The 102 in-place transposition requires equal dimensions for "
								 "the first and second dimensions.\n" );
				return -1;
			}
			dev_transpose_102_in_place<<< num_blocks, block_size >>>( output,
                                                                      size[0],
                                                                      size[2] );
		}
		else {
			switch(elements_per_thread)
			{
			case 1:
				dev_transpose_102_ept1<<< num_blocks, block_size >>>( output,
				                                                      input,
				                                                      size[0],
				                                                      size[1],
				                                                      size[2]);
				break;
			case 2:
				dev_transpose_102_ept2<<< num_blocks, block_size >>>( output,
				                                                      input,
				                                                      size[0],
				                                                      size[1],
				                                                      size[2]);
				break;
			case 4:
				dev_transpose_102_ept4<<< num_blocks, block_size >>>( output,
				                                                      input,
				                                                      size[0],
				                                                      size[1],
				                                                      size[2]);
				break;
			}
		}
	}
	else if( permutation[0] == 1 && permutation[1] == 2 && permutation[2] == 0 )
	{
		if( in_place )
		{
			if( size[0] != size[1] || size[0] != size[2] )
			{
				fprintf( stderr, "The 120 in-place transposition requires a cubic input.\n" );
				return -1;
			}
			dev_transpose_120_in_place<<< num_blocks, block_size >>>( output,
			                                                          size[0] );
		}
		else {
			switch(elements_per_thread)
			{
			case 1:
				dev_transpose_120_ept1<<< num_blocks, block_size >>>( output,
				                                                      input,
				                                                      size[0],
				                                                      size[1],
				                                                      size[2]);
				break;
			case 2:
				dev_transpose_120_ept2<<< num_blocks, block_size >>>( output,
				                                                      input,
				                                                      size[0],
				                                                      size[1],
				                                                      size[2]);
				break;
			case 4:
				dev_transpose_120_ept4<<< num_blocks, block_size >>>( output,
				                                                      input,
				                                                      size[0],
				                                                      size[1],
				                                                      size[2]);
				break;
			}
		}
	}
	else if( permutation[0] == 2 && permutation[1] == 0 && permutation[2] == 1 )
	{
		if( in_place )
		{
			if( size[0] != size[1] || size[0] != size[2] )
			{
				fprintf( stderr, "The 201 in-place transposition requires a cubic input.\n" );
				return -1;
			}
			dev_transpose_201_in_place<<< num_blocks, block_size >>>( output,
			                                                          size[0] );
		}
		else {
			switch(elements_per_thread)
			{
			case 1:
				dev_transpose_201_ept1<<< num_blocks, block_size >>>( output,
				                                                      input,
				                                                      size[0],
				                                                      size[1],
				                                                      size[2]);
				break;
			case 2:
				dev_transpose_201_ept2<<< num_blocks, block_size >>>( output,
				                                                      input,
				                                                      size[0],
				                                                      size[1],
				                                                      size[2]);
				break;
			case 4:
				dev_transpose_201_ept4<<< num_blocks, block_size >>>( output,
				                                                      input,
				                                                      size[0],
				                                                      size[1],
				                                                      size[2]);
				break;
			}
		}
	}
	else if( permutation[0] == 2 && permutation[1] == 1 && permutation[2] == 0 )
	{
		if( in_place )
		{
			if( size[0] != size[2] )
			{
				fprintf( stderr, "The 210 in-place transposition requires equal dimensions for "
								 "the first and third dimensions.\n" );
				return -1;
			}
			dev_transpose_210_in_place<<< num_blocks, block_size >>>( output,
			                                                          size[0],
			                                                          size[1] );
		}
		else {
			switch(elements_per_thread)
			{
			case 1:
				dev_transpose_210_ept1<<< num_blocks, block_size >>>( output,
				                                                      input,
				                                                      size[0],
				                                                      size[1],
				                                                      size[2]);
				break;
			case 2:
				dev_transpose_210_ept2<<< num_blocks, block_size >>>( output,
				                                                      input,
				                                                      size[0],
				                                                      size[1],
				                                                      size[2]);
				break;
			case 4:
				dev_transpose_210_ept4<<< num_blocks, block_size >>>( output,
				                                                      input,
				                                                      size[0],
				                                                      size[1],
				                                                      size[2]);
				break;
			}
		}
	}
	
	return 0;
}

/********************************************
 * Private functions                        *
 ********************************************/
static void set_grid_dims( const int* size,
                           int        d2,
                           dim3*      block_size,
                           dim3*      num_blocks,
                           int        elements_per_thread)
{
	block_size->x = TILE_SIZE;
	block_size->y = TILE_SIZE / elements_per_thread;
	block_size->z = 1;
	num_blocks->x = size[0] / TILE_SIZE;
	if( size[0] % TILE_SIZE != 0 )
		num_blocks->x++;
	if( d2 == 0 )
		d2 = 1;
	num_blocks->y = size[d2] / TILE_SIZE;
	if( size[d2] % TILE_SIZE != 0 )
		num_blocks->y++;
	num_blocks->z = size[(d2 == 1) ? 2 : 1];
}

static void set_grid_dims_cube( const int* size,
                                dim3*      block_size,
                                dim3*      num_blocks,
                                int        elements_per_thread )
{
	block_size->x = BRICK_SIZE;
	block_size->y = BRICK_SIZE;
	block_size->z = BRICK_SIZE / elements_per_thread;
	num_blocks->x = size[0] / BRICK_SIZE;
	if( size[0] % BRICK_SIZE != 0 )
		num_blocks->x++;
	num_blocks->y = size[1] / BRICK_SIZE;
	if( size[1] % BRICK_SIZE != 0 )
		num_blocks->y++;
	num_blocks->z = size[2] / BRICK_SIZE;
	if( size[2] % BRICK_SIZE != 0 )
		num_blocks->z++;
}

static int valid_parameters( int        in_place,
                             const int* size,
                             const int* permutation,
                             int        elements_per_thread )
{
	int dims[] = { 0, 0, 0 },
        i;

	if( in_place && elements_per_thread != 1 )
		return 0;
	if( size == NULL || permutation == NULL )
		return 0;
	if( size[0] < 2 || size[1] < 2 || size[2] < 2 )
		return 0;

	for( i = 0; i < 3; i++ )
	{
		if( permutation[i] < 0 || permutation[i] > 2 )
			return 0;
		else if( dims[permutation[i]] == 1 )
			return 0;
		else
			dims[permutation[i]] = 1;
	}

	return 1;
}


/********************************************
 * Public functions                         *
 ********************************************/

__global__
void dev_copy( std::complex<float>*       out,
               const std::complex<float>* in,
               int           np0,
               int           np1,
               int           elements_per_thread )
{
	int x, y, z,
	    ind,
	    i;

	x = threadIdx.x + TILE_SIZE * blockIdx.x;
	y = threadIdx.y + TILE_SIZE * blockIdx.y;
	z = blockIdx.z;

	if( x >= np0 || y >= np1 )
		return;

	ind = x + (y + z * np1) * np0;

	for( i = 0;
	     i < TILE_SIZE && y + i < np1;
	     i += TILE_SIZE / elements_per_thread )
	{
		out[ind + i*np0] = in[ind + i*np0];
	}
}

/********************************************
 * Public functions                         *
 ********************************************/

__global__
void dev_transpose_021_ept1( double*       out,
                             const double* in,
                             int           np0,
                             int           np1,
                             int           np2 )
{
	__shared__ double tile[TILE_SIZE][TILE_SIZE + 1];

	int x, y, z,
	    ind_in,
	    ind_out;

	int lx = threadIdx.x,
	    ly = threadIdx.y,
	    bx = blockIdx.x,
	    by = blockIdx.y;

	x = lx + TILE_SIZE * bx;
	y = ly + TILE_SIZE * by;
	z = blockIdx.z;

	ind_in = x + (y + z * np1) * np0;
	ind_out = x + (z + y * np2) * np0;

	if( x < np0 && y < np1 )
	{
		tile[lx][ly] = in[ind_in];
	}

	__syncthreads();

	if( x < np0 && y < np1	 )
	{
		out[ind_out] = tile[lx][ly];
	}
}


__global__
void dev_transpose_021_ept2( double*       out,
                             const double* in,
                             int           np0,
                             int           np1,
                             int           np2 )
{
	__shared__ double tile[TILE_SIZE][TILE_SIZE + 1];

	int x, y, z,
	    ind_in,
	    ind_out;

	int lx = threadIdx.x,
	    ly = threadIdx.y,
	    bx = blockIdx.x,
	    by = blockIdx.y;

	x = lx + TILE_SIZE * bx;
	y = ly + TILE_SIZE * by;
	z = blockIdx.z;

	ind_in = x + (y + z * np1) * np0;
	ind_out = x + (z + y * np2) * np0;

	if( x < np0 && y < np1 )
	{
		tile[lx][ly] = in[ind_in];
		if( y + 8 < np1 )
		{
			tile[lx][ly +  8] = in[ind_in +  8*np0];
		}
	}

	__syncthreads();

	if( x < np0 && y < np1 )
	{
		out[ind_out] = tile[lx][ly];
		if( y + 8 < np1 )
		{
			out[ind_out +  8*np0*np2] = tile[lx][ly + 8];
		}
	}
}


__global__
void dev_transpose_021_ept4( double*       out,
                             const double* in,
                             int           np0,
                             int           np1,
                             int           np2 )
{
__shared__ double tile[TILE_SIZE][TILE_SIZE + 1];

	int x, y, z,
	    ind_in,
	    ind_out;

	int lx = threadIdx.x,
	    ly = threadIdx.y,
	    bx = blockIdx.x,
	    by = blockIdx.y;

	x = lx + TILE_SIZE * bx;
	y = ly + TILE_SIZE * by;
	z = blockIdx.z;

	ind_in = x + (y + z * np1) * np0;
	ind_out = x + (z + y * np2) * np0;

	if( x < np0 && y < np1 )
	{
		tile[lx][ly] = in[ind_in];
		if( y + 4 < np1 )
		{
			tile[lx][ly +  4] = in[ind_in +  4*np0];
			if( y + 8 < np1 )
			{
				tile[lx][ly +  8] = in[ind_in +  8*np0];
				if( y + 12 < np1 )
				{
					tile[lx][ly +  12] = in[ind_in +  12*np0];
				}
			}
		}
	}

	__syncthreads();

	if( x < np0 && y < np1 )
	{
		out[ind_out] = tile[lx][ly];
		if( y + 4 < np1 )
		{
			out[ind_out +  4*np0*np2] = tile[lx][ly + 4];
			if( y + 8 < np1 )
			{
				out[ind_out +  8*np0*np2] = tile[lx][ly + 8];
				if( y + 12 < np1 )
				{
					out[ind_out +  12*np0*np2] = tile[lx][ly + 12];
				}
			}
		}
	}
}


__global__
void dev_transpose_021_in_place( double* data,
                                 int     np0,
                                 int     np1 )
{
	__shared__ double inf_tile[TILE_SIZE][TILE_SIZE + 1];
	__shared__ double sup_tile[TILE_SIZE][TILE_SIZE + 1];

	int x, y, z,
	    ind_inf,
	    ind_sup;

	int lx = threadIdx.x,
	    ly = threadIdx.y,
	    bx = blockIdx.x,
	    by = blockIdx.y;

	z = blockIdx.z;
	y = ly + TILE_SIZE * by;

	if( z > y ) // Thread in diagonal or upper triangle
		return;

	x = lx + TILE_SIZE * bx;

	if( x < np0 && y < np1 )
	{
		ind_inf = x + (y + z * np1) * np0;
		ind_sup = x + (z + y * np1) * np0;
		inf_tile[lx][ly] = data[ind_inf];
		sup_tile[lx][ly] = data[ind_sup];
	}

	__syncthreads();

	if( x < np0 && y < np1 )
	{
		data[ind_sup] = inf_tile[lx][ly];
		data[ind_inf] = sup_tile[lx][ly];
	}
}

/********************************************
 * Public functions                         *
 ********************************************/
/**
 * Kernel that performs a 1,0,2 transpose (order inversion) out of place.
 *
 * The grid of threads must be a 3-dimensional grid of 2-dimensional blocks.
 * The x dimension must match the innermost dimension of the input grid, and
 * the y dimension must match the innermost dimension of the transposed grid.
 * \param [out] out The output array.
 * \param [in] in The input array.
 * \param [in] np0 The size of the first dimension of the input array.
 * \param [in] np1 The size of the second dimension of the input array.
 * \param [in] np2 The size of the third dimension of the input array.
 */

 __global__
 void dev_transpose_102_ept1( double*       out,
                              const double* in,
                              int           np0,
                              int           np1,
                              int           np2 )
{
	__shared__ double tile[TILE_SIZE][TILE_SIZE + 1];

	int x_in, y_in, z,
	    x_out, y_out,
	    ind_in,
	    ind_out;

	int lx = threadIdx.x,
	    ly = threadIdx.y,
	    bx = blockIdx.x,
	    by = blockIdx.y;

	x_in = lx + TILE_SIZE * bx;
	y_in = ly + TILE_SIZE * by;

	z = blockIdx.z;

	x_out = ly + TILE_SIZE * bx;
	y_out = lx + TILE_SIZE * by;


	ind_in = x_in + (y_in + z * np1) * np0;
	ind_out = y_out + (x_out + z * np0) * np1;

	if( x_in < np0 && y_in < np1 )
	{
		tile[lx][ly] = in[ind_in];
	}

	__syncthreads();

	if( x_out < np0 && y_out < np1 )
	{
		out[ind_out] = tile[ly][lx];
	}
}

/**
 * Kernel that performs a 1,0,2 transpose (order inversion) out of place.
 *
 * The grid of threads must be a 3-dimensional grid of 2-dimensional blocks.
 * The x dimension must match the innermost dimension of the input grid, and
 * the y dimension must match the innermost dimension of the transposed grid.
 * \param [out] out The output array.
 * \param [in] in The input array.
 * \param [in] np0 The size of the first dimension of the input array.
 * \param [in] np1 The size of the second dimension of the input array.
 * \param [in] np2 The size of the third dimension of the input array.
 */

__global__
void dev_transpose_102_ept2( double*       out,
                             const double* in,
                             int           np0,
                             int           np1,
                             int           np2 )
{
	__shared__ double tile[TILE_SIZE][TILE_SIZE + 1];

	int x_in, y_in, z,
	    x_out, y_out,
	    ind_in,
	    ind_out;

	int lx = threadIdx.x,
	    ly = threadIdx.y,
	    bx = blockIdx.x,
	    by = blockIdx.y;

	x_in = lx + TILE_SIZE * bx;
	y_in = ly + TILE_SIZE * by;

	z = blockIdx.z;

	x_out = ly + TILE_SIZE * bx;
	y_out = lx + TILE_SIZE * by;


	ind_in = x_in + (y_in + z * np1) * np0;
	ind_out = y_out + (x_out + z * np0) * np1;

	if( x_in < np0 && y_in < np1 )
	{
		tile[lx][ly] = in[ind_in];
		if( y_in + 8 < np1 )
		{
			tile[lx][ly +  8] = in[ind_in +  8*np0];
		}
	}

	__syncthreads();

	if( x_out < np0 && y_out < np1 )
	{
		out[ind_out] = tile[ly][lx];
		if( x_out + 8 < np0 )
		{
			out[ind_out +  8*np1] = tile[ly + 8][lx];
		}
	}
}

/**
 * Kernel that performs a 1,0,2 transpose (order inversion) out of place.
 *
 * The grid of threads must be a 3-dimensional grid of 2-dimensional blocks.
 * The x dimension must match the innermost dimension of the input grid, and
 * the y dimension must match the innermost dimension of the transposed grid.
 * \param [out] out The output array.
 * \param [in] in The input array.
 * \param [in] np0 The size of the first dimension of the input array.
 * \param [in] np1 The size of the second dimension of the input array.
 * \param [in] np2 The size of the third dimension of the input array.
 */

__global__
void dev_transpose_102_ept4( double*       out,
                             const double* in,
                             int           np0,
                             int           np1,
                             int           np2 )
{
	__shared__ double tile[TILE_SIZE][TILE_SIZE + 1];

	int x_in, y_in, z,
	    x_out, y_out,
	    ind_in,
	    ind_out;

	int lx = threadIdx.x,
	    ly = threadIdx.y,
	    bx = blockIdx.x,
	    by = blockIdx.y;

	x_in = lx + TILE_SIZE * bx;
	y_in = ly + TILE_SIZE * by;

	z = blockIdx.z;

	x_out = ly + TILE_SIZE * bx;
	y_out = lx + TILE_SIZE * by;


	ind_in = x_in + (y_in + z * np1) * np0;
	ind_out = y_out + (x_out + z * np0) * np1;

	if( x_in < np0 && y_in < np1 )
	{
		tile[lx][ly] = in[ind_in];
		if( y_in + 4 < np1 )
		{
			tile[lx][ly +  4] = in[ind_in +  4*np0];
			if( y_in + 8 < np1 )
			{
				tile[lx][ly +  8] = in[ind_in +  8*np0];
				if( y_in + 12 < np1 )
				{
					tile[lx][ly +  12] = in[ind_in +  12*np0];
				}
			}
		}
	}

	__syncthreads();

	if( x_out < np0 && y_out < np1 )
	{
		out[ind_out] = tile[ly][lx];
		if( x_out + 4 < np0 )
		{
			out[ind_out +  4*np1] = tile[ly + 4][lx];
			if( x_out + 8 < np0 )
			{
				out[ind_out +  8*np1] = tile[ly + 8][lx];
				if( x_out + 12 < np0 )
				{
					out[ind_out +  12*np1] = tile[ly + 12][lx];
				}
			}
		}
	}
}

/**
 * Kernel that performs a 1,0,2 transpose (order inversion) in place.
 *
 * The grid of threads must be a 3-dimensional grid of 2-dimensional blocks.
 * The x dimension must match the innermost dimension of the input grid, and
 * the y dimension must match the innermost dimension of the transposed grid.
 * \param [in, out] in The data array.
 * \param [in] np0 The size of the first and last dimensions of the input array.
 * \param [in] np2 The size of the third dimension of the input array.
 */

__global__
void dev_transpose_102_in_place( double* data,
                                 int     np0,
                                 int     np2 )
{
	__shared__ double inf_tile[TILE_SIZE][TILE_SIZE + 1];
	__shared__ double sup_tile[TILE_SIZE][TILE_SIZE + 1];

	int x_inf, y_inf, z,
	    x_sup, y_sup,
	    ind_inf,
	    ind_sup;

	int lx = threadIdx.x,
	    ly = threadIdx.y,
	    bx = blockIdx.x,
	    by = blockIdx.y;

	if( bx > by ) // Block in upper triangle
		return;

	x_inf = lx + TILE_SIZE * bx;
	y_inf = ly + TILE_SIZE * by;
	z = blockIdx.z;

	if( x_inf < np0 && y_inf < np0 )
	{
		ind_inf = x_inf + (y_inf + z * np0) * np0;
		inf_tile[lx][ly] = data[ind_inf];
	}

	x_sup = ly + TILE_SIZE * bx;
	y_sup = lx + TILE_SIZE * by;
	if( bx < by ) // Block in lower triangle
	{
		if( x_sup < np0 && y_sup < np0 )
		{
			ind_sup = y_sup + (x_sup + z * np0) * np0;
			sup_tile[lx][ly] = data[ind_sup];
		}
	}
	else // Block in diagonal
		ind_sup = ind_inf;

	__syncthreads();

	if( x_sup < np0 && y_sup < np0 )
	{
		data[ind_sup] = inf_tile[ly][lx];
	}

	if( bx < by )
	{
		if( x_inf < np0 && y_inf < np0 )
		{
			data[ind_inf] = sup_tile[ly][lx];
		}
	}
}


/********************************************
 * Public functions                         *
 ********************************************/

__global__
void dev_transpose_120_ept1( double*       out,
                             const double* in,
                             int           np0,
                             int           np1,
                             int           np2 )
{

	__shared__ double tile[TILE_SIZE][TILE_SIZE + 1];

	int x_in, y_in, z,
	    x_out, y_out,
	    ind_in,
	    ind_out;

	int lx = threadIdx.x,
	    ly = threadIdx.y,
	    bx = blockIdx.x,
	    by = blockIdx.y;

	x_in = lx + TILE_SIZE * bx;
	y_in = ly + TILE_SIZE * by;

	z = blockIdx.z;

	x_out = ly + TILE_SIZE * bx;
	y_out = lx + TILE_SIZE * by;

	ind_in = x_in + (y_in + z * np1) * np0;
	ind_out = y_out + (z + x_out * np2) * np1;

	if( x_in < np0 && y_in < np1 )
	{
		tile[lx][ly] = in[ind_in];
	}

	__syncthreads();

	if( y_out < np1 && x_out < np0 )
	{
		out[ind_out] = tile[ly][lx];
	}
}


__global__
void dev_transpose_120_ept2( double*       out,
                             const double* in,
                             int           np0,
                             int           np1,
                             int           np2 )
{

	__shared__ double tile[TILE_SIZE][TILE_SIZE + 1];

	int x_in, y_in, z,
	    x_out, y_out,
	    ind_in,
	    ind_out;

	int lx = threadIdx.x,
	    ly = threadIdx.y,
	    bx = blockIdx.x,
	    by = blockIdx.y;

	x_in = lx + TILE_SIZE * bx;
	y_in = ly + TILE_SIZE * by;

	z = blockIdx.z;

	x_out = ly + TILE_SIZE * bx;
	y_out = lx + TILE_SIZE * by;

	ind_in = x_in + (y_in + z * np1) * np0;
	ind_out = y_out + (z + x_out * np2) * np1;

	if( x_in < np0 && y_in < np1 )
	{
		tile[lx][ly] = in[ind_in];
		if( y_in + 8 < np1 )
		{
			tile[lx][ly +  8] = in[ind_in +  8*np0];
		}
	}

	__syncthreads();

	if( y_out < np1 && x_out < np0 )
	{
		out[ind_out] = tile[ly][lx];
		if( x_out + 8 < np0 )
		{
			out[ind_out +  8*np1*np2] = tile[ly + 8][lx];
		}
	}

}


__global__
void dev_transpose_120_ept4( double*       out,
                             const double* in,
                             int           np0,
                             int           np1,
                             int           np2 )
{

	__shared__ double tile[TILE_SIZE][TILE_SIZE + 1];

	int x_in, y_in, z,
	    x_out, y_out,
	    ind_in,
	    ind_out;

	int lx = threadIdx.x,
	    ly = threadIdx.y,
	    bx = blockIdx.x,
	    by = blockIdx.y;

	x_in = lx + TILE_SIZE * bx;
	y_in = ly + TILE_SIZE * by;

	z = blockIdx.z;

	x_out = ly + TILE_SIZE * bx;
	y_out = lx + TILE_SIZE * by;

	ind_in = x_in + (y_in + z * np1) * np0;
	ind_out = y_out + (z + x_out * np2) * np1;

	if( x_in < np0 && y_in < np1 )
	{
		tile[lx][ly] = in[ind_in];
		if( y_in + 4 < np1 )
		{
			tile[lx][ly +  4] = in[ind_in +  4*np0];
			if( y_in + 8 < np1 )
			{
				tile[lx][ly +  8] = in[ind_in +  8*np0];
				if( y_in + 12 < np1 )
				{
					tile[lx][ly +  12] = in[ind_in +  12*np0];
				}
			}
		}
	}

	__syncthreads();

	if( y_out < np1 && x_out < np0 )
	{
		out[ind_out] = tile[ly][lx];
		if( x_out + 4 < np0 )
		{
			out[ind_out +  4*np1*np2] = tile[ly + 4][lx];
			if( x_out + 8 < np0 )
			{
				out[ind_out +  8*np1*np2] = tile[ly + 8][lx];
				if( x_out + 12 < np0 )
				{
					out[ind_out +  12*np1*np2] = tile[ly + 12][lx];
				}
			}
		}
	}
}


__global__
void dev_transpose_120_in_place( double* data,
                                 int     np0 )
{
	__shared__ double cube13[BRICK_SIZE][BRICK_SIZE][BRICK_SIZE + 1];
	__shared__ double cube2[BRICK_SIZE][BRICK_SIZE][BRICK_SIZE + 1];

	int x1, y1, z1,
	    x2, y2, z2,
	    x3, y3, z3,
	    ind1, ind2, ind3;

	int lx = threadIdx.x,
	    ly = threadIdx.y,
	    lz = threadIdx.z,
	    bx = blockIdx.x,
	    by = blockIdx.y,
	    bz = blockIdx.z;

	int diagonal = (bx == by && by == bz);

	if( bx > by || bx > bz ||
	    ((bx == by || bx == bz) && by > bz) )
		return;

	x1 = lx + BRICK_SIZE * bx;
	y1 = ly + BRICK_SIZE * by;
	z1 = lz + BRICK_SIZE * bz;

	x2 = ly + BRICK_SIZE * bx;
	y2 = lx + BRICK_SIZE * by;
	z2 = lz + BRICK_SIZE * bz;

	x3 = lz + BRICK_SIZE * bx;
	y3 = ly + BRICK_SIZE * by;
	z3 = lx + BRICK_SIZE * bz;

	ind1 = x1 + (y1 + z1 * np0) * np0;
	ind2 = y2 + (z2 + x2 * np0) * np0;
	ind3 = z3 + (x3 + y3 * np0) * np0;

	// Swap lx and ly to avoid the synchronization commented below.
	if( x1 < np0 && y1 < np0 && z1 < np0 )
		cube13[ly][lx][lz] = data[ind1];
	if( x2 < np0 && y2 < np0 && z2 < np0 )
		if( ! diagonal )
			cube2[lx][ly][lz] = data[ind2];

	__syncthreads();

	if( diagonal )
	{
		if( x1 < np0 && y1 < np0 && z1 < np0 )
			data[ind1] = cube13[lx][lz][ly];
	}
	else
	{
		if( x2 < np0 && y2 < np0 && z2 < np0 )
			data[ind2] = cube13[lx][ly][lz];

		//__syncthreads();
		if( x3 < np0 && y3 < np0 && z3 < np0 )
		{
			cube13[lx][ly][lz] = data[ind3];
			data[ind3] = cube2[ly][lz][lx];
		}
		__syncthreads();

		if( x1 < np0 && y1 < np0 && z1 < np0 )
			data[ind1] = cube13[lz][ly][lx];
	}
}


/********************************************
 * Public functions                         *
 ********************************************/

__global__
void dev_transpose_201_ept1( double*       out,
                             const double* in,
                             int           np0,
                             int           np1,
                             int           np2 )
{
__shared__ double tile[TILE_SIZE][TILE_SIZE + 1];

	int x_in, y, z_in,
	    x_out, z_out,
	    ind_in,
	    ind_out;

	int lx = threadIdx.x,
	    ly = threadIdx.y,
	    bx = blockIdx.x,
	    by = blockIdx.y;

	x_in = lx + TILE_SIZE * bx;
	z_in = ly + TILE_SIZE * by;

	y = blockIdx.z;

	x_out = ly + TILE_SIZE * bx;
	z_out = lx + TILE_SIZE * by;


	ind_in = x_in + (y + z_in * np1) * np0;
	ind_out = z_out + (x_out + y * np0) * np2;

	if( x_in < np0 && z_in < np2 )
	{
		tile[lx][ly] = in[ind_in];
	}

	__syncthreads();

	if( x_out < np0 && z_out < np2 )
	{
		out[ind_out] = tile[ly][lx];
	}
}


__global__
void dev_transpose_201_ept2( double*       out,
                              const double* in,
                              int           np0,
                              int           np1,
                              int           np2 )
{
__shared__ double tile[TILE_SIZE][TILE_SIZE + 1];

	int x_in, y, z_in,
	    x_out, z_out,
	    ind_in,
	    ind_out;

	int lx = threadIdx.x,
	    ly = threadIdx.y,
	    bx = blockIdx.x,
	    by = blockIdx.y;

	x_in = lx + TILE_SIZE * bx;
	z_in = ly + TILE_SIZE * by;

	y = blockIdx.z;

	x_out = ly + TILE_SIZE * bx;
	z_out = lx + TILE_SIZE * by;


	ind_in = x_in + (y + z_in * np1) * np0;
	ind_out = z_out + (x_out + y * np0) * np2;

	if( x_in < np0 && z_in < np2 )
	{
		tile[lx][ly] = in[ind_in];
		if( z_in + 8 < np2 )
		{
			tile[lx][ly +  8] = in[ind_in +  8*np0*np1];
		}
	}

	__syncthreads();

	if( x_out < np0 && z_out < np2 )
	{
		out[ind_out] = tile[ly][lx];
		if( x_out + 8 < np0 )
		{
			out[ind_out +  8*np2] = tile[ly + 8][lx];
		}
	}
}


__global__
void dev_transpose_201_ept4( double*       out,
                             const double* in,
                             int           np0,
                             int           np1,
                             int           np2 )
{
	__shared__ double tile[TILE_SIZE][TILE_SIZE + 1];

	int x_in, y, z_in,
	    x_out, z_out,
	    ind_in,
	    ind_out;

	int lx = threadIdx.x,
	    ly = threadIdx.y,
	    bx = blockIdx.x,
	    by = blockIdx.y;

	x_in = lx + TILE_SIZE * bx;
	z_in = ly + TILE_SIZE * by;

	y = blockIdx.z;

	x_out = ly + TILE_SIZE * bx;
	z_out = lx + TILE_SIZE * by;


	ind_in = x_in + (y + z_in * np1) * np0;
	ind_out = z_out + (x_out + y * np0) * np2;

	if( x_in < np0 && z_in < np2 )
	{
		tile[lx][ly] = in[ind_in];
		if( z_in + 4 < np2 )
		{
			tile[lx][ly +  4] = in[ind_in +  4*np0*np1];
			if( z_in + 8 < np2 )
			{
				tile[lx][ly +  8] = in[ind_in +  8*np0*np1];
				if( z_in + 12 < np2 )
				{
					tile[lx][ly +  12] = in[ind_in +  12*np0*np1];
				}
			}
		}
	}

	__syncthreads();

	if( x_out < np0 && z_out < np2 )
	{
		out[ind_out] = tile[ly][lx];
		if( x_out + 4 < np0 )
		{
			out[ind_out +  4*np2] = tile[ly + 4][lx];
			if( x_out + 8 < np0 )
			{
				out[ind_out +  8*np2] = tile[ly + 8][lx];
				if( x_out + 12 < np0 )
				{
					out[ind_out +  12*np2] = tile[ly + 12][lx];
				}
			}
		}
	}
}


__global__
void dev_transpose_201_in_place( double* data,
                                 int     np0 )
{
	__shared__ double cube12[BRICK_SIZE][BRICK_SIZE][BRICK_SIZE + 1];
	__shared__ double cube3[BRICK_SIZE][BRICK_SIZE][BRICK_SIZE + 1];

	int x1, y1, z1,
	    x2, y2, z2,
	    x3, y3, z3,
	    ind1, ind2, ind3;

	int lx = threadIdx.x,
	    ly = threadIdx.y,
	    lz = threadIdx.z,
	    bx = blockIdx.x,
	    by = blockIdx.y,
	    bz = blockIdx.z;

	int diagonal = (bx == by && by == bz);

	if( bx > by || bx > bz ||
	    ((bx == by || bx == bz) && by > bz) )
		return;

	x1 = lx + BRICK_SIZE * bx;
	y1 = ly + BRICK_SIZE * by;
	z1 = lz + BRICK_SIZE * bz;

	x2 = ly + BRICK_SIZE * bx;
	y2 = lx + BRICK_SIZE * by;
	z2 = lz + BRICK_SIZE * bz;

	x3 = lz + BRICK_SIZE * bx;
	y3 = ly + BRICK_SIZE * by;
	z3 = lx + BRICK_SIZE * bz;

	ind1 = x1 + (y1 + z1 * np0) * np0;
	ind2 = y2 + (z2 + x2 * np0) * np0;
	ind3 = z3 + (x3 + y3 * np0) * np0;

	// Swap lx and lz to avoid the synchronization commented below.
	if( x1 < np0 && y1 < np0 && z1 < np0 )
		cube12[lz][ly][lx] = data[ind1];
	if( x3 < np0 && y3 < np0 && z3 < np0 )
		if( ! diagonal )
			cube3[lx][ly][lz] = data[ind3];

	__syncthreads();

	if( diagonal )
	{
		if( x1 < np0 && y1 < np0 && z1 < np0 )
			data[ind1] = cube12[lx][lz][ly];
	}
	else
	{
		if( x3 < np0 && y3 < np0 && z3 < np0 )
			data[ind3] = cube12[lx][ly][lz];

		//__syncthreads();
		if( x2 < np0 && y2 < np0 && z2 < np0 )
		{
			cube12[lx][ly][lz] = data[ind2];
			data[ind2] = cube3[lz][lx][ly];
		}
		__syncthreads();

		if( x1 < np0 && y1 < np0 && z1 < np0 )
			data[ind1] = cube12[ly][lx][lz];
	}
}

 
/********************************************
 * Public functions                         *
 ********************************************/

__global__
void dev_transpose_210_ept1( double*       out,
                             const double* in,
                             int           np0,
                             int           np1,
                             int           np2 )
{

	__shared__ double tile[TILE_SIZE][TILE_SIZE + 1];
	
	int x_in, y, z_in,
	    x_out, z_out,
	    ind_in,
	    ind_out;
	
	int lx = threadIdx.x,
	    ly = threadIdx.y,
	    bx = blockIdx.x,
	    by = blockIdx.y;

	x_in = lx + TILE_SIZE * bx;
	z_in = ly + TILE_SIZE * by;

	y = blockIdx.z;

	x_out = ly + TILE_SIZE * bx;
	z_out = lx + TILE_SIZE * by;


	ind_in = x_in + (y + z_in * np1) * np0;
	ind_out = z_out + (y + x_out * np1) * np2;

	if( x_in < np0 && z_in < np2 )
	{
			tile[lx][ly] = in[ind_in];
	}	

	__syncthreads();

	if( z_out < np2 && x_out < np0 )
	{
		out[ind_out] = tile[ly][lx];
	}

}


__global__
void dev_transpose_210_ept2( double*       out,
                             const double* in,
                             int           np0,
                             int           np1,
                             int           np2 )
{

	__shared__ double tile[TILE_SIZE][TILE_SIZE + 1];
	
	int x_in, y, z_in,
	    x_out, z_out,
	    ind_in,
	    ind_out;
	
	int lx = threadIdx.x,
	    ly = threadIdx.y,
	    bx = blockIdx.x,
	    by = blockIdx.y;
		
	x_in = lx + TILE_SIZE * bx;
	z_in = ly + TILE_SIZE * by;

	y = blockIdx.z;

	x_out = ly + TILE_SIZE * bx;
	z_out = lx + TILE_SIZE * by;

	ind_in = x_in + (y + z_in * np1) * np0;
	ind_out = z_out + (y + x_out * np1) * np2;

	if( x_in < np0 && z_in < np2 )
	{
		tile[lx][ly] = in[ind_in];
		if( z_in + 8 < np2 )
		{
			tile[lx][ly +  8] = in[ind_in +  8*np0*np1];
		}
	}
	
	__syncthreads();

	if( z_out < np2 && x_out < np0 )
	{
		out[ind_out] = tile[ly][lx];
		if( x_out + 8 < np0 )
		{
			out[ind_out +  8*np2*np1] = tile[ly + 8][lx];
		}
	}

}


__global__
void dev_transpose_210_ept4( double*       out,
                             const double* in,
                             int           np0,
                             int           np1,
                             int           np2 )
{

	__shared__ double tile[TILE_SIZE][TILE_SIZE + 1];
	
	int x_in, y, z_in,
	    x_out, z_out,
	    ind_in,
	    ind_out;
	
	int lx = threadIdx.x,
	    ly = threadIdx.y,
	    bx = blockIdx.x,
	    by = blockIdx.y;
		
	x_in = lx + TILE_SIZE * bx;
	z_in = ly + TILE_SIZE * by;

	y = blockIdx.z;

	x_out = ly + TILE_SIZE * bx;
	z_out = lx + TILE_SIZE * by;

	ind_in = x_in + (y + z_in * np1) * np0;
	ind_out = z_out + (y + x_out * np1) * np2;

	if( x_in < np0 && z_in < np2 )
	{
		tile[lx][ly] = in[ind_in];
		if( z_in + 4 < np2 )
		{
			tile[lx][ly +  4] = in[ind_in +  4*np0*np1];
			if( z_in + 8 < np2 )
			{
				tile[lx][ly +  8] = in[ind_in +  8*np0*np1];
				if( z_in + 12 < np2 )
				{
					tile[lx][ly + 12] = in[ind_in + 12*np0*np1];
				}
			}
		}
	}
	
	__syncthreads();
	
	if( z_out < np2 && x_out < np0 )
	{
		out[ind_out] = tile[ly][lx];
		if( x_out + 4 < np0 )
		{
			out[ind_out +  4*np2*np1] = tile[ly +  4][lx];
			if( x_out + 8 < np0 )
			{
				out[ind_out +  8*np2*np1] = tile[ly +  8][lx];
				if( x_out + 12 < np0 )
				{
					out[ind_out + 12*np2*np1] = tile[ly + 12][lx];
				}
			}
		}
	}
}


__global__
void dev_transpose_210_in_place( double* data,
                                 int     np0,
                                 int     np1 )
{
	__shared__ double inf_tile[TILE_SIZE][TILE_SIZE + 1];
	__shared__ double sup_tile[TILE_SIZE][TILE_SIZE + 1];
	
	int x_inf, y, z_inf,
	    x_sup, z_sup,
	    ind_inf,
	    ind_sup;

	int lx = threadIdx.x,
	    ly = threadIdx.y,
	    bx = blockIdx.x,
	    by = blockIdx.y;

	if( bx > by ) // Block in upper triangle
		return;

	x_inf = lx + TILE_SIZE * bx;
	z_inf = ly + TILE_SIZE * by;
	y = blockIdx.z;

	if( x_inf < np0 && z_inf < np0 )
	{
		ind_inf = x_inf + (y + z_inf * np1) * np0;
		inf_tile[lx][ly] = data[ind_inf];
	}

	x_sup = ly + TILE_SIZE * bx;
	z_sup = lx + TILE_SIZE * by;
	if( bx < by ) // Block in lower triangle
	{
		if( x_sup < np0 && z_sup < np0 )
		{
			ind_sup = z_sup + (y + x_sup * np1) * np0;
			sup_tile[lx][ly] = data[ind_sup];
		}
	}
	else // Block in diagonal
		ind_sup = ind_inf;

	__syncthreads();

	if( x_sup < np0 && z_sup < np0 )
	{
		data[ind_sup] = inf_tile[ly][lx];
	}

	if( bx < by )
	{
		if( x_inf < np0 && z_inf < np0 )
		{
			data[ind_inf] = sup_tile[ly][lx];
		}
	}
}
/*  This file is part of cuTranspose.

    cuTranspose is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    cuTranspose is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with cuTranspose.  If not, see <http://www.gnu.org/licenses/>.

    Copyright 2016 Ibai Gurrutxaga, Javier Muguerza, Jose L. Jodra.
*/

int cut_transpose3d( double*       output,
                     const double* input,
                     const int*    size,
                     const int*    permutation,
                     int           elements_per_thread,
                     gpuStream_t stream)

{
	dim3 num_blocks,
	     block_size;
	int in_place = output == input;

	if( !valid_parameters( in_place, size, permutation, elements_per_thread ) )
		return -1;

	if(( permutation[0] == 1 && permutation[1] == 2 && permutation[2] == 0 && in_place ) ||
	   ( permutation[0] == 2 && permutation[1] == 0 && permutation[2] == 1 && in_place ))
		set_grid_dims_cube( size,
		                    &block_size,
		                    &num_blocks,
		                    elements_per_thread );
	else
		set_grid_dims( size,
		               permutation[0],
		               &block_size,
		               &num_blocks,
		               elements_per_thread );

	if( permutation[0] == 0 && permutation[1] == 1 && permutation[2] == 2 )
	{
			dev_copy<<< num_blocks, block_size >>>( output,
			                                        input,
			                                        size[0],
			                                        size[1], 
			                                        elements_per_thread );
	}
	else if( permutation[0] == 0 && permutation[1] == 2 && permutation[2] == 1 )
	{
		if( in_place )
		{
			if( size[1] != size[2] )
			{
				fprintf( stderr, "This in-place transposition requires equal dimensions for "
								 "the second and third dimensions.\n" );
				return -1;
			}
			dev_transpose_021_in_place<<< num_blocks, block_size >>>( output,
			                                                          size[0],
			                                                          size[1] );
		}
		else {
			switch(elements_per_thread)
			{
			case 1:
				dev_transpose_021_ept1<<< num_blocks, block_size >>>( output,
				                                                      input,
				                                                      size[0],
				                                                      size[1],
				                                                      size[2]);
				break;
			case 2:
				dev_transpose_021_ept2<<< num_blocks, block_size >>>( output,
				                                                      input,
				                                                      size[0],
				                                                      size[1],
				                                                      size[2]);
				break;
			case 4:
				dev_transpose_021_ept4<<< num_blocks, block_size >>>( output,
				                                                      input,
				                                                      size[0],
				                                                      size[1],
				                                                      size[2]);
				break;
			}
		}
	}
	else if( permutation[0] == 1 && permutation[1] == 0 && permutation[2] == 2 )
	{
		if( in_place )
		{
			if( size[0] != size[1] )
			{
				fprintf( stderr, "The 102 in-place transposition requires equal dimensions for "
								 "the first and second dimensions.\n" );
				return -1;
			}
			dev_transpose_102_in_place<<< num_blocks, block_size >>>( output,
                                                                      size[0],
                                                                      size[2] );
		}
		else {
			switch(elements_per_thread)
			{
			case 1:
				dev_transpose_102_ept1<<< num_blocks, block_size >>>( output,
				                                                      input,
				                                                      size[0],
				                                                      size[1],
				                                                      size[2]);
				break;
			case 2:
				dev_transpose_102_ept2<<< num_blocks, block_size >>>( output,
				                                                      input,
				                                                      size[0],
				                                                      size[1],
				                                                      size[2]);
				break;
			case 4:
				dev_transpose_102_ept4<<< num_blocks, block_size >>>( output,
				                                                      input,
				                                                      size[0],
				                                                      size[1],
				                                                      size[2]);
				break;
			}
		}
	}
	else if( permutation[0] == 1 && permutation[1] == 2 && permutation[2] == 0 )
	{
		if( in_place )
		{
			if( size[0] != size[1] || size[0] != size[2] )
			{
				fprintf( stderr, "The 120 in-place transposition requires a cubic input.\n" );
				return -1;
			}
			dev_transpose_120_in_place<<< num_blocks, block_size >>>( output,
			                                                          size[0] );
		}
		else {
			switch(elements_per_thread)
			{
			case 1:
				dev_transpose_120_ept1<<< num_blocks, block_size >>>( output,
				                                                      input,
				                                                      size[0],
				                                                      size[1],
				                                                      size[2]);
				break;
			case 2:
				dev_transpose_120_ept2<<< num_blocks, block_size >>>( output,
				                                                      input,
				                                                      size[0],
				                                                      size[1],
				                                                      size[2]);
				break;
			case 4:
				dev_transpose_120_ept4<<< num_blocks, block_size >>>( output,
				                                                      input,
				                                                      size[0],
				                                                      size[1],
				                                                      size[2]);
				break;
			}
		}
	}
	else if( permutation[0] == 2 && permutation[1] == 0 && permutation[2] == 1 )
	{
		if( in_place )
		{
			if( size[0] != size[1] || size[0] != size[2] )
			{
				fprintf( stderr, "The 201 in-place transposition requires a cubic input.\n" );
				return -1;
			}
			dev_transpose_201_in_place<<< num_blocks, block_size >>>( output,
			                                                          size[0] );
		}
		else {
			switch(elements_per_thread)
			{
			case 1:
				dev_transpose_201_ept1<<< num_blocks, block_size >>>( output,
				                                                      input,
				                                                      size[0],
				                                                      size[1],
				                                                      size[2]);
				break;
			case 2:
				dev_transpose_201_ept2<<< num_blocks, block_size >>>( output,
				                                                      input,
				                                                      size[0],
				                                                      size[1],
				                                                      size[2]);
				break;
			case 4:
				dev_transpose_201_ept4<<< num_blocks, block_size >>>( output,
				                                                      input,
				                                                      size[0],
				                                                      size[1],
				                                                      size[2]);
				break;
			}
		}
	}
	else if( permutation[0] == 2 && permutation[1] == 1 && permutation[2] == 0 )
	{
		if( in_place )
		{
			if( size[0] != size[2] )
			{
				fprintf( stderr, "The 210 in-place transposition requires equal dimensions for "
								 "the first and third dimensions.\n" );
				return -1;
			}
			dev_transpose_210_in_place<<< num_blocks, block_size >>>( output,
			                                                          size[0],
			                                                          size[1] );
		}
		else {
			switch(elements_per_thread)
			{
			case 1:
				dev_transpose_210_ept1<<< num_blocks, block_size >>>( output,
				                                                      input,
				                                                      size[0],
				                                                      size[1],
				                                                      size[2]);
				break;
			case 2:
				dev_transpose_210_ept2<<< num_blocks, block_size >>>( output,
				                                                      input,
				                                                      size[0],
				                                                      size[1],
				                                                      size[2]);
				break;
			case 4:
				dev_transpose_210_ept4<<< num_blocks, block_size >>>( output,
				                                                      input,
				                                                      size[0],
				                                                      size[1],
				                                                      size[2]);
				break;
			}
		}
	}
	
	return 0;
}

