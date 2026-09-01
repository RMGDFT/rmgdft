/*
 *
 * Copyright (c) 2013, Emil Briggs
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the <organization> nor the
 *       names of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.

 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * 
 * @file
 *
 *
 * @section DESCRIPTION
 * Used to access grid and node related data.
 */

#include "rmg_grid.h"

namespace rmg
{
    /// Used to set up global coarse grid dimensions, MPI node dimensions and the ratio of the fine grid to the coarse in each coordinate dimension.
    /// @param newNX_GRID New global coarse grid X dimension
    /// @param newNY_GRID New global coarse grid Y dimension
    /// @param newNZ_GRID New global coarse grid Z dimension
    /// @param newPE_X New MPI grid X dimension
    /// @param newPE_Y New MPI grid Y dimension
    /// @param newPE_Z New MPI grid Z dimension
    /// @param newFG_RATIO New ratio of fine grid to coarse
    grid::grid(int newNX_GRID, int newNY_GRID, int newNZ_GRID, int newPE_X, int newPE_Y, int newPE_Z, int new_rank, int newFG_RATIO)
    {

        grid::NX_GRID = newNX_GRID;
        grid::NY_GRID = newNY_GRID;
        grid::NZ_GRID = newNZ_GRID;

        grid::PE_X = newPE_X;
        grid::PE_Y = newPE_Y;
        grid::PE_Z = newPE_Z;
        grid::NPES = newPE_X * newPE_Y * newPE_Z;

        grid::default_FG_RATIO = newFG_RATIO;

//        grid::set_rank(new_rank);
    }

    void grid::set_rank(int new_rank, MPI_Comm comm)
    {
        int rem;

        grid::rank = new_rank;
        grid::comm = comm;

        /// Have each node figure out who it's neighbors are
        int ii, jj, kk;
        grid::pe2xyz (grid::rank, &ii, &jj, &kk);
        neighbors[NB_N] = grid::xyz2pe (ii, (jj + 1) % grid::PE_Y, kk);
        neighbors[NB_S] = grid::xyz2pe (ii, (jj - 1 + grid::PE_Y) % get_PE_Y(), kk);
        neighbors[NB_E] = grid::xyz2pe ((ii + 1) % grid::PE_X, jj, kk);
        neighbors[NB_W] = grid::xyz2pe ((ii - 1 + grid::PE_X) % grid::PE_X, jj, kk);
        neighbors[NB_U] = grid::xyz2pe (ii, jj, (kk + 1) % grid::PE_Z);
        neighbors[NB_D] = grid::xyz2pe (ii, jj, (kk - 1 + grid::PE_Z) % grid::PE_Z);

        // Compute grid sizes for each node.
        grid::PX0_GRID = grid::NX_GRID / grid::PE_X;
        rem = grid::NX_GRID % grid::PE_X;
        if(rem && (ii < rem)) grid::PX0_GRID++;

        grid::PY0_GRID = grid::NY_GRID / grid::PE_Y;
        rem = grid::NY_GRID % grid::PE_Y;
        if(rem && (jj < rem)) grid::PY0_GRID++;

        grid::PZ0_GRID = grid::NZ_GRID / grid::PE_Z;
        rem = grid::NZ_GRID % grid::PE_Z;
        if(rem && (kk < rem)) grid::PZ0_GRID++;

        // Adjust if needed
        grid::find_node_sizes(rank, grid::NX_GRID, grid::NY_GRID, grid::NZ_GRID, &this->PX0_GRID, &this->PY0_GRID, &this->PZ0_GRID);

        grid::P0_BASIS = grid::PX0_GRID * grid::PY0_GRID * grid::PZ0_GRID;

        // Now compute the global grid offset of the first point of the coarse and fine node grids
        grid::find_node_offsets(rank, grid::NX_GRID, grid::NY_GRID, grid::NZ_GRID,
                          &this->PX_OFFSET, &this->PY_OFFSET, &this->PZ_OFFSET);

                 

    }

    void grid::find_node_sizes(int rank, int nxgrid, int nygrid, int nzgrid, int *pxsize, int *pysize, int *pzsize)
    {

	int ii, jj, kk;
	int ix, iy, iz, mfac;

	grid::pe2xyz (rank, &ii, &jj, &kk);

	mfac = nxgrid / grid::NX_GRID;
	*pxsize = mfac * (grid::NX_GRID / grid::PE_X);
	ix = grid::NX_GRID % grid::PE_X;
	if(ii < ix) *pxsize += mfac;

	mfac = nygrid / grid::NY_GRID;
	*pysize = mfac * (grid::NY_GRID / grid::PE_Y);
	iy = grid::NY_GRID % grid::PE_Y;
	if(jj < iy) *pysize += mfac;

	mfac = nzgrid / grid::NZ_GRID;
	*pzsize = mfac * (grid::NZ_GRID / grid::PE_Z);
	iz = grid::NZ_GRID % grid::PE_Z;
	if(kk < iz) *pzsize += mfac;
    }

    void grid::find_node_offsets(int rank, int nxgrid, int nygrid, int nzgrid, int *pxoffset, int *pyoffset, int *pzoffset)
    {

	int ii, jj, kk;
	int idx, ix, iy, iz, ioffset, mfac;

	grid::pe2xyz (rank, &ii, &jj, &kk);

	// Now compute the global grid offset of the first point of the node grid
	mfac = nxgrid / grid::NX_GRID;
	*pxoffset = mfac * ii*(grid::NX_GRID / grid::PE_X);
        ix = grid::NX_GRID % grid::PE_X;
	ioffset = 0;
	for(idx = 1;idx <= ii;idx++) {
	    if(idx <= ix) ioffset++;
	}
	ioffset *= mfac;
	*pxoffset = *pxoffset + ioffset;


	mfac = nygrid / grid::NY_GRID;
	*pyoffset = mfac * jj*(grid::NY_GRID / grid::PE_Y);
        iy = grid::NY_GRID % grid::PE_Y;
	ioffset = 0;
	for(idx = 1;idx <= jj;idx++) {
	    if(idx <= iy) ioffset++;
	}
	ioffset *= mfac;
	*pyoffset += ioffset;

	mfac = nzgrid / grid::NZ_GRID;
	*pzoffset = mfac * kk*(grid::NZ_GRID / grid::PE_Z);
        iz = grid::NZ_GRID % grid::PE_Z;
	ioffset = 0;
	for(idx = 1;idx <= kk;idx++) {
	    if(idx <= iz) ioffset++;
	}
	ioffset *= mfac;
	*pzoffset += ioffset;

    }

    int grid::get_default_FG_RATIO(void)
    {
        return grid::default_FG_RATIO;
    }
    int grid::get_NPES(void)
    {
	return grid::NPES;
    }
    int grid::get_PE_X(void)
    {
	return grid::PE_X;
    }
    int grid::get_PE_Y(void)
    {
	return grid::PE_Y;
    }
    int grid::get_PE_Z(void)
    {
	return grid::PE_Z;
    }

    int grid::get_NX_GRID(int density)
    {
	return density * grid::NX_GRID;
    }
    int grid::get_NY_GRID(int density)
    {
	return density * grid::NY_GRID;
    }
    int grid::get_NZ_GRID(int density)
    {
	return density * grid::NZ_GRID;
    }

    double grid::get_hxgrid(int density)
    {
        return 1.0 / ((double)(density * grid::NX_GRID));
    }
    double grid::get_hygrid(int density)
    {
        return 1.0 / ((double)(density * grid::NY_GRID));
    }
    double grid::get_hzgrid(int density)
    {
        return 1.0 / ((double)(density * grid::NZ_GRID));
    }


    int grid::get_PX0_GRID(int density)
    {
	return density * grid::PX0_GRID;
    }
    int grid::get_PY0_GRID(int density)
    {
	return density * grid::PY0_GRID;
    }
    int grid::get_PZ0_GRID(int density)
    {
	return density * grid::PZ0_GRID;
    }
    int grid::get_PX_OFFSET(int density)
    {
         return density * grid::PX_OFFSET;
    }
    int grid::get_PY_OFFSET(int density)
    {
	return density * grid::PY_OFFSET;
    }
    int grid::get_PZ_OFFSET(int density)
    {
	return density * grid::PZ_OFFSET;
    }
    size_t grid::get_P0_BASIS(int density)
    {
	return (size_t)(density * density * density) * grid::P0_BASIS;
    }
    size_t grid::get_GLOBAL_BASIS(int density)
    {
	return (size_t)(density * density * density) * (size_t)grid::NX_GRID * (size_t)grid::NY_GRID * (size_t)grid::NZ_GRID;
    }
    void grid::pe2xyz(int pe, int *x, int *y, int *z)
    {

        *x = pe;
        *z = *x % grid::PE_Z;
        *x /= grid::PE_Z;
        *y = *x % grid::PE_Y;
        *x /= grid::PE_Y;

        if (*x >= grid::PE_X)
            *x -= grid::PE_X;
        if (*x >= grid::PE_X)
            *x -= grid::PE_X;

    }                             

    int grid::xyz2pe(int x, int y, int z)
    {
        return  x * grid::PE_Y * grid::PE_Z + y * grid::PE_Z + z;
    }

    // Returns a pointer to the neighbors structure which contains the rank
    // of neighboring processors in three-dimensional space.
    int *grid::get_neighbors(void)
    {
	return grid::neighbors;
    }

    int grid::get_rank(void)
    {
        return grid::rank;
    }

    bool grid::is_face_pe(void)
    {
        return grid::face_pe;
    }

    bool grid::is_edge_pe(void)
    {
        return grid::edge_pe;
    }

    bool grid::is_corner_pe(void)
    {
        return grid::corner_pe;
    }
}
