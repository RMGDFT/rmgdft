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

#include "BaseGrid.h"


using namespace std;



    /// Used to set up global coarse grid dimensions, MPI node dimensions and the ratio of the fine grid to the coarse in each coordinate dimension.
    /// @param newNX_GRID New global coarse grid X dimension
    /// @param newNY_GRID New global coarse grid Y dimension
    /// @param newNZ_GRID New global coarse grid Z dimension
    /// @param newPE_X New MPI grid X dimension
    /// @param newPE_Y New MPI grid Y dimension
    /// @param newPE_Z New MPI grid Z dimension
    /// @param newFG_RATIO New ratio of fine grid to coarse
    BaseGrid::BaseGrid(int newNX_GRID, int newNY_GRID, int newNZ_GRID, int newPE_X, int newPE_Y, int newPE_Z, int newFG_RATIO)
    {

        BaseGrid::NX_GRID = newNX_GRID;
        BaseGrid::NY_GRID = newNY_GRID;
        BaseGrid::NZ_GRID = newNZ_GRID;

        BaseGrid::PE_X = newPE_X;
        BaseGrid::PE_Y = newPE_Y;
        BaseGrid::PE_Z = newPE_Z;

        BaseGrid::default_FG_RATIO = newFG_RATIO;
        BaseGrid::neighbor_first = 0;
    }

    void BaseGrid::set_nodes(int newgridpe)
    {
        int rem;

        BaseGrid::gridpe = newgridpe;

        /// Have each node figure out who it's neighbors are
        int ii, jj, kk;
        BaseGrid::pe2xyz (BaseGrid::gridpe, &ii, &jj, &kk);
        neighbors[NB_N] = BaseGrid::xyz2pe (ii, (jj + 1) % BaseGrid::PE_Y, kk);
        neighbors[NB_S] = BaseGrid::xyz2pe (ii, (jj - 1 + BaseGrid::PE_Y) % get_PE_Y(), kk);
        neighbors[NB_E] = BaseGrid::xyz2pe ((ii + 1) % BaseGrid::PE_X, jj, kk);
        neighbors[NB_W] = BaseGrid::xyz2pe ((ii - 1 + BaseGrid::PE_X) % BaseGrid::PE_X, jj, kk);
        neighbors[NB_U] = BaseGrid::xyz2pe (ii, jj, (kk + 1) % BaseGrid::PE_Z);
        neighbors[NB_D] = BaseGrid::xyz2pe (ii, jj, (kk - 1 + BaseGrid::PE_Z) % BaseGrid::PE_Z);
        BaseGrid::neighbor_first = 1;

        // Compute grid sizes for each node.
        BaseGrid::PX0_GRID = BaseGrid::NX_GRID / BaseGrid::PE_X;
        rem = BaseGrid::NX_GRID % BaseGrid::PE_X;
        if(rem && (ii < rem)) BaseGrid::PX0_GRID++;

        BaseGrid::PY0_GRID = BaseGrid::NY_GRID / BaseGrid::PE_Y;
        rem = BaseGrid::NY_GRID % BaseGrid::PE_Y;
        if(rem && (jj < rem)) BaseGrid::PY0_GRID++;

        BaseGrid::PZ0_GRID = BaseGrid::NZ_GRID / BaseGrid::PE_Z;
        rem = BaseGrid::NZ_GRID % BaseGrid::PE_Z;
        if(rem && (kk < rem)) BaseGrid::PZ0_GRID++;

        // Adjust if needed
        BaseGrid::find_node_sizes(gridpe, BaseGrid::NX_GRID, BaseGrid::NY_GRID, BaseGrid::NZ_GRID, &this->PX0_GRID, &this->PY0_GRID, &this->PZ0_GRID);

        BaseGrid::P0_BASIS = BaseGrid::PX0_GRID * BaseGrid::PY0_GRID * BaseGrid::PZ0_GRID;

        // Now compute the global grid offset of the first point of the coarse and fine node grids
        BaseGrid::find_node_offsets(gridpe, BaseGrid::NX_GRID, BaseGrid::NY_GRID, BaseGrid::NZ_GRID,
                          &this->PX_OFFSET, &this->PY_OFFSET, &this->PZ_OFFSET);

                 

    }

    void BaseGrid::find_node_sizes(int gridpe, int nxgrid, int nygrid, int nzgrid, int *pxsize, int *pysize, int *pzsize)
    {

	int ii, jj, kk;
	int ix, iy, iz, mfac;

	BaseGrid::pe2xyz (gridpe, &ii, &jj, &kk);

	mfac = nxgrid / BaseGrid::NX_GRID;
	*pxsize = mfac * (BaseGrid::NX_GRID / BaseGrid::PE_X);
	ix = BaseGrid::NX_GRID % BaseGrid::PE_X;
	if(ii < ix) *pxsize += mfac;

	mfac = nygrid / BaseGrid::NY_GRID;
	*pysize = mfac * (BaseGrid::NY_GRID / BaseGrid::PE_Y);
	iy = BaseGrid::NY_GRID % BaseGrid::PE_Y;
	if(jj < iy) *pysize += mfac;

	mfac = nzgrid / BaseGrid::NZ_GRID;
	*pzsize = mfac * (BaseGrid::NZ_GRID / BaseGrid::PE_Z);
	iz = BaseGrid::NZ_GRID % BaseGrid::PE_Z;
	if(kk < iz) *pzsize += mfac;
    }

    void BaseGrid::find_node_offsets(int gridpe, int nxgrid, int nygrid, int nzgrid, int *pxoffset, int *pyoffset, int *pzoffset)
    {

	int ii, jj, kk;
	int idx, ix, iy, iz, ioffset, mfac;

	BaseGrid::pe2xyz (gridpe, &ii, &jj, &kk);

	// Now compute the global grid offset of the first point of the node grid
	mfac = nxgrid / BaseGrid::NX_GRID;
	*pxoffset = mfac * ii*(BaseGrid::NX_GRID / BaseGrid::PE_X);
        ix = BaseGrid::NX_GRID % BaseGrid::PE_X;
	ioffset = 0;
	for(idx = 1;idx <= ii;idx++) {
	    if(idx <= ix) ioffset++;
	}
	ioffset *= mfac;
	*pxoffset = *pxoffset + ioffset;


	mfac = nygrid / BaseGrid::NY_GRID;
	*pyoffset = mfac * jj*(BaseGrid::NY_GRID / BaseGrid::PE_Y);
        iy = BaseGrid::NY_GRID % BaseGrid::PE_Y;
	ioffset = 0;
	for(idx = 1;idx <= jj;idx++) {
	    if(idx <= iy) ioffset++;
	}
	ioffset *= mfac;
	*pyoffset += ioffset;

	mfac = nzgrid / BaseGrid::NZ_GRID;
	*pzoffset = mfac * kk*(BaseGrid::NZ_GRID / BaseGrid::PE_Z);
        iz = BaseGrid::NZ_GRID % BaseGrid::PE_Z;
	ioffset = 0;
	for(idx = 1;idx <= kk;idx++) {
	    if(idx <= iz) ioffset++;
	}
	ioffset *= mfac;
	*pzoffset += ioffset;

    }

    int BaseGrid::get_default_FG_RATIO(void)
    {
        return BaseGrid::default_FG_RATIO;
    }
    int BaseGrid::get_PE_X(void)
    {
	return BaseGrid::PE_X;
    }
    int BaseGrid::get_PE_Y(void)
    {
	return BaseGrid::PE_Y;
    }
    int BaseGrid::get_PE_Z(void)
    {
	return BaseGrid::PE_Z;
    }

    int BaseGrid::get_NX_GRID(int density)
    {
	return density * BaseGrid::NX_GRID;
    }
    int BaseGrid::get_NY_GRID(int density)
    {
	return density * BaseGrid::NY_GRID;
    }
    int BaseGrid::get_NZ_GRID(int density)
    {
	return density * BaseGrid::NZ_GRID;
    }

    double BaseGrid::get_hxgrid(int density)
    {
        return 1.0 / ((double)(density * BaseGrid::NX_GRID));
    }
    double BaseGrid::get_hygrid(int density)
    {
        return 1.0 / ((double)(density * BaseGrid::NY_GRID));
    }
    double BaseGrid::get_hzgrid(int density)
    {
        return 1.0 / ((double)(density * BaseGrid::NZ_GRID));
    }


    int BaseGrid::get_PX0_GRID(int density)
    {
	return density * BaseGrid::PX0_GRID;
    }
    int BaseGrid::get_PY0_GRID(int density)
    {
	return density * BaseGrid::PY0_GRID;
    }
    int BaseGrid::get_PZ0_GRID(int density)
    {
	return density * BaseGrid::PZ0_GRID;
    }
    int BaseGrid::get_PX_OFFSET(int density)
    {
	return BaseGrid::PX_OFFSET;
    }
    int BaseGrid::get_PY_OFFSET(int density)
    {
	return BaseGrid::PY_OFFSET;
    }
    int BaseGrid::get_PZ_OFFSET(int density)
    {
	return BaseGrid::PZ_OFFSET;
    }
    int BaseGrid::get_P0_BASIS(int density)
    {
	return density * density * density * BaseGrid::P0_BASIS;
    }
    int BaseGrid::get_GLOBAL_BASIS(int density)
    {
	return density * density * density * BaseGrid::NX_GRID * BaseGrid::NY_GRID * BaseGrid::NZ_GRID;
    }
    void BaseGrid::pe2xyz(int pe, int *x, int *y, int *z)
    {

        *x = pe;
        *z = *x % BaseGrid::PE_Z;
        *x /= BaseGrid::PE_Z;
        *y = *x % BaseGrid::PE_Y;
        *x /= BaseGrid::PE_Y;

        if (*x >= BaseGrid::PE_X)
            *x -= BaseGrid::PE_X;
        if (*x >= BaseGrid::PE_X)
            *x -= BaseGrid::PE_X;

    }                             

    int BaseGrid::xyz2pe(int x, int y, int z)
    {
        return  x * BaseGrid::PE_Y * BaseGrid::PE_Z + y * BaseGrid::PE_Z + z;
    }

    // Returns a pointer to the neighbors structure which contains the rank
    // of neighboring processors in three-dimensional space.
    int *BaseGrid::get_neighbors(void)
    {
	if(!BaseGrid::neighbor_first)
	    rmg_error_handler (__FILE__, __LINE__, "Neighbor list not initialized. Please call set_neighbors first");

	return BaseGrid::neighbors;
    }

    int BaseGrid::get_gridpe(void)
    {
	if(!BaseGrid::neighbor_first)
	    rmg_error_handler (__FILE__, __LINE__, "Neighbor list not initialized. Please call set_neighbors first");
        return BaseGrid::gridpe;
    }


