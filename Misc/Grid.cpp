/* This C++ class is currently tracking with Misc/grid.c until the C to C++ conversion is further along */

#include <iostream>
#include <cstdio>
#include "rmgtypes.h"
#include "grid.h"
#include "Grid.h"
#include "rmg_error.h"

extern "C" {
int find_node_offsets(int gridpe, int nxgrid, int nygrid, int nzgrid,
                      int *pxoffset, int *pyoffset, int *pzoffset);
int find_node_sizes(int gridpe, int nxgrid, int nygrid, int nzgrid,
                      int *pxsize, int *pysize, int *pzsize);
}

using namespace std;



    void Grid::set_grids(int newgridpe, int ii, int jj, int kk)
    {
        int rem;

        if(grid_first) return;

        gridpe = newgridpe;

        // Compute grid sizes for each node.
        Grid::PX0_GRID = NX_GRID / PE_X;
        rem = NX_GRID % PE_X;
        if(rem && (ii < rem)) Grid::PX0_GRID++;

        Grid::PY0_GRID = NY_GRID / PE_Y;
        rem = NY_GRID % PE_Y;
        if(rem && (jj < rem)) Grid::PY0_GRID++;

        Grid::PZ0_GRID = NZ_GRID / PE_Z;
        rem = NZ_GRID % PE_Z;
        if(rem && (kk < rem)) Grid::PZ0_GRID++;

        find_node_sizes(gridpe, NX_GRID, NY_GRID, NZ_GRID, &Grid::PX0_GRID, &Grid::PY0_GRID, &Grid::PZ0_GRID);
        find_node_sizes(gridpe, FNX_GRID, FNY_GRID, FNZ_GRID, &Grid::FPX0_GRID, &Grid::FPY0_GRID, &Grid::FPZ0_GRID);

        Grid::P0_BASIS = Grid::PX0_GRID * Grid::PY0_GRID * Grid::PZ0_GRID;
        Grid::FP0_BASIS = Grid::FPX0_GRID * Grid::FPY0_GRID * Grid::FPZ0_GRID;

        // Now compute the global grid offset of the first point of the coarse and fine node grids
        find_node_offsets(gridpe, NX_GRID, NY_GRID, NZ_GRID,
                          &Grid::PX_OFFSET, &Grid::PY_OFFSET, &Grid::PZ_OFFSET);

        find_node_offsets(gridpe, FNX_GRID, FNY_GRID, FNZ_GRID,
                          &Grid::FPX_OFFSET, &Grid::FPY_OFFSET, &Grid::FPZ_OFFSET);

        Grid::grid_first = 1;
    }

    void Grid::set_neighbors(int *list)
    {
        int idx;

        if(Grid::neighbor_first) return;

        for(idx = 0;idx < 6;idx++)
            Grid::neighbors[idx] = list[idx];

        Grid::neighbor_first = 1;

    }

    int Grid::get_PX0_GRID(void)
    {
	if(!Grid::grid_first)
	    rmg_error_handler("Grids not initialized. Please call set_grids first");
	return Grid::PX0_GRID;
    }
    int Grid::get_PY0_GRID(void)
    {
	if(!Grid::grid_first)
	    rmg_error_handler("Grids not initialized. Please call set_grids first");
	return Grid::PY0_GRID;
    }
    int Grid::get_PZ0_GRID(void)
    {
	if(!Grid::grid_first)
	    rmg_error_handler("Grids not initialized. Please call set_grids first");
	return Grid::PZ0_GRID;
    }
    int Grid::get_PX_OFFSET(void)
    {
	if(!Grid::grid_first)
	    rmg_error_handler("Grids not initialized. Please call set_grids first");
	return Grid::PX_OFFSET;
    }
    int Grid::get_PY_OFFSET(void)
    {
	if(!Grid::grid_first)
	    rmg_error_handler("Grids not initialized. Please call set_grids first");
	return Grid::PY_OFFSET;
    }
    int Grid::get_PZ_OFFSET(void)
    {
	if(!Grid::grid_first)
	    rmg_error_handler("Grids not initialized. Please call set_grids first");
	return Grid::PZ_OFFSET;
    }
    int Grid::get_FPX_OFFSET(void)
    {
	if(!Grid::grid_first)
	    rmg_error_handler("Grids not initialized. Please call set_grids first");
	return Grid::FPX_OFFSET;
    }
    int Grid::get_FPY_OFFSET(void)
    {
	if(!Grid::grid_first)
	    rmg_error_handler("Grids not initialized. Please call set_grids first");
	return Grid::FPY_OFFSET;
    }
    int Grid::get_FPZ_OFFSET(void)
    {
	if(!Grid::grid_first)
	    rmg_error_handler("Grids not initialized. Please call set_grids first");
	return Grid::FPZ_OFFSET;
    }
    int Grid::get_P0_BASIS(void)
    {
	if(!Grid::grid_first)
	    rmg_error_handler("Grids not initialized. Please call set_grids first");
	return Grid::P0_BASIS;
    }
    int Grid::get_FP0_BASIS(void)
    {
	if(!Grid::grid_first)
	    rmg_error_handler("Grids not initialized. Please call set_grids first");
	return Grid::FP0_BASIS;
    }
    int Grid::get_FPX0_GRID(void)
    {
	if(!Grid::grid_first)
	    rmg_error_handler("Grids not initialized. Please call set_grids first");
	return Grid::FPX0_GRID;
    }
    int Grid::get_FPY0_GRID(void)
    {
	if(!Grid::grid_first)
	    rmg_error_handler("Grids not initialized. Please call set_grids first");
	return Grid::FPY0_GRID;
    }
    int Grid::get_FPZ0_GRID(void)
    {
	if(!Grid::grid_first)
	    rmg_error_handler("Grids not initialized. Please call set_grids first");
	return Grid::FPZ0_GRID;
    }
    int Grid::get_ibrav_type(void)
    {
	if(!Grid::grid_first)
	    rmg_error_handler("Grids not initialized. Please call set_grids first");
	return ibrav;
    }
    void Grid::set_ibrav_type(int value)
    {
	Grid::ibrav = value;
    }
    void Grid::set_anisotropy(rmg_double_t a)
    {
        if(anisotropy_first) return;
        Grid::anisotropy = a;
        anisotropy_first = 1;
    }
    rmg_double_t Grid::get_anisotropy(void)
    {
	if(!Grid::anisotropy_first)
	    rmg_error_handler("Anisotropy not initialized. Please call set_anisotropy first");
	return Grid::anisotropy;
    }

    // Returns a pointer to the neighbors structure which contains the rank
    // of neighboring processors in three-dimensional space.
    int *Grid::get_neighbors(void)
    {
	if(!neighbor_first)
	    rmg_error_handler("Neighbor list not initialized. Please call set_neighbors first");

	return Grid::neighbors;
    }

/* Grid sizes on each PE */
int Grid::PX0_GRID;
int Grid::PY0_GRID;
int Grid::PZ0_GRID;

/* Grid offsets on each PE */
int Grid::PX_OFFSET;
int Grid::PY_OFFSET;
int Grid::PZ_OFFSET;

/* Basis size on each PE */
int Grid::P0_BASIS;

/* Fine grid sizes on each PE */
int Grid::FPX0_GRID;
int Grid::FPY0_GRID;
int Grid::FPZ0_GRID;

/* Fine Grid offsets on each PE */
int Grid::FPX_OFFSET;
int Grid::FPY_OFFSET;
int Grid::FPZ_OFFSET;

/* Fine grid basis size on each PE */
int Grid::FP0_BASIS;

/* Grid bravais lattice type */
int Grid::ibrav;

/* MPI stuff */
int Grid::gridpe;
int Grid::neighbors[6];

/* Grid anisotropy defined as the ratio of hmaxgrid to hmingrid. A value larger than 1.05 can lead to convergence problems. */
rmg_double_t Grid::anisotropy;

int Grid::neighbor_first=0;
int Grid::grid_first=0;
int Grid::anisotropy_first=0;


// C interfaces during transition
extern "C" void set_grids(int newgridpe, int ii, int jj, int kk)
{
  Grid G;
  G.set_grids(newgridpe, ii, jj, kk);
}
extern "C" void set_neighbors(int *newneighbors)
{
  Grid G;
  G.set_neighbors(newneighbors);
}
extern "C" void set_anisotropy(rmg_double_t newanisotropy)
{
  Grid G;
  G.set_anisotropy(newanisotropy);
}

extern "C" void set_ibrav_type(int newtype)
{
  Grid G;
  G.set_ibrav_type(newtype);
}
extern "C" int get_ibrav_type(int newtype)
{
  Grid G;
  return G.get_ibrav_type();
}

extern "C" int get_PX0_GRID(void)
{
  //Grid G;
  //return G.get_PX0_GRID();
  return Grid::PX0_GRID;
}
extern "C" int get_PY0_GRID(void)
{
  Grid G;
  return G.get_PY0_GRID();
}
extern "C" int get_PZ0_GRID(void)
{
  Grid G;
  return G.get_PZ0_GRID();
}
extern "C" int get_PX_OFFSET(void)
{
  Grid G;
  return G.get_PX_OFFSET();
}
extern "C" int get_PY_OFFSET(void)
{
  Grid G;
  return G.get_PY_OFFSET();
}
extern "C" int get_PZ_OFFSET(void)
{
  Grid G;
  return G.get_PZ_OFFSET();
}
extern "C" int get_FPX_OFFSET(void)
{
  Grid G;
  return G.get_FPX_OFFSET();
}
extern "C" int get_FPY_OFFSET(void)
{
  Grid G;
  return G.get_FPY_OFFSET();
}
extern "C" int get_FPZ_OFFSET(void)
{
  Grid G;
  return G.get_FPZ_OFFSET();
}
extern "C" int get_P0_BASIS(void)
{
  Grid G;
  return G.get_P0_BASIS();
}
extern "C" int get_FP0_BASIS(void)
{
  Grid G;
  return G.get_FP0_BASIS();
}
extern "C" int get_FPX0_GRID(void)
{
  Grid G;
  return G.get_FPX0_GRID();
}
extern "C" int get_FPY0_GRID(void)
{
  Grid G;
  return G.get_FPY0_GRID();
}
extern "C" int get_FPZ0_GRID(void)
{
  Grid G;
  return G.get_FPZ0_GRID();
}
extern "C" rmg_double_t get_anisotropy(void)
{
  Grid G;
  return G.get_anisotropy();
}
extern "C" int *get_neighbors(void)
{
  Grid G;
  return G.get_neighbors();
}
