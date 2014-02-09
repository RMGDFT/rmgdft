/*
    Used to access grid and node related data.

*/

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



    void Grid::set_grids(int newNX_GRID, int newNY_GRID, int newNZ_GRID, int newPE_X, int newPE_Y, int newPE_Z, int newFG_NX, int newFG_NY, int newFG_NZ)
    {

        Grid::NX_GRID = newNX_GRID;
        Grid::NY_GRID = newNY_GRID;
        Grid::NZ_GRID = newNZ_GRID;

        Grid::PE_X = newPE_X;
        Grid::PE_Y = newPE_Y;
        Grid::PE_Z = newPE_Z;

        Grid::FG_NX = newFG_NX;
        Grid::FG_NY = newFG_NY;
        Grid::FG_NZ = newFG_NZ;

        Grid::FNX_GRID = newNX_GRID * newFG_NX;
        Grid::FNY_GRID = newNY_GRID * newFG_NY;
        Grid::FNZ_GRID = newNZ_GRID * newFG_NZ;

        Grid::grid_first = 1;
    }

    void Grid::set_nodes(int newgridpe, int ii, int jj, int kk)
    {
        int rem;

       if(!Grid::grid_first)
            rmg_error_handler("Grids must be initialized before nodes. Please call set_grids first");

        Grid::gridpe = newgridpe;

        // Compute grid sizes for each node.
        Grid::PX0_GRID = Grid::NX_GRID / Grid::PE_X;
        rem = Grid::NX_GRID % Grid::PE_X;
        if(rem && (ii < rem)) Grid::PX0_GRID++;

        Grid::PY0_GRID = Grid::NY_GRID / Grid::PE_Y;
        rem = Grid::NY_GRID % Grid::PE_Y;
        if(rem && (jj < rem)) Grid::PY0_GRID++;

        Grid::PZ0_GRID = Grid::NZ_GRID / Grid::PE_Z;
        rem = Grid::NZ_GRID % Grid::PE_Z;
        if(rem && (kk < rem)) Grid::PZ0_GRID++;

        // Adjust if needed
        find_node_sizes(gridpe, Grid::NX_GRID, Grid::NY_GRID, Grid::NZ_GRID, &Grid::PX0_GRID, &Grid::PY0_GRID, &Grid::PZ0_GRID);
        find_node_sizes(gridpe, Grid::FNX_GRID, Grid::FNY_GRID, Grid::FNZ_GRID, &Grid::FPX0_GRID, &Grid::FPY0_GRID, &Grid::FPZ0_GRID);

        Grid::P0_BASIS = Grid::PX0_GRID * Grid::PY0_GRID * Grid::PZ0_GRID;
        Grid::FP0_BASIS = Grid::FPX0_GRID * Grid::FPY0_GRID * Grid::FPZ0_GRID;

        // Now compute the global grid offset of the first point of the coarse and fine node grids
        find_node_offsets(gridpe, Grid::NX_GRID, Grid::NY_GRID, Grid::NZ_GRID,
                          &Grid::PX_OFFSET, &Grid::PY_OFFSET, &Grid::PZ_OFFSET);

        find_node_offsets(gridpe, Grid::FNX_GRID, Grid::FNY_GRID, Grid::FNZ_GRID,
                          &Grid::FPX_OFFSET, &Grid::FPY_OFFSET, &Grid::FPZ_OFFSET);


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
	if(!Grid::ibrav_first)
	    rmg_error_handler("Lattice type not initialized. Please call set_ibrav first");
	return ibrav;
    }
    void Grid::set_ibrav_type(int value)
    {
	Grid::ibrav = value;
        Grid::ibrav_first = 1;
    }
    void Grid::set_anisotropy(rmg_double_t a)
    {
        if(Grid::anisotropy_first) return;
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
	if(!Grid::neighbor_first)
	    rmg_error_handler("Neighbor list not initialized. Please call set_neighbors first");

	return Grid::neighbors;
    }


/* Global coarse grid dimensions */
int Grid::NX_GRID;
int Grid::NY_GRID;
int Grid::NZ_GRID;

/* Fine grid/coarse grid ratio */
int Grid::FG_NX;
int Grid::FG_NY;
int Grid::FG_NZ;

/* Global fine grid dimensions */
int Grid::FNX_GRID;
int Grid::FNY_GRID;
int Grid::FNZ_GRID;

/* Node (PE) dimensions */
int Grid::PE_X;
int Grid::PE_Y;
int Grid::PE_Z;

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

int Grid::ibrav_first=0;
int Grid::neighbor_first=0;
int Grid::grid_first=0;
int Grid::anisotropy_first=0;



// C interfaces for use during transition
extern "C" int get_PE_X(void)
{
  Grid G;
  return G.PE_X;
}
extern "C" int get_PE_Y(void)
{
  Grid G;
  return G.PE_Y;
}
extern "C" int get_PE_Z(void)
{
  Grid G;
  return G.PE_Z;
}
extern "C" int get_NX_GRID(void)
{
  Grid G;
  return G.NX_GRID;
}
extern "C" int get_NY_GRID(void)
{
  Grid G;
  return G.NY_GRID;
}
extern "C" int get_NZ_GRID(void)
{
  Grid G;
  return G.NZ_GRID;
}
extern "C" int get_FNX_GRID(void)
{
  Grid G;
  return G.FNX_GRID;
}
extern "C" int get_FNY_GRID(void)
{
  Grid G;
  return G.FNY_GRID;
}
extern "C" int get_FNZ_GRID(void)
{
  Grid G;
  return G.FNZ_GRID;
}
extern "C" int get_FG_NX(void)
{
  Grid G;
  return G.FG_NX;
}
extern "C" int get_FG_NY(void)
{
  Grid G;
  return G.FG_NY;
}
extern "C" int get_FG_NZ(void)
{
  Grid G;
  return G.FG_NZ;
}
extern "C" void set_grids(int newNX_GRID, int newNY_GRID, int newNZ_GRID, int newPE_X, int newPE_Y, int newPE_Z, int newFG_NX, int newFG_NY, int newFG_NZ)
{
  Grid G;
  G.set_grids(newNX_GRID, newNY_GRID, newNZ_GRID, newPE_X, newPE_Y, newPE_Z, newFG_NX, newFG_NY, newFG_NZ);
}
extern "C" void set_nodes(int newgridpe, int ii, int jj, int kk)
{
  Grid G;
  G.set_nodes(newgridpe, ii, jj, kk);
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
  Grid G;
  return G.get_PX0_GRID();
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
