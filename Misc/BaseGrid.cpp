/*
    Used to access grid and node related data.

*/

#include "BaseGrid.h"


using namespace std;



    void BaseGrid::set_grids(int newNX_GRID, int newNY_GRID, int newNZ_GRID, int newPE_X, int newPE_Y, int newPE_Z, int newFG_NX, int newFG_NY, int newFG_NZ)
    {

        BaseGrid::NX_GRID = newNX_GRID;
        BaseGrid::NY_GRID = newNY_GRID;
        BaseGrid::NZ_GRID = newNZ_GRID;

        BaseGrid::PE_X = newPE_X;
        BaseGrid::PE_Y = newPE_Y;
        BaseGrid::PE_Z = newPE_Z;

        BaseGrid::FG_NX = newFG_NX;
        BaseGrid::FG_NY = newFG_NY;
        BaseGrid::FG_NZ = newFG_NZ;

        BaseGrid::FNX_GRID = newNX_GRID * newFG_NX;
        BaseGrid::FNY_GRID = newNY_GRID * newFG_NY;
        BaseGrid::FNZ_GRID = newNZ_GRID * newFG_NZ;

        BaseGrid::grid_first = 1;
    }

    void BaseGrid::set_nodes(int newgridpe, int ii, int jj, int kk)
    {
        int rem;

       if(!BaseGrid::grid_first)
            rmg_error_handler("Grids must be initialized before nodes. Please call set_grids first");

        BaseGrid::gridpe = newgridpe;

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
        BaseGrid::find_node_sizes(gridpe, BaseGrid::NX_GRID, BaseGrid::NY_GRID, BaseGrid::NZ_GRID, &BaseGrid::PX0_GRID, &BaseGrid::PY0_GRID, &BaseGrid::PZ0_GRID);
        BaseGrid::find_node_sizes(gridpe, BaseGrid::FNX_GRID, BaseGrid::FNY_GRID, BaseGrid::FNZ_GRID, &BaseGrid::FPX0_GRID, &BaseGrid::FPY0_GRID, &BaseGrid::FPZ0_GRID);

        BaseGrid::P0_BASIS = BaseGrid::PX0_GRID * BaseGrid::PY0_GRID * BaseGrid::PZ0_GRID;
        BaseGrid::FP0_BASIS = BaseGrid::FPX0_GRID * BaseGrid::FPY0_GRID * BaseGrid::FPZ0_GRID;

        // Now compute the global grid offset of the first point of the coarse and fine node grids
        BaseGrid::find_node_offsets(gridpe, BaseGrid::NX_GRID, BaseGrid::NY_GRID, BaseGrid::NZ_GRID,
                          &BaseGrid::PX_OFFSET, &BaseGrid::PY_OFFSET, &BaseGrid::PZ_OFFSET);

        BaseGrid::find_node_offsets(gridpe, BaseGrid::FNX_GRID, BaseGrid::FNY_GRID, BaseGrid::FNZ_GRID,
                          &BaseGrid::FPX_OFFSET, &BaseGrid::FPY_OFFSET, &BaseGrid::FPZ_OFFSET);


    }

    void BaseGrid::set_neighbors(int *list)
    {
        int idx;

        if(BaseGrid::neighbor_first) return;

        for(idx = 0;idx < 6;idx++)
            BaseGrid::neighbors[idx] = list[idx];

        BaseGrid::neighbor_first = 1;

    }

    int BaseGrid::find_node_sizes(int gridpe, int nxgrid, int nygrid, int nzgrid, int *pxsize, int *pysize, int *pzsize)
    {

	int ii, jj, kk;
	int idx, ix, iy, iz, mfac;

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

    int BaseGrid::find_node_offsets(int gridpe, int nxgrid, int nygrid, int nzgrid, int *pxoffset, int *pyoffset, int *pzoffset)
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

    int BaseGrid::get_PX0_GRID(void)
    {
	if(!BaseGrid::grid_first)
	    rmg_error_handler("Grids not initialized. Please call set_grids first");
	return BaseGrid::PX0_GRID;
    }
    int BaseGrid::get_PY0_GRID(void)
    {
	if(!BaseGrid::grid_first)
	    rmg_error_handler("Grids not initialized. Please call set_grids first");
	return BaseGrid::PY0_GRID;
    }
    int BaseGrid::get_PZ0_GRID(void)
    {
	if(!BaseGrid::grid_first)
	    rmg_error_handler("Grids not initialized. Please call set_grids first");
	return BaseGrid::PZ0_GRID;
    }
    int BaseGrid::get_PX_OFFSET(void)
    {
	if(!BaseGrid::grid_first)
	    rmg_error_handler("Grids not initialized. Please call set_grids first");
	return BaseGrid::PX_OFFSET;
    }
    int BaseGrid::get_PY_OFFSET(void)
    {
	if(!BaseGrid::grid_first)
	    rmg_error_handler("Grids not initialized. Please call set_grids first");
	return BaseGrid::PY_OFFSET;
    }
    int BaseGrid::get_PZ_OFFSET(void)
    {
	if(!BaseGrid::grid_first)
	    rmg_error_handler("Grids not initialized. Please call set_grids first");
	return BaseGrid::PZ_OFFSET;
    }
    int BaseGrid::get_FPX_OFFSET(void)
    {
	if(!BaseGrid::grid_first)
	    rmg_error_handler("Grids not initialized. Please call set_grids first");
	return BaseGrid::FPX_OFFSET;
    }
    int BaseGrid::get_FPY_OFFSET(void)
    {
	if(!BaseGrid::grid_first)
	    rmg_error_handler("Grids not initialized. Please call set_grids first");
	return BaseGrid::FPY_OFFSET;
    }
    int BaseGrid::get_FPZ_OFFSET(void)
    {
	if(!BaseGrid::grid_first)
	    rmg_error_handler("Grids not initialized. Please call set_grids first");
	return BaseGrid::FPZ_OFFSET;
    }
    int BaseGrid::get_P0_BASIS(void)
    {
	if(!BaseGrid::grid_first)
	    rmg_error_handler("Grids not initialized. Please call set_grids first");
	return BaseGrid::P0_BASIS;
    }
    int BaseGrid::get_FP0_BASIS(void)
    {
	if(!BaseGrid::grid_first)
	    rmg_error_handler("Grids not initialized. Please call set_grids first");
	return BaseGrid::FP0_BASIS;
    }
    int BaseGrid::get_FPX0_GRID(void)
    {
	if(!BaseGrid::grid_first)
	    rmg_error_handler("Grids not initialized. Please call set_grids first");
	return BaseGrid::FPX0_GRID;
    }
    int BaseGrid::get_FPY0_GRID(void)
    {
	if(!BaseGrid::grid_first)
	    rmg_error_handler("Grids not initialized. Please call set_grids first");
	return BaseGrid::FPY0_GRID;
    }
    int BaseGrid::get_FPZ0_GRID(void)
    {
	if(!BaseGrid::grid_first)
	    rmg_error_handler("Grids not initialized. Please call set_grids first");
	return BaseGrid::FPZ0_GRID;
    }
    int BaseGrid::get_ibrav_type(void)
    {
	if(!BaseGrid::ibrav_first)
	    rmg_error_handler("Lattice type not initialized. Please call set_ibrav first");
	return ibrav;
    }
    void BaseGrid::set_ibrav_type(int value)
    {
	BaseGrid::ibrav = value;
        BaseGrid::ibrav_first = 1;
    }
    void BaseGrid::set_anisotropy(rmg_double_t a)
    {
        if(BaseGrid::anisotropy_first) return;
        BaseGrid::anisotropy = a;
        anisotropy_first = 1;
    }
    rmg_double_t BaseGrid::get_anisotropy(void)
    {
	if(!BaseGrid::anisotropy_first)
	    rmg_error_handler("Anisotropy not initialized. Please call set_anisotropy first");
	return BaseGrid::anisotropy;
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

    // Returns a pointer to the neighbors structure which contains the rank
    // of neighboring processors in three-dimensional space.
    int *BaseGrid::get_neighbors(void)
    {
	if(!BaseGrid::neighbor_first)
	    rmg_error_handler("Neighbor list not initialized. Please call set_neighbors first");

	return BaseGrid::neighbors;
    }


/* Global coarse grid dimensions */
int BaseGrid::NX_GRID;
int BaseGrid::NY_GRID;
int BaseGrid::NZ_GRID;

/* Fine grid/coarse grid ratio */
int BaseGrid::FG_NX;
int BaseGrid::FG_NY;
int BaseGrid::FG_NZ;

/* Global fine grid dimensions */
int BaseGrid::FNX_GRID;
int BaseGrid::FNY_GRID;
int BaseGrid::FNZ_GRID;

/* Node (PE) dimensions */
int BaseGrid::PE_X;
int BaseGrid::PE_Y;
int BaseGrid::PE_Z;

/* Grid sizes on each PE */
int BaseGrid::PX0_GRID;
int BaseGrid::PY0_GRID;
int BaseGrid::PZ0_GRID;

/* Grid offsets on each PE */
int BaseGrid::PX_OFFSET;
int BaseGrid::PY_OFFSET;
int BaseGrid::PZ_OFFSET;

/* Basis size on each PE */
int BaseGrid::P0_BASIS;

/* Fine grid sizes on each PE */
int BaseGrid::FPX0_GRID;
int BaseGrid::FPY0_GRID;
int BaseGrid::FPZ0_GRID;

/* Fine Grid offsets on each PE */
int BaseGrid::FPX_OFFSET;
int BaseGrid::FPY_OFFSET;
int BaseGrid::FPZ_OFFSET;

/* Fine grid basis size on each PE */
int BaseGrid::FP0_BASIS;

/* Grid bravais lattice type */
int BaseGrid::ibrav;

/* MPI stuff */
int BaseGrid::gridpe;
int BaseGrid::neighbors[6];

/* Grid anisotropy defined as the ratio of hmaxgrid to hmingrid. A value larger than 1.05 can lead to convergence problems. */
rmg_double_t BaseGrid::anisotropy;

int BaseGrid::ibrav_first=0;
int BaseGrid::neighbor_first=0;
int BaseGrid::grid_first=0;
int BaseGrid::anisotropy_first=0;



// C interfaces for use during transition
extern "C" int get_PE_X(void)
{
  BaseGrid G;
  return G.PE_X;
}
extern "C" int get_PE_Y(void)
{
  BaseGrid G;
  return G.PE_Y;
}
extern "C" int get_PE_Z(void)
{
  BaseGrid G;
  return G.PE_Z;
}
extern "C" int get_NX_GRID(void)
{
  BaseGrid G;
  return G.NX_GRID;
}
extern "C" int get_NY_GRID(void)
{
  BaseGrid G;
  return G.NY_GRID;
}
extern "C" int get_NZ_GRID(void)
{
  BaseGrid G;
  return G.NZ_GRID;
}
extern "C" int get_FNX_GRID(void)
{
  BaseGrid G;
  return G.FNX_GRID;
}
extern "C" int get_FNY_GRID(void)
{
  BaseGrid G;
  return G.FNY_GRID;
}
extern "C" int get_FNZ_GRID(void)
{
  BaseGrid G;
  return G.FNZ_GRID;
}
extern "C" int get_FG_NX(void)
{
  BaseGrid G;
  return G.FG_NX;
}
extern "C" int get_FG_NY(void)
{
  BaseGrid G;
  return G.FG_NY;
}
extern "C" int get_FG_NZ(void)
{
  BaseGrid G;
  return G.FG_NZ;
}
extern "C" void set_grids(int newNX_GRID, int newNY_GRID, int newNZ_GRID, int newPE_X, int newPE_Y, int newPE_Z, int newFG_NX, int newFG_NY, int newFG_NZ)
{
  BaseGrid G;
  G.set_grids(newNX_GRID, newNY_GRID, newNZ_GRID, newPE_X, newPE_Y, newPE_Z, newFG_NX, newFG_NY, newFG_NZ);
}
extern "C" void set_nodes(int newgridpe, int ii, int jj, int kk)
{
  BaseGrid G;
  G.set_nodes(newgridpe, ii, jj, kk);
}
extern "C" void set_neighbors(int *newneighbors)
{
  BaseGrid G;
  G.set_neighbors(newneighbors);
}
extern "C" void set_anisotropy(rmg_double_t newanisotropy)
{
  BaseGrid G;
  G.set_anisotropy(newanisotropy);
}

extern "C" void set_ibrav_type(int newtype)
{
  BaseGrid G;
  G.set_ibrav_type(newtype);
}
extern "C" int get_ibrav_type(int newtype)
{
  BaseGrid G;
  return G.get_ibrav_type();
}

extern "C" int get_PX0_GRID(void)
{
  BaseGrid G;
  return G.get_PX0_GRID();
}
extern "C" int get_PY0_GRID(void)
{
  BaseGrid G;
  return G.get_PY0_GRID();
}
extern "C" int get_PZ0_GRID(void)
{
  BaseGrid G;
  return G.get_PZ0_GRID();
}
extern "C" int get_PX_OFFSET(void)
{
  BaseGrid G;
  return G.get_PX_OFFSET();
}
extern "C" int get_PY_OFFSET(void)
{
  BaseGrid G;
  return G.get_PY_OFFSET();
}
extern "C" int get_PZ_OFFSET(void)
{
  BaseGrid G;
  return G.get_PZ_OFFSET();
}
extern "C" int get_FPX_OFFSET(void)
{
  BaseGrid G;
  return G.get_FPX_OFFSET();
}
extern "C" int get_FPY_OFFSET(void)
{
  BaseGrid G;
  return G.get_FPY_OFFSET();
}
extern "C" int get_FPZ_OFFSET(void)
{
  BaseGrid G;
  return G.get_FPZ_OFFSET();
}
extern "C" int get_P0_BASIS(void)
{
  BaseGrid G;
  return G.get_P0_BASIS();
}
extern "C" int get_FP0_BASIS(void)
{
  BaseGrid G;
  return G.get_FP0_BASIS();
}
extern "C" int get_FPX0_GRID(void)
{
  BaseGrid G;
  return G.get_FPX0_GRID();
}
extern "C" int get_FPY0_GRID(void)
{
  BaseGrid G;
  return G.get_FPY0_GRID();
}
extern "C" int get_FPZ0_GRID(void)
{
  BaseGrid G;
  return G.get_FPZ0_GRID();
}
extern "C" rmg_double_t get_anisotropy(void)
{
  BaseGrid G;
  return G.get_anisotropy();
}
extern "C" void pe2xyz(int pe, int *x, int *y, int *z)
{
  BaseGrid G;
  G.pe2xyz(pe, x, y, z);
}
extern "C" int *get_neighbors(void)
{
  BaseGrid G;
  return G.get_neighbors();
}
extern "C" int find_node_sizes(int gridpe, int nxgrid, int nygrid, int nzgrid, int *pxsize, int *pysize, int *pzsize)
{
  BaseGrid G;
  return G.find_node_sizes(gridpe, nxgrid, nygrid, nzgrid, pxsize, pysize, pzsize);
}
extern "C" int find_node_offsets(int gridpe, int nxgrid, int nygrid, int nzgrid, int *pxoffset, int *pyoffset, int *pzoffset)
{
  BaseGrid G;
  return G.find_node_offsets(gridpe, nxgrid, nygrid, nzgrid, pxoffset, pyoffset, pzoffset);

}
