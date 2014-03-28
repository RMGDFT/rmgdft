#include "BaseGrid.h"


/// C interface function
extern "C" int get_PE_X(void)
{
  BaseGrid G;
  return G.get_PE_X();
}
/// C interface function
extern "C" int get_PE_Y(void)
{
  BaseGrid G;
  return G.get_PE_Y();
}
/// C interface function
extern "C" int get_PE_Z(void)
{
  BaseGrid G;
  return G.get_PE_Z();
}
/// C interface function
extern "C" int get_NX_GRID(void)
{
  BaseGrid G;
  return G.get_NX_GRID(1);
}
/// C interface function
extern "C" int get_NY_GRID(void)
{
  BaseGrid G;
  return G.get_NY_GRID(1);
}
/// C interface function
extern "C" int get_NZ_GRID(void)
{
  BaseGrid G;
  return G.get_NZ_GRID(1);
}
/// C interface function
extern "C" double get_hxgrid(void)
{
    BaseGrid G;
    return G.get_hxgrid(1);
}
/// C interface function
extern "C" double get_hygrid(void)
{
    BaseGrid G;
    return G.get_hygrid(1);
}
/// C interface function
extern "C" double get_hzgrid(void)
{
    BaseGrid G;
    return G.get_hzgrid(1);
}
/// C interface function
extern "C" double get_hxxgrid(void)
{
    BaseGrid G;
    return G.get_hxgrid(G.default_FG_RATIO);
}
/// C interface function
extern "C" double get_hyygrid(void)
{
    BaseGrid G;
    return G.get_hygrid(G.default_FG_RATIO);
}
/// C interface function
extern "C" double get_hzzgrid(void)
{
    BaseGrid G;
    return G.get_hzgrid(G.default_FG_RATIO);
}
/// C interface function
extern "C" int get_FNX_GRID(void)
{
  BaseGrid G;
  return G.get_NX_GRID(G.default_FG_RATIO);
}
/// C interface function
extern "C" int get_FNY_GRID(void)
{
  BaseGrid G;
  return G.get_NY_GRID(G.default_FG_RATIO);
}
/// C interface function
extern "C" int get_FNZ_GRID(void)
{
  BaseGrid G;
  return G.get_NZ_GRID(G.default_FG_RATIO);
}
/// C interface function
extern "C" int get_FG_RATIO(void)
{
  BaseGrid G;
  return G.default_FG_RATIO;
}
/// C interface function
extern "C" void set_grids(int newNX_GRID, int newNY_GRID, int newNZ_GRID, int newPE_X, int newPE_Y, int newPE_Z, int newFG_RATIO)
{
  BaseGrid G;
  G.set_grids(newNX_GRID, newNY_GRID, newNZ_GRID, newPE_X, newPE_Y, newPE_Z, newFG_RATIO);
}
/// C interface function
extern "C" void set_nodes(int newgridpe)
{
  BaseGrid G;
  G.set_nodes(newgridpe);
}
/// C interface function
extern "C" void set_anisotropy(double newanisotropy)
{
  BaseGrid G;
  G.set_anisotropy(newanisotropy);
}
/// C interface function
extern "C" int get_PX0_GRID(void)
{
  BaseGrid G;
  return G.get_PX0_GRID(1);
}
/// C interface function
extern "C" int get_PY0_GRID(void)
{
  BaseGrid G;
  return G.get_PY0_GRID(1);
}
/// C interface function
extern "C" int get_PZ0_GRID(void)
{
  BaseGrid G;
  return G.get_PZ0_GRID(1);
}
/// C interface function
extern "C" int get_PX_OFFSET(void)
{
  BaseGrid G;
  return G.get_PX_OFFSET(1);
}
/// C interface function
extern "C" int get_PY_OFFSET(void)
{
  BaseGrid G;
  return G.get_PY_OFFSET(1);
}
/// C interface function
extern "C" int get_PZ_OFFSET(void)
{
  BaseGrid G;
  return G.get_PZ_OFFSET(1);
}
/// C interface function
extern "C" int get_FPX_OFFSET(void)
{
  BaseGrid G;
  int FPX_OFFSET, FPY_OFFSET, FPZ_OFFSET;
  G.find_node_offsets(G.get_gridpe(), G.get_NX_GRID(G.default_FG_RATIO), G.get_NY_GRID(G.default_FG_RATIO), G.get_NZ_GRID(G.default_FG_RATIO),
                          &FPX_OFFSET, &FPY_OFFSET, &FPZ_OFFSET);
  return FPX_OFFSET;
}
/// C interface function
extern "C" int get_FPY_OFFSET(void)
{
  BaseGrid G;
  int FPX_OFFSET, FPY_OFFSET, FPZ_OFFSET;
  G.find_node_offsets(G.get_gridpe(), G.get_NX_GRID(G.default_FG_RATIO), G.get_NY_GRID(G.default_FG_RATIO), G.get_NZ_GRID(G.default_FG_RATIO),
                          &FPX_OFFSET, &FPY_OFFSET, &FPZ_OFFSET);
  return FPY_OFFSET;
}
/// C interface function
extern "C" int get_FPZ_OFFSET(void)
{
  BaseGrid G;
  int FPX_OFFSET, FPY_OFFSET, FPZ_OFFSET;
  G.find_node_offsets(G.get_gridpe(), G.get_NX_GRID(G.default_FG_RATIO), G.get_NY_GRID(G.default_FG_RATIO), G.get_NZ_GRID(G.default_FG_RATIO),
                          &FPX_OFFSET, &FPY_OFFSET, &FPZ_OFFSET);
  return FPZ_OFFSET;
}
/// C interface function
extern "C" int get_P0_BASIS(void)
{
  BaseGrid G;
  return G.get_P0_BASIS(1);
}
/// C interface function
extern "C" int get_FP0_BASIS(void)
{
  BaseGrid G;
  return G.get_P0_BASIS(G.default_FG_RATIO);
}
/// C interface function
extern "C" int get_FPX0_GRID(void)
{
  BaseGrid G;
  return G.get_PX0_GRID(G.default_FG_RATIO);
}
/// C interface function
extern "C" int get_FPY0_GRID(void)
{
  BaseGrid G;
  return G.get_PY0_GRID(G.default_FG_RATIO);
}
/// C interface function
extern "C" int get_FPZ0_GRID(void)
{
  BaseGrid G;
  return G.get_PZ0_GRID(G.default_FG_RATIO);
}
/// C interface function
extern "C" double get_anisotropy(void)
{
  BaseGrid G;
  return G.get_anisotropy();
}
/// C interface function
extern "C" void pe2xyz(int pe, int *x, int *y, int *z)
{
  BaseGrid G;
  G.pe2xyz(pe, x, y, z);
}
/// C interface function
extern "C" int *get_neighbors(void)
{
  BaseGrid G;
  return G.get_neighbors();
}
/// C interface function
extern "C" void find_node_sizes(int gridpe, int nxgrid, int nygrid, int nzgrid, int *pxsize, int *pysize, int *pzsize)
{
  BaseGrid G;
  G.find_node_sizes(gridpe, nxgrid, nygrid, nzgrid, pxsize, pysize, pzsize);
}
/// C interface function
extern "C" void find_node_offsets(int gridpe, int nxgrid, int nygrid, int nzgrid, int *pxoffset, int *pyoffset, int *pzoffset)
{
  BaseGrid G;
  G.find_node_offsets(gridpe, nxgrid, nygrid, nzgrid, pxoffset, pyoffset, pzoffset);

}
