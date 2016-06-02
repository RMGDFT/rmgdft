/*
 *
 * Copyright 2014 The RMG Project Developers. See the COPYRIGHT file 
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

#include "transition.h"
#include "packfuncs.h"

/// C interface function
extern "C" int get_PE_X(void)
{
  return Rmg_G->get_PE_X();
}
/// C interface function
extern "C" int get_PE_Y(void)
{
  return Rmg_G->get_PE_Y();
}
/// C interface function
extern "C" int get_PE_Z(void)
{
  return Rmg_G->get_PE_Z();
}
/// C interface function
extern "C" int get_NX_GRID(void)
{
  return Rmg_G->get_NX_GRID(1);
}
/// C interface function
extern "C" int get_NY_GRID(void)
{
  return Rmg_G->get_NY_GRID(1);
}
/// C interface function
extern "C" int get_NZ_GRID(void)
{
  return Rmg_G->get_NZ_GRID(1);
}
/// C interface function
extern "C" double get_hxgrid(void)
{
    return Rmg_G->get_hxgrid(1);
}
/// C interface function
extern "C" double get_hygrid(void)
{
    return Rmg_G->get_hygrid(1);
}
/// C interface function
extern "C" double get_hzgrid(void)
{
    return Rmg_G->get_hzgrid(1);
}
/// C interface function
extern "C" double get_hxxgrid(void)
{
    return Rmg_G->get_hxgrid(Rmg_G->default_FG_RATIO);
}
/// C interface function
extern "C" double get_hyygrid(void)
{
    return Rmg_G->get_hygrid(Rmg_G->default_FG_RATIO);
}
/// C interface function
extern "C" double get_hzzgrid(void)
{
    return Rmg_G->get_hzgrid(Rmg_G->default_FG_RATIO);
}
/// C interface function
extern "C" int get_FNX_GRID(void)
{
  return Rmg_G->get_NX_GRID(Rmg_G->default_FG_RATIO);
}
/// C interface function
extern "C" int get_FNY_GRID(void)
{
  return Rmg_G->get_NY_GRID(Rmg_G->default_FG_RATIO);
}
/// C interface function
extern "C" int get_FNZ_GRID(void)
{
  return Rmg_G->get_NZ_GRID(Rmg_G->default_FG_RATIO);
}
/// C interface function
extern "C" int get_FG_RATIO(void)
{
  return Rmg_G->default_FG_RATIO;
}
/// C interface function
extern "C" void set_grids(int newNX_GRID, int newNY_GRID, int newNZ_GRID, int newPE_X, int newPE_Y, int newPE_Z, int newFG_RATIO)
{
  Rmg_G = new BaseGrid(newNX_GRID, newNY_GRID, newNZ_GRID, newPE_X, newPE_Y, newPE_Z, 0, newFG_RATIO);
}
/// C interface function
extern "C" void set_rank(int newrank, MPI_Comm comm)
{
  Rmg_G->set_rank(newrank, comm);
}
/// C interface function
extern "C" int get_PX0_GRID(void)
{
  return Rmg_G->get_PX0_GRID(1);
}
/// C interface function
extern "C" int get_PY0_GRID(void)
{
  return Rmg_G->get_PY0_GRID(1);
}
/// C interface function
extern "C" int get_PZ0_GRID(void)
{
  return Rmg_G->get_PZ0_GRID(1);
}
/// C interface function
extern "C" int get_PX_OFFSET(void)
{
  return Rmg_G->get_PX_OFFSET(1);
}
/// C interface function
extern "C" int get_PY_OFFSET(void)
{
  return Rmg_G->get_PY_OFFSET(1);
}
/// C interface function
extern "C" int get_PZ_OFFSET(void)
{
  return Rmg_G->get_PZ_OFFSET(1);
}
/// C interface function
extern "C" int get_FPX_OFFSET(void)
{
  int FPX_OFFSET, FPY_OFFSET, FPZ_OFFSET;
  Rmg_G->find_node_offsets(Rmg_G->get_rank(), Rmg_G->get_NX_GRID(Rmg_G->default_FG_RATIO), Rmg_G->get_NY_GRID(Rmg_G->default_FG_RATIO), Rmg_G->get_NZ_GRID(Rmg_G->default_FG_RATIO),
                          &FPX_OFFSET, &FPY_OFFSET, &FPZ_OFFSET);
  return FPX_OFFSET;
}
/// C interface function
extern "C" int get_FPY_OFFSET(void)
{
  int FPX_OFFSET, FPY_OFFSET, FPZ_OFFSET;
  Rmg_G->find_node_offsets(Rmg_G->get_rank(), Rmg_G->get_NX_GRID(Rmg_G->default_FG_RATIO), Rmg_G->get_NY_GRID(Rmg_G->default_FG_RATIO), Rmg_G->get_NZ_GRID(Rmg_G->default_FG_RATIO),
                          &FPX_OFFSET, &FPY_OFFSET, &FPZ_OFFSET);
  return FPY_OFFSET;
}
/// C interface function
extern "C" int get_FPZ_OFFSET(void)
{
  int FPX_OFFSET, FPY_OFFSET, FPZ_OFFSET;
  Rmg_G->find_node_offsets(Rmg_G->get_rank(), Rmg_G->get_NX_GRID(Rmg_G->default_FG_RATIO), Rmg_G->get_NY_GRID(Rmg_G->default_FG_RATIO), Rmg_G->get_NZ_GRID(Rmg_G->default_FG_RATIO),
                          &FPX_OFFSET, &FPY_OFFSET, &FPZ_OFFSET);
  return FPZ_OFFSET;
}
/// C interface function
extern "C" int get_P0_BASIS(void)
{
  return Rmg_G->get_P0_BASIS(1);
}
/// C interface function
extern "C" int get_FP0_BASIS(void)
{
  return Rmg_G->get_P0_BASIS(Rmg_G->default_FG_RATIO);
}
/// C interface function
extern "C" int get_FPX0_GRID(void)
{
  return Rmg_G->get_PX0_GRID(Rmg_G->default_FG_RATIO);
}
/// C interface function
extern "C" int get_FPY0_GRID(void)
{
  return Rmg_G->get_PY0_GRID(Rmg_G->default_FG_RATIO);
}
/// C interface function
extern "C" int get_FPZ0_GRID(void)
{
  return Rmg_G->get_PZ0_GRID(Rmg_G->default_FG_RATIO);
}
/// C interface function
extern "C" void pe2xyz(int pe, int *x, int *y, int *z)
{
  Rmg_G->pe2xyz(pe, x, y, z);
}
/// C interface function
extern "C" int *get_neighbors(void)
{
  return Rmg_G->get_neighbors();
}
/// C interface function
extern "C" void find_node_sizes(int gridpe, int nxgrid, int nygrid, int nzgrid, int *pxsize, int *pysize, int *pzsize)
{
  Rmg_G->find_node_sizes(gridpe, nxgrid, nygrid, nzgrid, pxsize, pysize, pzsize);
}
/// C interface function
extern "C" void find_node_offsets(int gridpe, int nxgrid, int nygrid, int nzgrid, int *pxoffset, int *pyoffset, int *pzoffset)
{
  Rmg_G->find_node_offsets(gridpe, nxgrid, nygrid, nzgrid, pxoffset, pyoffset, pzoffset);

}
extern "C" void pack_vhdtos (double * s, double * d, int dimx, int dimy, int dimz, int boundaryflag)
{
    CPP_pack_dtos(Rmg_G, s, d, dimx, dimy, dimz, boundaryflag);
}
extern "C" void pack_vhstod (double * s, double * d, int dimx, int dimy, int dimz, int boundaryflag)
{
    CPP_pack_stod (Rmg_G, s, d, dimx, dimy, dimz, boundaryflag);
}

